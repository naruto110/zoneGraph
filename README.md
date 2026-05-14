# zoneGraph 项目说明

本文档用于快速恢复项目上下文。当前工程围绕动态图邻接表在 ZNS SSD 上的存储、更新、GC 和恢复进行实验。代码里同时保留了原始 payload-only 路径和新增实验路径，后续修改时应优先保持这些路径可对比。

## 1. 总体思路

每个顶点维护一个邻接表。小度顶点优先放在 DRAM 的 `cache_buf`，超过阈值后迁移到 PMA。PMA 的物理块可以暂存在 DRAM，也可以刷写到 ZNS 上。

核心数据流：

```text
dataset edge
  -> graph_insert_directed_edge()
  -> cache_buf 或 PMA block
  -> dumpToZone*()
  -> ZNS physical block / ZNS record
  -> graph_recovery() 或 recover_graph_from_zns_records()
```

默认旧路径：

```text
[payload: 64 * int32_t]
pm_blk.global_wp = payload 起始地址
分配步长 = ZONE_BLK_SIZE
```

ZNS-native record 路径：

```text
[ZNSBlockHeader][payload: 64 * int32_t][padding 到 512B]
pm_blk.global_wp = record 起始地址
分配步长 = ZNS_RECORD_SIZE
```

不加 `--zns-record` 时，仍走旧路径。加 `--zns-record` 后才使用带 header 的 record 格式。

record 模式下会预留三个特殊 zone：

```text
0 .. global_num_queues - 1: data zone，用于 DATA record
global_empty_zone:         GC 搬迁使用的空白 zone
global_redo_zone:          redo log zone，只追加边更新日志
global_meta_zone:          metadata zone，只追加 COMMIT record
```

当前约定：

```cpp
global_meta_zone = zns.nr_zones - 1;
global_redo_zone = zns.nr_zones - 2;
global_empty_zone = zns.nr_zones - 3;
global_num_queues = zns.nr_zones - 3;
```

因此 `--zns-record` 至少需要 4 个 zone：data、GC empty、redo log、metadata。

### 1.1 ZNS zone block 划分与 PMA 映射

项目显式把每个 ZNS zone 划分成固定大小的可写 block：

```text
payload-only 模式:
  block size = ZONE_BLK_SIZE = PM_BLK_SIZE * sizeof(int32_t)

record 模式:
  block/record size = ZNS_RECORD_SIZE = 512B
```

每个顶点的 PMA 邻接表被切分为多个逻辑 PMA block：

```text
Vertex::blk_list[pma_blk_idx] -> pm_blk
pm_blk.global_wp             -> ZNS physical offset
```

也就是说，系统维护的是：

```text
(src_vid, pma_blk_idx) -> ZNS zone block / ZNS record
```

更新时不原地覆盖旧 block，而是把受影响的 PMA block 作为新版本追加写到 ZNS，并更新 `pm_blk.global_wp`。旧版本通过 `global_zone_vec` 中的 zone/block metadata 失效，后续由 GC 搬迁或回收。record 模式下，`ZNSBlockHeader` 会持久化 `src_vid` 和 `pma_blk_idx`，因此恢复时可以扫描 ZNS records 重建该映射，而不依赖掉电前的内存状态。

## 2. 主要文件

`main.cpp`

程序入口、命令行参数解析、设备打开、zone reset、构图模式选择、额外插入/删除、打印、恢复入口。

`global.h`

全局常量、核心结构体和全局变量。重要结构包括：

- `Vertex`：顶点状态，包含 `degree`、`cache_buf`、`blk_list`、`pma`、`update_count`、`vclass`。
- `pm_blk`：顶点 PMA 的一个逻辑块描述，包含 `global_wp`、可选 DRAM `buf`，以及内存态 block meta。
- `ZoneStatus` / `blk_info`：运行时 zone/block 到顶点逻辑块的映射。
- `ZNSBlockHeader`：ZNS-native recovery 使用的持久化 record header。
- `VertexClass`：`cold / hot / hub`。
- `GCMode`：三种 GC 模式。

`utils.h`

ZNS 设备打开、zone 信息、块分配队列、底层读写、zone metadata 标记、旧 GC 和 record 读写函数。关键函数：

- `zns_open()`
- `generate_block_baseline()`
- `generate_block_opt()`
- `requestMemoryBlocks()`
- `writeBlockToZone()` / `readBlockFromZone()`
- `writeZNSRecord()` / `readZNSRecordPayload()`
- `beginZNSEpoch()` / `commitZNSEpoch()`
- `readZNSRecordChecked()`
- `mark_zone()` / `invalidBlock()` / `test_block()`
- `cacheBlock()`

`pma.h`

PMA 查找、插入、删除、扩容、刷盘、恢复和实验策略。关键函数：

- `binary_pma_find()`
- `pma_insert()`
- `pma_delete()`
- `pma_compact()`
- `graph_insert_directed_edge()`
- `graph_delete_directed_edge()`
- `dumpToZone()`
- `dumpToZoneParallel()`
- `dumpToZoneZoneSerial()`
- `garbage_collection_by_mode()`
- `persist_cache_vertices_as_records()`
- `recover_graph_from_zns_records()`

`benchmark.h`

图计算 benchmark，包括 PageRank、BFS、CC 等。读取邻接时已按 `blk_list.empty()` 判断 cache/PMA 路径，避免小度但已迁移 PMA 的顶点被错误当成 cache-only 顶点。

`Makefile`

构建 `main_pma_check_new`。

## 3. 邻居存储和更新逻辑

每个顶点的邻居有两种状态。

### 3.1 cache_buf

当顶点度数小于 `CACHE_BLK_SIZE` 时，邻居直接放在：

```cpp
Vertex::cache_buf[CACHE_BLK_SIZE]
```

当前 demo 阈值为：

```cpp
#define CACHE_BLK_SIZE 4
```

这是 `cache_buf` 的数组上限。运行时可以用参数调整实际缓存阈值：

```bash
--cache-threshold <0..4>
--no-cache
```

`--no-cache` 等价于 `--cache-threshold 0`，顶点第一次插入就迁入 PMA 路径，适合作为无小度缓存的恢复/写放大对比实验。默认阈值仍为 4。恢复实验中，保守 redo 回放可以补齐 cache-only 顶点；如果关闭 cache，则 `committed_op_id` 的全图高水位语义更容易成立。

### 3.2 PMA block

当 `cache_buf` 满后，已有 cache 邻居整体迁移到 PMA：

```cpp
Vertex::blk_list
Vertex::pma
```

PMA 逻辑上仍以 `PM_BLK_SIZE=64` 个 `int32_t` 为一个块。空位用 `-1` 表示。

每个 `pm_blk` 额外维护一组内存态 block meta：

```cpp
int32_t min_val;
int32_t max_val;
int32_t valid_count;
uint64_t valid_bitmap;
```

含义：

```text
min_val / max_val:
  当前 block 内有效邻居值范围。当前已维护，后续可用于 block-level 查询剪枝。

valid_count:
  当前 block 内非 -1 slot 数。为 0 时查询可以直接跳过整个 block。

valid_bitmap:
  bit i 表示 block 内第 i 个 slot 是否有效。由于 PM_BLK_SIZE=64，一个 uint64_t 正好覆盖一个 block。
```

这些 meta 不持久化到 ZNS，也不进入 `ZNSBlockHeader`。它们是由 block payload 派生出的内存索引，恢复时会重新读取 block 内容并重建。

插入流程：

```text
graph_insert_directed_edge()
  -> cache_buf 未满：cache_insert_unique()
  -> cache_buf 满：创建第一个 PMA block，把 cache 邻居迁入 PMA
  -> PMA 顶点：binary_pma_find() + pma_insert()
```

删除流程：

```text
graph_delete_directed_edge()
  -> cache-only：cache_delete()
  -> PMA：pma_delete()
  -> 可选 --compact-delete：pma_compact()
```

删除边后，如果 dirty block 数达到 `cache_blk_threshold`，由 main 根据当前模式刷回 ZNS。

### 3.3 PMA 查询与 block meta

PMA 查询入口仍是：

```cpp
binary_pma_find()
binary_pma_find_recovery()
```

当前查询优化集中在 `find_prev_valid_index()`：

```text
1. 根据 index 定位 PMA block。
2. 如果 block.valid_count == 0，直接跳到前一个 block，不读取 ZNS。
3. 如果 block 有有效元素但 buf == nullptr，读取整个 block 到 DRAM，并刷新 block meta。
4. 使用 valid_bitmap 在 block 内找到不大于当前 index 的最近有效 slot。
```

这替代了旧实现中跨 slot 向左扫描 `-1` 空洞的逻辑。当前实际用于查询加速的是 `valid_count` 和 `valid_bitmap`；`min_val / max_val` 已同步维护，但还没有接入 block-level 二分或范围剪枝。

block meta 的刷新策略是“从 block 内容重算”，而不是在每次元素移动时逐 slot 增量维护。主要刷新时机：

```text
cacheBlock():
  从 ZNS 读入 block 后刷新。

pma_insert() / pma_delete():
  PMA 内容变化后刷新受影响或已加载 block。

spread() / resize() / pma_compact():
  PMA 重排或扩容后刷新受影响 block。

dumpToZone() / dumpToZoneZoneSerial():
  写回 ZNS 前根据 DRAM buf 刷新；写完后 buf 释放，但 meta 保留。

garbage_collection*() / class-aware GC:
  搬迁有效 block 时根据读出的 payload 刷新。

recover_graph_from_zns_records():
  record 恢复选择 latest DATA record 后，读取 payload 重建 meta。

graph_recovery():
  旧路径根据 global_zone_vec 扫描有效块，读取 payload 后重建 meta。
```

## 4. 写入 ZNS 的两条路径

### 4.1 旧路径：payload-only

默认路径，不加 `--zns-record`。

```text
writeBlockToZone(fd, payload, PM_BLK_SIZE, wp)
readBlockFromZone(fd, payload, PM_BLK_SIZE, wp)
```

`pm_blk.global_wp` 指向 payload 起始地址。`mark_zone()` 用 `ZONE_BLK_SIZE` 语义计算 zone 内 block index。

这个路径用于和新增恢复方案对比。

旧路径是不考虑 recovery 的 baseline：

```text
不写 ZNSBlockHeader
不写 COMMIT record
不写 redo log
不预留 redo/meta zone
```

也就是说，普通运行命令不会产生 redo log：

```bash
sudo ./main_pma_check_new 0 1 parallel2
```

### 4.2 新路径：ZNS record

启用参数：

```bash
--zns-record
```

record 格式定义在 `global.h`：

```cpp
struct ZNSBlockHeader {
    uint32_t magic;
    uint32_t version;
    uint32_t record_type;   // DATA 或 COMMIT
    uint32_t payload_bytes;
    uint64_t epoch;
    int32_t src_vid;
    int32_t pma_blk_idx;
    uint64_t seq;
    uint64_t first_seq;     // COMMIT 使用
    uint64_t last_seq;      // COMMIT 使用
    uint32_t data_record_count;
    uint32_t aggregate_checksum;
    uint32_t valid_count;
    uint32_t pma_ele;
    uint32_t pma_cap;
    uint64_t committed_op_id; // COMMIT 使用
    uint32_t checksum;
    uint32_t tail_magic;
};
```

写入时：

```text
beginZNSEpoch()
writeZNSRecord()
  -> 生成 DATA header
  -> 写入当前 epoch
  -> 计算 checksum
  -> 拼成 512B record
  -> pwrite(fd, record, 512, record_wp)
commitZNSEpoch()
  -> 向 metadata zone 追加 COMMIT record
  -> COMMIT 记录 epoch、seq 范围、data record 数、汇总 checksum 和 committed_op_id
```

读取时：

```text
readZNSRecordPayload()
  -> pread 整条 512B record
  -> 校验 magic/version/checksum
  -> 拷贝 payload
```

record 路径下，`pm_blk.global_wp` 指向 record 起始地址，不是 payload 起始地址。

record 路径是 recovery 版本。加 `--zns-record` 后才会启用 DATA/COMMIT record 和独立 redo log zone：

```bash
sudo ./main_pma_check_new 0 1 parallel2 --zns-record
```

恢复命令 `--recover-zns` 会自动打开 record/recovery 路径，但 redo replay 期间不会再次追加 redo log。

### 4.3 Epoch Commit 事务边界

record 模式下，一个 `dumpToZone*()` 或 `persist_cache_vertices_as_records()` 会被当成一个 epoch：

```text
beginZNSEpoch()
  DATA epoch=N src=a blk=x
  DATA epoch=N src=b blk=y
  DATA epoch=N src=c blk=z
commitZNSEpoch()
  COMMIT epoch=N
```

DATA record 表示“某个顶点 PMA block 的一个版本”。COMMIT record 不属于任何顶点，只表示该 epoch 的 DATA record 已完整写入并可参与恢复。

COMMIT record 写入最后一个 metadata zone，包含：

```text
epoch
first_seq / last_seq
data_record_count
aggregate_checksum
committed_op_id
```

`aggregate_checksum` 是顺序无关的，可以覆盖 DATA record 跨 zone 分布时物理扫描顺序和 seq 顺序不一致的情况。当前 `dumpToZoneParallel()` 为保证 ZNS zone 内顺序写，已收敛到 zone-serial writer；插入阶段仍可并行，物理写入阶段按 zone 顺序提交。

恢复时只有通过以下条件的 epoch 才会生效：

```text
存在合法 COMMIT record
DATA record 数量等于 data_record_count
DATA seq 范围等于 first_seq / last_seq
DATA checksum 汇总等于 aggregate_checksum
```

如果掉电发生在 COMMIT 写入前，或者 COMMIT / DATA 任一部分不完整导致校验不匹配，该 epoch 会被整体忽略。

## 5. Zone-Class 策略

启用参数：

```bash
--zone-class
```

顶点分类：

```cpp
VERTEX_COLD
VERTEX_HOT
VERTEX_HUB
```

分类函数：

```cpp
classify_vertex(Vertex *v)
```

当前规则：

- `hub`：`degree >= hub_degree_threshold` 或 `blk_list.size() >= hub_block_threshold`
- `hot`：`update_count >= hot_update_threshold`
- 其他为 `cold`

zone 划分在 `generate_block_opt()` 中完成：

- 前一部分 zone：hot
- 中间 zone：cold
- 最后一部分 zone：hub

刷盘时通过 `requestClassifiedBlock()` 优先从对应 class zone 分配 block。对应 class 空间不足时会 fallback 到其他 class zone，以保证实验能继续跑。

## 6. GC 模式

三个模式通过参数选择：

```bash
--gc-global
--gc-class-aware
--gc-class-local
```

### 6.1 `--gc-global`

不开 `--zone-class` 时走旧 `garbage_collection()`。

开 `--zone-class` 时，在所有 class zone 中选择 victim zone，使用 record/class-aware 的搬迁逻辑，避免落回旧单队列造成无法重新分配。

### 6.2 `--gc-class-aware`

全局选择 victim zone，但搬迁有效 block 时按顶点当前 class 重新分配目标 zone。

关键点：搬迁后必须更新：

```cpp
v->blk_list[blkIndex].global_wp
mark_zone(src, blkIndex, new_wp)
```

### 6.3 `--gc-class-local`

优先只在当前请求 class 的 zone 集合中选 victim。比如 hot zone 满了，先对 hot zone 做 GC。若失败，再 fallback 到全局 class-aware GC。

### 6.4 GC 并发与 record commit

GC 与查询/更新的并发边界：

```text
1. GC 先持有 global_gc_mutex，保证同一时刻只有一个 GC。
2. 搬迁某个有效 block 前，GC 尝试获取该 block 所属顶点的 unique_lock。
3. 如果顶点正被查询、更新或刷盘占用，GC 放弃本轮 victim，不 reset zone。
4. 只有成功拿到顶点锁后，GC 才会更新 v->blk_list[blk_idx].global_wp 和 zone metadata。
```

record 模式下，GC 搬迁也有 epoch 边界：

```text
beginZNSEpoch()
  写搬迁后的 DATA record
commitZNSEpoch()
  COMMIT 成功后才 reset victim zone
```

如果 GC 中途因为 busy vertex 或分配失败放弃，会先提交已经成功搬迁的 DATA records，再返回失败；victim zone 不会 reset。因此恢复时不会出现“旧 zone 已 reset、新 DATA record 未 COMMIT 而被忽略”的丢块窗口。

## 7. ZNS-Native 恢复

启用写入 record：

```bash
sudo ./main_pma_check_new 0 1 new --zns-record --print-pma
```

从 ZNS record 恢复：

```bash
sudo ./main_pma_check_new 0 1 --recover-zns --print-pma
```

`--recover-zns` 会自动打开 `--zns-record`，并且不会执行：

```cpp
zbd_reset_zones()
```

恢复逻辑：

```text
recover_graph_from_zns_records()
  -> 扫描 metadata zone，找合法 COMMIT record
  -> 扫描 data zone，统计已提交 epoch 的 DATA record
  -> 校验 data_record_count / aggregate_checksum / seq 范围
  -> 只使用通过校验的 epoch
  -> 得到最新 committed_op_id 高水位（当前仅用于统计和后续优化）
  -> key = (src_vid, pma_blk_idx)
  -> 同一 key 保留已提交 epoch 中最新的版本
  -> 重建 global_vectex_vec
  -> 重建 v->blk_list[blk_idx].global_wp
  -> 重建 degree / pma_ele / pma_cap
  -> mark_zone()
  -> 扫描 redo log zone
  -> 保守顺序回放全量 redo log，重建未进入 ZNS 快照的 DRAM-only 小度顶点更新
```

恢复完成后会输出两类统计。第一行是概要：

```text
ZNS-native recovery records: <latest_records>
  committed_epochs: <committed_epochs>
  valid_epochs: <valid_epochs>
  committed_op_id: <committed_op_id>
  redo_records: <redo_records>
  max_vid: <max_vid>
  time: <seconds>
```

第二行是可用于实验分析的详细指标：

```text
ZNS recovery stats:
  meta_slots_scanned=<n>
  commit_records_valid=<n>
  committed_epochs=<n>
  valid_epochs=<n>
  data_slots_scanned=<n>
  data_records_valid=<n>
  data_records_committed_epoch=<n>
  latest_records=<n>
  redo_slots_scanned=<n>
  redo_records_valid=<n>
  redo_records_replayed=<n>
  rebuilt_vertices=<n>
  rebuilt_blocks=<n>
  rebuilt_edges_before_redo=<n>
  committed_op_id=<n>
  max_seq=<n>
  max_epoch=<n>
  max_op_id=<n>
  recovery_time_sec=<seconds>
```

各项含义：

```text
meta_slots_scanned:
  扫描 metadata zone 的 record slot 数，反映 metadata zone 扫描成本。

commit_records_valid:
  checksum/header 合法且类型为 COMMIT 的 record 数。

committed_epochs:
  metadata zone 中出现过的 epoch 数；同一 epoch 多条 COMMIT 时保留 seq 最大的一条。

valid_epochs:
  通过 COMMIT 校验的 epoch 数。只有这些 epoch 的 DATA record 会参与恢复。

data_slots_scanned:
  扫描 data zones 的 record slot 数，反映全盘 DATA 扫描成本。

data_records_valid:
  header/checksum 合法且类型为 DATA 的 record 数。

data_records_committed_epoch:
  属于已出现 COMMIT epoch 的 DATA record 数。之后还会用数量、seq 范围和 aggregate checksum 过滤出 valid epoch。

latest_records:
  按 (src_vid, pma_blk_idx) 去重后保留下来的最新 DATA record 数，也就是恢复快照中的有效块版本数。

redo_slots_scanned:
  扫描 redo log zone 的 record slot 数，反映 redo zone 扫描成本。

redo_records_valid:
  header/checksum 合法的 redo log record 数。

redo_records_replayed:
  实际回放的 redo log 数。当前为保守全量回放，主要用于补齐 cache-only 小度顶点更新。

rebuilt_vertices:
  恢复过程中重建的 Vertex 数。

rebuilt_blocks:
  从 latest_records 恢复出的 PMA/cache record 块数。

rebuilt_edges_before_redo:
  redo 回放前，仅由 ZNS DATA record 恢复出的邻边数量。

committed_op_id:
  COMMIT record 中记录的最大更新高水位。当前主要用于统计和后续 checkpoint 优化。

max_seq:
  扫描到的最大 ZNS DATA/COMMIT record 序号。

max_epoch:
  扫描到的最大 epoch。

max_op_id:
  扫描到的最大 redo log op_id。

recovery_time_sec:
  recover_graph_from_zns_records() 的总耗时。
```

一致性语义：

- 新版本写完整后，即使旧版本仍然在盘上，恢复时也会选择 `seq` 最大的 record。
- 掉电发生在 DATA record 写入过程中时，checksum/tail magic 不匹配的 record 会被恢复扫描忽略。
- 掉电发生在一个 epoch 的 COMMIT 写入前时，该 epoch 的 DATA record 整批不会参与恢复。
- 掉电发生在 COMMIT 完整写入后时，恢复会用 COMMIT 中的数量、seq 范围和汇总 checksum 验证该 epoch 是否完整。
- redo log 使用独立预留 zone，当前恢复采用保守全量回放策略。
- 由于小度顶点可能长期只存在 DRAM cache，只有在 checkpoint 明确持久化 cache-only 顶点后，`committed_op_id` 才能作为跳过 redo log 的安全高水位；当前实现采用保守策略，按顺序回放 redo log，利用 INSERT/DELETE 的幂等性保证正确性。
- 当前实现已有 metadata zone 的 epoch commit，但还没有 checkpoint summary，恢复仍需要扫描 data zone。
- `ZNS_BLOCK_VERSION` 当前为 3；旧版 record 不会被新版 `--recover-zns` 识别。

## 8. 并发查询与更新

当前代码已经把图计算查询路径从“会污染更新缓存的读路径”中拆出来。`benchmark.h` 中的 PageRank、BFS、CC 现在通过：

```cpp
for_each_neighbor_readonly(Vertex *v, Fn &&fn)
```

统一遍历邻接。这个 helper 的规则是：

```text
1. 对当前顶点持有 shared lock，防止同一顶点的边更新或刷盘同时改 degree/cache_buf/blk_list/buf/global_wp。
2. 如果顶点仍是小度 cache-only 状态，只读该顶点版本中的 degree 和 cache_buf。
3. 如果顶点已经迁入 PMA，遍历 blk_list。
4. 对已经在 DRAM 的 block，直接只读 buf。
5. 对只在 ZNS 上的 block，读到查询私有临时 buffer，处理完立即复用/释放。
6. 查询路径不调用 cacheBlock()，因此不会 invalid 旧 ZNS block，不会增加 global_cache_blk，也不会触发 dumpToZone()。
```

对应地，边插入/删除仍使用 `Vertex::mutex` 的 `unique_lock`。`dumpToZone()`、`dumpToZoneZoneSerial()` 和 `persist_cache_vertices_as_records()` 在释放 `buf`、更新 `global_wp` 或把 cache-only 顶点持久化为 record 前，也会持有该顶点的 `unique_lock`。GC 搬迁有效 block 时会先持有 `global_gc_mutex` 串行化 GC，再对被搬迁 block 所属顶点尝试获取 `unique_lock`；如果顶点正被查询、更新或刷盘占用，本轮 GC 会放弃当前 victim，不 reset zone。这样避免 GC 与查询/更新同时改同一顶点的 `blk_list/global_wp`，也避免自锁。

因此当前实现保证的是：

```text
单个顶点读取期间的内存安全和结构一致性：
  查询不会读到同一顶点被更新/刷盘改到一半的 degree/cache_buf/blk_list/buf/global_wp。

查询读 ZNS block 不产生写放大：
  查询从 ZNS 读入的 block 是临时只读 buffer，用完即丢，不进入 dirty cache。

GC 与查询/更新的顶点级互斥：
  GC 只有拿到目标顶点 unique_lock 后才会搬迁该顶点 block 并更新 global_wp。
```

当前实现还不是完整的全图 MVCC snapshot。也就是说，长查询可能在不同顶点上看到不同提交时刻的状态，但每个顶点自身的读取是互斥保护的。完整快照一致性需要后续引入：

```text
global_graph_epoch:
  查询开始时记录 snapshot_epoch。

VertexVersion:
  版本化 degree/cache_buf/PMA 元数据/blk_list，用于覆盖低度顶点和 cache -> PMA 迁移。

BlockVersion:
  版本化 PMA block 的 buf/global_wp/meta，被更新 block 走 copy-on-write。

active query epochs:
  旧 VertexVersion、BlockVersion 和旧 ZNS physical block 延迟回收。
```

各种场景的目标语义如下：

```text
低度顶点更新:
  复制 VertexVersion，旧查询读旧 cache_buf，新查询读新 cache_buf。

低度顶点迁移 PMA:
  旧 VertexVersion 保持 cache-only，新 VertexVersion 指向 PMA block。

DRAM cached block 更新:
  复制旧 buf，在新 buf 上修改，旧 buf 保留给旧查询。

ZNS-only block 更新:
  从旧 global_wp 读 payload，生成新 BlockVersion，旧 global_wp 延迟 invalid/GC。

PMA spread/resize/compact:
  所有被写入的 block 都生成新 BlockVersion，未改 block 共享旧版本。

无向边更新:
  u->v 和 v->u 应使用同一个 commit_epoch 原子发布，避免查询看到半条无向边。
```

因此，当前提交解决的是查询读路径的安全性和读写缓存分离；如果实验需要“一个长查询看到一个固定全图快照”，下一步需要按上面的 `VertexVersion + BlockVersion + epoch` 方案继续实现。

## 9. 命令行参数

基本格式：

```bash
sudo ./main_pma_check_new <dataset_id> <thread_num> [mode] [options]
```

构图模式：

```text
baseline
new
parallel
parallel2
```

常用参数：

```text
--print-pma
--insert <src> <dst>
--delete <src> <dst>
--compact-delete
--no-cache
--cache-threshold <n>
--zone-writer
--zone-class
--gc-global
--gc-class-aware
--gc-class-local
--zns-record
--recover-zns
```

设备路径固定在：

```cpp
#define NVME_DEVICE_PATH "/dev/nvme0n1"
```

## 10. 常用验证命令

编译：

```bash
make
```

旧路径 demo：

```bash
sudo ./main_pma_check_new 0 1 new --print-pma
```

无 recovery baseline：

```bash
sudo ./main_pma_check_new 0 1 parallel2 --print-pma
```

该命令不写 redo log，也不会预留 redo/meta zone。

zone-class 路径：

```bash
sudo ./main_pma_check_new 0 1 new --zone-class --print-pma
```

三种 GC 参数形式：

```bash
sudo ./main_pma_check_new 0 1 new --zone-class --gc-global --print-pma
sudo ./main_pma_check_new 0 1 new --zone-class --gc-class-aware --print-pma
sudo ./main_pma_check_new 0 1 new --zone-class --gc-class-local --print-pma
```

删除边：

```bash
sudo ./main_pma_check_new 0 1 new --delete 1 8 --print-pma
sudo ./main_pma_check_new 0 1 new --delete 1 8 --compact-delete --print-pma
```

ZNS-native record 写入和恢复：

```bash
sudo ./main_pma_check_new 0 1 new --zns-record --print-pma
sudo ./main_pma_check_new 0 1 --recover-zns --print-pma
```

小度 cache 恢复对比实验：

```bash
# 方案 1：关闭 cache_buf，顶点第一次插入就进入 PMA/ZNS 路径
sudo ./main_pma_check_new 0 1 parallel2 --zns-record --no-cache --print-pma
sudo ./main_pma_check_new 0 1 --recover-zns --print-pma

# 方案 2：保留默认小度 DRAM cache，依靠 redo log 保守回放补齐 cache-only 更新
sudo ./main_pma_check_new 0 1 parallel2 --zns-record --print-pma
sudo ./main_pma_check_new 0 1 --recover-zns --print-pma

# 可选：显式设置 cache 阈值，0 等价于 --no-cache
sudo ./main_pma_check_new 0 1 parallel2 --zns-record --cache-threshold 0 --print-pma
sudo ./main_pma_check_new 0 1 parallel2 --zns-record --cache-threshold 4 --print-pma
```

两个方案的核心差异：

```text
方案 1 no-cache:
  边更新 -> redo log -> PMA dirty block -> ZNS DATA record -> COMMIT
  优点：恢复高水位语义简单，cache-only 丢失问题基本不存在
  缺点：小度顶点更早占用 ZNS block，写放大和空间开销更高

方案 2 cache + redo:
  边更新 -> redo log -> cache_buf 或 PMA dirty block
  cache-only 顶点可能不在 ZNS DATA record 中
  恢复时：先恢复 ZNS 块级快照，再顺序回放 redo log
  优点：减少小度顶点写放大
  缺点：恢复需要更多 redo 查询/回放
```

## 11. 当前创新点

### 11.1 ZNS-aware PMA Block Mapping

系统将 ZNS zone 切分为固定大小的可写 block/record，并把每个顶点的 PMA 邻接表映射到一组 ZNS block 上。每个 PMA block 通过 `pm_blk.global_wp` 指向当前最新物理版本；record 模式下还通过 `src_vid` 和 `pma_blk_idx` 将该映射持久化。这样 PMA 插入、删除和重平衡产生的局部随机更新可以转化为 ZNS 友好的 append-only block rewrite，同时为恢复和 GC 提供明确的逻辑-物理映射。

### 11.2 小度顶点 DRAM cache + PMA 延迟刷盘

小度顶点不立即占用 ZNS block，减少小邻接表写放大。超过阈值后再迁移到 PMA。

### 11.3 PMA Block Meta 查询加速

每个 `pm_blk` 维护内存态 `valid_count`、`valid_bitmap`、`min_val` 和 `max_val`。查询时可以跳过空 block，并用 bitmap 快速定位 block 内最近有效 slot，减少 PMA 稀疏布局中反复扫描 `-1` 空洞的开销。meta 不持久化，恢复时从 block payload 重建，因此 payload-only 和 ZNS record 两条路径都能使用同一套查询辅助信息。

### 11.4 删除边延迟写回

删除边先在 cache/PMA DRAM block 中修改。达到 `global_cache_blk >= cache_blk_threshold` 后再刷回 ZNS，避免每个删除都立即读写 ZNS。

### 11.5 Vertex-Class to Zone-Class

根据顶点更新热度和度数将顶点映射到 hot/cold/hub zone，减少冷热数据混放，为后续 class-aware GC 和磨损/写放大分析提供实验基础。

### 11.6 三种 GC 策略对比

保留全局 GC，同时提供 class-aware 和 class-local GC：

- 全局清理：作为基线。
- class-aware 搬迁：减少恢复后布局混乱。
- class-local 优先：hot zone 满时优先清理 hot zone，而不是立即全盘 GC。

### 11.7 ZNS-native record 恢复

每条 DATA record 自描述，恢复不依赖运行时 `global_zone_vec`。record 模式下最后一个 zone 作为 metadata zone，批量写入完成后追加 COMMIT record。恢复只使用通过 COMMIT 校验的 epoch，再通过 `epoch/seq` 选择最新版本，为掉电恢复提供基本事务边界。

### 11.8 独立 redo log zone + 保守回放

record 模式下预留独立 redo log zone，边插入/删除先追加 `EdgeLogRecord`，再修改 DRAM 图结构。恢复时先恢复已提交的 ZNS 块级快照，再按 `op_id` 顺序回放 redo log。该设计覆盖 cache-only 小度顶点未写入 ZNS 的情况，同时支持 `--no-cache` 与默认 cache 策略的恢复开销/写放大对比。

## 12. 当前限制和注意事项

1. `--recover-zns` 仍需要扫描 data zone，没有 checkpoint summary，盘大时会慢。
2. `ZNS_RECORD_SIZE` 当前固定为 512B，payload 仍是 64 个 `int32_t`，剩余空间用于 header 和 padding。
3. record 恢复目前重建的是 block 指针、degree 和 PMA 元数据，不会重新加载原 dataset；增量更新依赖 redo log zone。
4. 旧路径和 record 路径的 `global_wp` 语义不同，新增代码必须通过 `zns_alloc_unit_size()` 和 record-aware 函数处理，不要手写 `ZONE_BLK_SIZE`。
5. `--recover-zns` 不 reset zone；普通构图命令会 reset `/dev/nvme0n1`。
6. 当前 redo 恢复为了覆盖 cache-only 顶点采用全量回放；如果后续实现完整 cache checkpoint，可改为只回放 `op_id > committed_op_id`。
7. GC 搬迁时需要始终更新顶点 `blk_list` 和 zone metadata，否则恢复和后续查询都会读到旧位置。
8. record header 当前是 version 3，升级前写入的旧版 record 需要旧程序恢复，或后续补兼容读取逻辑。
9. `pm_blk` 的 block meta 是内存态派生索引，不持久化。新增会改变 block payload 的路径时，需要调用 `refresh_pm_blk_meta()` 或确保恢复/读入时重建。
10. 当前并发查询实现提供顶点级读写互斥和只读 ZNS block streaming，不提供全图固定 snapshot；需要严格 snapshot 的实验应继续实现 `VertexVersion + BlockVersion + global_graph_epoch`。
11. GC 若遇到 busy vertex 会放弃当前 victim，不 reset zone；这保证并发安全，但可能导致本次空间回收失败，调用方需要有重试或错误处理。
12. record 模式下 GC 搬迁已经在 reset 前提交 DATA epoch，但还没有实现 MVCC 旧版本延迟 reset；如果未来引入全图 snapshot 查询，需要把 reset 延迟到没有活跃查询引用旧 epoch 之后。

## 13. 下次继续开发建议

下一步可以优先做：

1. 给 `--recover-zns` 增加更详细的统计输出：扫描 record 数、无效 header 数、checksum fail 数、每个 zone live record 数。
2. 增加 `ZoneSummary`，避免每次恢复扫描全部 data zone。
3. 增加 crash 模拟开关，例如写 header 后退出、写 payload 后退出，验证 checksum 和 seq 选择逻辑。
4. 对 record 模式下 GC 后恢复做专项验证，确认搬迁后的 `seq` 和 `global_wp` 正确。
5. 把 hot/cold/hub 阈值做成命令行参数，方便实验扫描。
6. 将 `min_val / max_val` 接入 block-level 查询剪枝，先在 block meta 上定位候选 block，再读取 payload。
7. 实现完整并发快照：`global_graph_epoch`、`VertexVersion`、`BlockVersion`、活跃查询 epoch 表、旧 ZNS block 延迟 invalid/GC。
