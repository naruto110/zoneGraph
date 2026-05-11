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

## 2. 主要文件

`main.cpp`

程序入口、命令行参数解析、设备打开、zone reset、构图模式选择、额外插入/删除、打印、恢复入口。

`global.h`

全局常量、核心结构体和全局变量。重要结构包括：

- `Vertex`：顶点状态，包含 `degree`、`cache_buf`、`blk_list`、`pma`、`update_count`、`vclass`。
- `pm_blk`：顶点 PMA 的一个物理块描述，包含 `global_wp` 和可选 DRAM `buf`。
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

这是为了让 `data.txt` 这样的小数据集也能触发刷写到 ZNS，便于验证。恢复实验中，`--zns-record` 会把 cache-only 顶点也持久化为 record，否则这些顶点不会出现在盘上。

### 3.2 PMA block

当 `cache_buf` 满后，已有 cache 邻居整体迁移到 PMA：

```cpp
Vertex::blk_list
Vertex::pma
```

PMA 逻辑上仍以 `PM_BLK_SIZE=64` 个 `int32_t` 为一个块。空位用 `-1` 表示。

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

## 4. 写入 ZNS 的两条路径

### 4.1 旧路径：payload-only

默认路径，不加 `--zns-record`。

```text
writeBlockToZone(fd, payload, PM_BLK_SIZE, wp)
readBlockFromZone(fd, payload, PM_BLK_SIZE, wp)
```

`pm_blk.global_wp` 指向 payload 起始地址。`mark_zone()` 用 `ZONE_BLK_SIZE` 语义计算 zone 内 block index。

这个路径用于和新增恢复方案对比。

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
    int32_t src_vid;
    int32_t pma_blk_idx;
    uint64_t seq;
    uint32_t valid_count;
    uint32_t pma_ele;
    uint32_t pma_cap;
    uint32_t checksum;
};
```

写入时：

```text
writeZNSRecord()
  -> 生成 header
  -> 计算 checksum
  -> 拼成 512B record
  -> pwrite(fd, record, 512, record_wp)
```

读取时：

```text
readZNSRecordPayload()
  -> pread 整条 512B record
  -> 校验 magic/version/checksum
  -> 拷贝 payload
```

record 路径下，`pm_blk.global_wp` 指向 record 起始地址，不是 payload 起始地址。

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
  -> 顺序扫描所有 zone
  -> 按 ZNS_RECORD_SIZE 读取 record
  -> 校验 magic/version/checksum
  -> key = (src_vid, pma_blk_idx)
  -> 同一 key 保留 seq 最大的版本
  -> 重建 global_vectex_vec
  -> 重建 v->blk_list[blk_idx].global_wp
  -> 重建 degree / pma_ele / pma_cap
  -> mark_zone()
```

一致性语义：

- 新版本写完整后，即使旧版本仍然在盘上，恢复时也会选择 `seq` 最大的 record。
- 如果写入不完整或 checksum 不匹配，恢复时会忽略或报错。
- 当前实现是第一版扫描恢复，还没有单独 metadata zone 或 checkpoint summary。

## 8. 命令行参数

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

## 9. 常用验证命令

编译：

```bash
make
```

旧路径 demo：

```bash
sudo ./main_pma_check_new 0 1 new --print-pma
```

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

## 10. 当前创新点

### 10.1 小度顶点 DRAM cache + PMA 延迟刷盘

小度顶点不立即占用 ZNS block，减少小邻接表写放大。超过阈值后再迁移到 PMA。

### 10.2 删除边延迟写回

删除边先在 cache/PMA DRAM block 中修改。达到 `global_cache_blk >= cache_blk_threshold` 后再刷回 ZNS，避免每个删除都立即读写 ZNS。

### 10.3 Vertex-Class to Zone-Class

根据顶点更新热度和度数将顶点映射到 hot/cold/hub zone，减少冷热数据混放，为后续 class-aware GC 和磨损/写放大分析提供实验基础。

### 10.4 三种 GC 策略对比

保留全局 GC，同时提供 class-aware 和 class-local GC：

- 全局清理：作为基线。
- class-aware 搬迁：减少恢复后布局混乱。
- class-local 优先：hot zone 满时优先清理 hot zone，而不是立即全盘 GC。

### 10.5 ZNS-native record 恢复

每条持久化 record 自描述，恢复不依赖运行时 `global_zone_vec`。通过 `seq` 选择最新版本，为 crash consistency 提供基础。

## 11. 当前限制和注意事项

1. `--recover-zns` 是第一版全盘扫描恢复，没有 checkpoint summary，盘大时会慢。
2. `ZNS_RECORD_SIZE` 当前固定为 512B，payload 仍是 64 个 `int32_t`，剩余空间用于 header 和 padding。
3. record 恢复目前重建的是 block 指针、degree 和 PMA 元数据，不会重新加载原 dataset。
4. 旧路径和 record 路径的 `global_wp` 语义不同，新增代码必须通过 `zns_alloc_unit_size()` 和 record-aware 函数处理，不要手写 `ZONE_BLK_SIZE`。
5. `--recover-zns` 不 reset zone；普通构图命令会 reset `/dev/nvme0n1`。
6. 当前并行插入路径仍主要用于已有实验，后续若要做严格一致性实验，建议先以 `new` 串行路径验证。
7. GC 搬迁时需要始终更新顶点 `blk_list` 和 zone metadata，否则恢复和后续查询都会读到旧位置。

## 12. 下次继续开发建议

下一步可以优先做：

1. 给 `--recover-zns` 增加更详细的统计输出：扫描 record 数、无效 header 数、checksum fail 数、每个 zone live record 数。
2. 增加 `ZoneSummary` 或 metadata zone，避免每次恢复全盘扫描。
3. 增加 crash 模拟开关，例如写 header 后退出、写 payload 后退出，验证 checksum 和 seq 选择逻辑。
4. 对 record 模式下 GC 后恢复做专项验证，确认搬迁后的 `seq` 和 `global_wp` 正确。
5. 把 hot/cold/hub 阈值做成命令行参数，方便实验扫描。

