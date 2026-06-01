# ZoneGraph README

ZoneGraph 是一个面向 ZNS SSD 的动态图存储原型系统。它面向频繁边插入、边删除和图查询的动态图场景，目标是在保持图更新效率的同时，将物理写入组织成 ZNS 友好的顺序块写。

系统核心思想：

```text
边级图更新
  -> DRAM 中的 degree-aware 维护
  -> PMA 风格邻接块
  -> ZNS 上的 append-only record
```

本文档只保留系统框架、实验设置和写作相关说明，不再维护函数级代码导览。

## 1. 系统概览

ZoneGraph 以顶点为单位组织邻接数据。每个顶点根据度数选择不同表示：

- 低度顶点：邻居保留在 DRAM cache 中；
- 高度顶点：邻居迁移到 PMA 风格的稀疏邻接块；
- 邻接块：作为固定大小的逻辑块进行维护；
- 持久化：dirty block 以 append-only ZNS record 写入设备。

核心映射关系：

```text
(vertex id, logical adjacency block id) -> latest ZNS record
```

更新不会原地覆盖旧 record，而是追加写新版本。旧版本变为 invalid block，后续由 GC 回收。

## 2. 主要设计

### 2.1 Degree-Aware 邻接管理

真实图中大量顶点度数很小。ZoneGraph 避免将这些低度顶点过早物化为 ZNS block，而是优先保留在 DRAM 中。这样可以降低 dump 次数、减少 ZNS 写入，并改善低度顶点查询路径。

高度顶点则使用 PMA 风格邻接块，使插入、删除和查询都能在较稳定的块级结构上执行。

### 2.2 PMA 风格块级维护

高度顶点的邻接表被划分成带空隙的 PMA block。边插入和删除优先影响局部块或局部窗口，避免每次更新都重写完整邻接表。

这个设计将大量边级更新转化为较少的块级 rewrite，是 ZoneGraph 连接动态图语义和 ZNS 写入约束的关键。

### 2.3 ZNS 友好的追加写

ZoneGraph 将修改后的邻接块作为 ZNS record 追加写入。每个 record 自带逻辑身份信息，因此系统可以在运行时和 recovery 时恢复逻辑块到物理位置的映射。

### 2.4 Block-Level GC

当邻接块被更新后，旧物理版本失效。GC 以 zone 为回收单位，搬迁仍然 live 的邻接块，并 reset victim zone。

当前实验比较：

- Basic GC；
- Greedy GC；
- ATZ-GC。

ATZ-GC 通过冷热感知减少 cold live block 被反复搬迁。

### 2.5 Restart Recovery

Recovery 从 ZNS record 重建 DRAM 中的图 metadata 和逻辑到物理映射，不需要重新加载原始图文件。

当前 recovery 实验关注 restart metadata recovery，统计恢复总时间并进行正确性验证。它不是完整掉电一致性协议实验。

### 2.6 Edge Deletion

删除边时，ZoneGraph 先在 DRAM/PMA 邻接结构中完成逻辑删除，物理持久化延后到块级 dump。删除实验按固定 directed delete attempts 统计时间。

## 3. 实验集合

### 3.1 Baseline 对比

对比系统：

- ZoneGraph；
- LiveGraph；
- LSMGraph；
- MBFGraph。

主要数据集：

- `com-DBLP`
- `com-Youtube`
- `twitch-gamers`
- `as-Skitter`
- `soc-Pokec`
- `soc-LiveJournal1`

指标：

- graph maintenance time；
- BFS time；
- PageRank time；
- WCC time；
- ZoneGraph before/after final dump query time。

### 3.2 Cache Block Threshold 消融

变量：

```text
50000
100000
200000
1000000
```

目的：

- 分析 dump 频率对 maintenance/query 的影响；
- 说明过大的阈值会使系统偏向 DRAM-only；
- 支撑主实验中 cache block threshold 的选择。

### 3.3 Degree-Aware 消融

变量：

```text
0
8
16
```

含义：

- `0`：关闭低度顶点 cache；
- `8`：较小低度 cache；
- `16`：默认设计。

该实验用于验证低度顶点 DRAM cache 的必要性。

### 3.4 Graph Maintenance 并行度消融

变量：

```text
1
4
8
64
```

机器最大逻辑线程数为 8。`64` 用于保持与主实验参数一致，并观察过量线程下的表现。

### 3.5 GC 策略实验

比较：

- Basic；
- Greedy；
- ATZ-GC。

指标：

- update total time；
- GC count；
- GC total time；
- moved live blocks；
- write amplification。

Skitter 是最干净的 update 阶段 GC 对比；Pokec 和 LiveJournal 可作为大图补充。

### 3.6 Recovery 实验

Recovery 实验只需要报告：

- recovery total time；
- validation time/result；
- recovery vs load。

该实验不需要运行 graph query。

### 3.7 删除效率实验

删除 workload：

- 每个数据集固定 `1,000,000` 条无向边；
- 每条无向边执行两个方向删除；
- 即 `2,000,000` 次 directed delete attempts。

对比说明：

- ZoneGraph：在线逻辑删除，并发执行；
- LiveGraph：事务删除；
- LSMGraph：tombstone-style 删除；
- MBFGraph：文件级 marker deletion，应单独标注。

## 4. 实验约定

### Directed / Undirected

主实验默认将输入边视为无向边：

```text
1 input edge -> 2 directed edges
```

删除实验同样遵循：

```text
1 undirected delete -> 2 directed delete attempts
```

### Final Dump

ZoneGraph 查询结果区分：

- before final dump；
- after final dump。

before final dump 可能仍有低度顶点停留在 DRAM；after final dump 更接近完整持久化后的查询路径。

### GC 实验边界

普通 baseline maintenance/query 实验不是 GC 压力实验。GC 效果应使用专门的 skew update workload 说明。

### MBFGraph 删除语义

MBFGraph 的删除是对二进制边文件写 deletion marker，和 ZoneGraph/LiveGraph/LSMGraph 的在线图结构删除不同。论文中需要单独标注。

## 5. 可支撑的主要结论

1. Degree-aware 邻接管理能显著减少低度顶点带来的 ZNS 写入和 dump 开销。
2. PMA 风格邻接块能在更新效率和查询友好性之间取得平衡。
3. ZNS record 将随机图更新转化为 append-only 块写。
4. ATZ-GC 在 skewed update 场景下降低 live block 搬迁量和 GC 时间。
5. Recovery 可以从 ZNS record 重建 DRAM metadata，不需要重新加载原始图。
6. ZoneGraph 删除效率高，因为删除先作用于 DRAM/PMA 逻辑结构，物理持久化被延后到块级 dump。

## 6. 限制说明

- 当前查询和结构更新不并发执行；查询阶段会获得稳定图视图。
- Recovery 实验是 restart metadata recovery，不等同于完整掉电一致性验证。
- MBFGraph 删除结果不应和在线删除系统直接混同。
- 过大的 cache block threshold 会使系统更接近 DRAM-only，只适合作为上界消融，不适合作为公平主结果。

