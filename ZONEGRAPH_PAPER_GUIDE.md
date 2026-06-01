# ZoneGraph 论文写作指南

本文档用于指导论文写作，重点说明 ZoneGraph 的系统框架、设计动机、创新点和实验叙事。这里不展开代码实现细节，也不按函数介绍系统。

## 1. 核心定位

ZoneGraph 是一个面向 ZNS SSD 的动态图存储系统。它解决的核心矛盾是：

- 动态图更新是细粒度的，通常以边为单位发生；
- ZNS SSD 需要顺序写入，并且更适合块级、zone 级管理；
- 如果每条边更新都直接落盘，会造成严重写放大和 zone 利用率问题。

ZoneGraph 的核心思想是将逻辑更新粒度和物理持久化粒度解耦：

```text
边级插入/删除
  -> DRAM 中的 degree-aware 维护
  -> PMA 风格邻接块
  -> ZNS 上的 append-only block record
```

论文中应强调：ZoneGraph 不是把每条边直接写到 SSD，而是先在内存中吸收、合并、重组更新，再以 ZNS 友好的块粒度持久化。

## 2. 系统框架

ZoneGraph 以顶点为基本组织单位。每个顶点维护自己的邻接对象，根据顶点度数采用不同表示：

- 低度顶点使用 DRAM 中的紧凑邻接缓存；
- 高度顶点迁移到 PMA 风格的稀疏邻接块；
- PMA 邻接表被切分为固定大小的逻辑块；
- 被修改的逻辑块在达到阈值后追加写入 ZNS。

系统维护的核心映射是：

```text
(vertex id, logical adjacency block id) -> latest ZNS physical record
```

这个映射支撑了三件事：

- 更新时不需要原地覆盖 ZNS 数据；
- GC 可以在块粒度搬迁 live block；
- recovery 可以通过扫描 ZNS record 重建 DRAM metadata。

## 3. 创新点一：Degree-Aware 邻接管理

### 动机

真实图通常具有明显的度数倾斜：大量顶点度数很小，少量顶点拥有很大的邻接表。统一处理所有顶点会带来两个问题：

- 低度顶点过早落盘会产生大量小块写入和无效 ZNS 占用；
- 高度顶点如果只用简单 append/list，会影响更新、删除和查询效率。

### 设计

ZoneGraph 将低度顶点和高度顶点分开处理：

- 低度顶点优先停留在 DRAM cache 中；
- 当顶点度数超过阈值后，再迁移到 PMA 邻接块；
- 高度顶点使用带空隙的 PMA 布局，使插入和删除尽量局部化。

### 写作重点

不要把这部分简单写成“用了 cache”。更准确的表述是：

> ZoneGraph uses degree-aware adjacency placement to avoid materializing small vertices as ZNS blocks while still providing structured, update-friendly storage for high-degree vertices.

### 实验支撑

对应 degree-aware 消融实验：

- `min_degree=0`：关闭低度顶点 cache；
- `min_degree=8`：较小 cache 阈值；
- `min_degree=16`：默认设计。

实验应说明：关闭或削弱低度顶点 cache 会显著增加 dump 次数、dump 时间和维护成本，尤其在 LiveJournal 这类大图上退化明显。

## 4. 创新点二：PMA 风格的块级邻接维护

### 动机

动态图系统需要频繁插入和删除边，但图查询又希望邻接表具有较好的顺序扫描性能。简单 append 会让查询和去重变慢；完全紧凑的有序数组又会让更新代价过高。

### 设计

ZoneGraph 对高度顶点使用 PMA 风格的稀疏邻接布局：

- 邻接表被划分为固定大小的逻辑块；
- 块内和相邻窗口保留空隙；
- 插入/删除优先修改局部块或局部窗口；
- 多个细粒度更新最终合并为少量块级 rewrite。

### 写作重点

PMA 在 ZoneGraph 中的作用是连接图语义和 ZNS 物理约束：

```text
图层面：边级更新
内存层面：PMA 局部维护
存储层面：块级追加写
```

论文可以强调：

> ZoneGraph preserves edge-level update semantics while persisting updates at block granularity.

## 5. 创新点三：面向 ZNS 的块级持久化

### 动机

ZNS SSD 不支持任意位置覆盖写。邻接块更新后，旧物理位置不能原地修改，只能追加写新版本，并将旧版本标记为无效。

### 设计

ZoneGraph 将每个被修改的邻接块写成自描述的 ZNS record。record 中包含恢复逻辑映射所需的元信息，例如所属顶点和逻辑块编号。

因此，系统中的持久化单元不是边，而是邻接块：

```text
logical adjacency block -> ZNS record
```

### 写作重点

这部分应突出“逻辑图结构”和“物理 ZNS 布局”的解耦：

- 图更新只改变逻辑邻接块；
- 物理层通过追加写产生新版本；
- 旧版本成为 GC 可回收对象；
- recovery 可通过 record 重建逻辑到物理的映射。

## 6. 创新点四：Graph-Aware Zone GC

### 动机

随着块更新不断追加写入，ZNS zone 中会混合 live block 和 invalid block。回收一个 zone 需要搬迁其中仍然有效的 block，然后 reset zone。GC 策略如果选择不当，会搬迁大量冷数据，导致写放大。

### 设计

ZoneGraph 的 GC 以图块为单位工作：

- 选择 victim zone；
- 搬迁其中仍然 live 的邻接块；
- 更新逻辑块到新物理位置的映射；
- reset 原 victim zone。

实验中比较三种策略：

- Basic GC：朴素 victim 选择；
- Greedy GC：优先选择 invalid block 多的 zone；
- ATZ-GC：结合冷热信息，减少 cold live block 被反复搬迁。

### 写作重点

ATZ-GC 的核心不是主动把所有数据重新分类搬迁，而是非对称地处理冷热数据：

- hot 数据未来更可能被更新，自然会迁移到新的写入区域；
- cold 数据应该尽量避免被反复 GC 搬迁；
- victim 选择需要同时考虑 invalid 空间收益和 live block 温度。

### 实验支撑

GC 实验应重点报告：

- update total time；
- GC count；
- GC total time；
- moved live blocks；
- write amplification。

Skitter 是最干净的 update-phase GC 对比；Pokec 和 LiveJournal 可作为大图补充，说明 GC 压力可能前移到 load 阶段。

## 7. 创新点五：基于 ZNS Record 的 Restart Recovery

### 动机

系统重启后，DRAM 中的顶点 metadata、PMA metadata 和逻辑到物理映射都会丢失。传统做法可能需要重新加载原始图数据，代价很高。

### 设计

ZoneGraph 的 ZNS record 是自描述的。Recovery 阶段扫描持久化 record，选择每个逻辑邻接块的最新有效版本，并重建：

- 顶点邻接 metadata；
- PMA block metadata；
- 逻辑块到物理 record 的映射。

### 写作边界

Recovery 实验应描述为 restart metadata recovery，而不是完整掉电一致性协议。当前实验关注：

- 清理 DRAM metadata；
- 从 ZNS record 重建图结构；
- 验证恢复后图状态；
- 统计恢复总时间。

最终结果中只需要强调 recovery total time 和 validation result。不要把未进入最终实验路径的日志机制作为主要贡献。

## 8. 删除实验的写法

ZoneGraph 的删除也是一种逻辑更新：

```text
先在 DRAM/PMA 邻接结构中删除边
后续再通过块级 dump 持久化
```

删除实验应按固定 directed delete attempts 统计时间。建议写清楚：

- 1 条无向边删除等价于 2 次 directed delete；
- ZoneGraph 使用并发删除；
- LiveGraph 是事务删除；
- LSMGraph 是 tombstone-style 删除；
- MBFGraph 是文件级 marker deletion，不是在线图结构删除，应单独标注。

删除后 query correctness 可以作为额外验证，不应混入 delete latency。

## 9. 推荐实验结构

论文实验部分建议按如下顺序组织：

1. 与 baseline 的 graph maintenance 和 query 对比；
2. cache block threshold 消融；
3. degree-aware 消融；
4. maintenance 并行度消融；
5. GC 策略实验；
6. recovery 时间实验；
7. 删除效率实验。

每组实验应明确区分：

- maintenance time；
- dump time；
- query time；
- GC time；
- recovery time；
- delete time。

避免把不同阶段混成一个总数，除非表头明确说明。

## 10. 推荐论文结构

1. Introduction
   - 动态图更新和 ZNS 顺序写约束之间的矛盾；
   - ZoneGraph 的核心思想：边级更新、块级持久化。

2. Background and Motivation
   - 动态图邻接维护；
   - PMA sparse array；
   - ZNS append/reset 模型。

3. System Overview
   - degree-aware 顶点表示；
   - PMA 邻接块；
   - ZNS record；
   - 逻辑到物理映射。

4. Graph Maintenance
   - 插入、删除、低度 cache、PMA block update、dump。

5. Garbage Collection
   - invalid block 产生；
   - victim selection；
   - ATZ-GC 冷热感知。

6. Recovery
   - self-describing records；
   - metadata rebuild；
   - validation。

7. Evaluation
   - baseline 对比；
   - 消融；
   - GC；
   - recovery；
   - delete。

## 11. 术语表

建议全文统一使用以下术语：

- directed edge；
- undirected edge as two directed edges；
- low-degree vertex cache；
- PMA adjacency block；
- dirty block；
- block dump；
- ZNS record；
- logical adjacency block；
- physical ZNS location；
- logical-to-physical mapping；
- invalid block；
- victim zone；
- moved live block；
- restart recovery。
