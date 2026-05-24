# zoneGraph pending merge 小规模验证

`validation_edges.txt` 是一个小规模无向图输入，用于观察：

- cache-only 顶点超过阈值后迁入 PMA；
- 后续边先进入 pending；
- pending 达到阈值后 merge 到 PMA；
- pending 批量写入导致目标 block 空位不够时，触发 PMA 局部 rebalance 或 resize；
- dirty PMA block 刷盘后仍能读出同样的 PMA 内容。

当前数据集包含：

```text
0 1
0 2
...
0 114
```

默认 `CACHE_BLK_SIZE=16`，所以顶点 0 会先迁入 PMA。当前默认：

```cpp
pending_merge_threshold = 2
```

前面的 pending 批次会逐步填充 PMA block；当某个 block 只剩 1 个空 slot，而 pending 里有 2 条落到同一 block 的边时，调试输出会显示 merge 前后的状态。

## 建议命令

推荐用串行 `new` 模式，输出顺序最容易读：

```bash
sudo ./main_pma_check_new 1 1 new --debug-pending-merge
```

如果想显式指定阈值：

```bash
sudo ./main_pma_check_new 1 1 new \
  --pending-merge-threshold 2 \
  --debug-pending-merge
```

重点看三段标题：

```text
==== before threshold pending merge ...
==== after threshold pending merge ...
==== after dirty block flush ...
```

第一段能看到 PMA 已有值加 pending 值的 `effective_neighbors`，以及 pending 落在哪个 block；第二段能看到 pending 被清空、PMA block 经 rebalance/resize 后的结果；第三段表示 dirty block 已写到 ZNS，再从 PMA/ZNS 路径读出的结果。
