# Experiment Commands

Default SSD experiments. zoneGraph resets zones internally on each normal run, so no cleanup is needed for zoneGraph.

## Dataset Names

SNAP text datasets are under:

```bash
/home/femu/project/datasets/snap/edgelist
```

Common datasets:

```bash
DATASET=as-Skitter
DATASET=com-Youtube
DATASET=soc-LiveJournal1
DATASET=com-Orkut
DATASET=soc-Pokec
```

zoneGraph dataset ids currently used:

```text
17  com-Youtube
19  as-Skitter
21  soc-LiveJournal1
```

## zoneGraph

LiveJournal, 8GB-memory-oriented DRAM cache:

```bash
cd /home/femu/project/zoneGraph

sudo ./main_pma_check_new 21 64 parallel2 \
  --cache-block-threshold 8000000 \
  --pending-merge-threshold 16000000
```

LiveJournal, smaller but not tiny thresholds:

```bash
cd /home/femu/project/zoneGraph

sudo ./main_pma_check_new 21 64 parallel2 \
  --cache-block-threshold 4000000 \
  --pending-merge-threshold 8000000
```

Skitter:

```bash
cd /home/femu/project/zoneGraph

sudo ./main_pma_check_new 19 64 parallel2 \
  --cache-block-threshold 200000 \
  --pending-merge-threshold 4000000
```

YouTube:

```bash
cd /home/femu/project/zoneGraph

sudo ./main_pma_check_new 17 64 parallel2 \
  --cache-block-threshold 200000 \
  --pending-merge-threshold 4000000
```

## LSMGraph

Template:

```bash
DATASET=soc-LiveJournal1

cd /home/femu/project/LSMGraph

./build/bin/bench_stream_load \
  -dataset_path=/home/femu/project/datasets/snap/edgelist/${DATASET}.txt \
  -db_path=/tmp/lsmgraph_default_${DATASET}_db \
  -mmap_path=/tmp/lsmgraph_default_${DATASET}_block.mmap \
  -support_mulversion=1 \
  -memtable_num=2

rm -rf /tmp/lsmgraph_default_${DATASET}_db \
  /tmp/lsmgraph_default_${DATASET}_block.mmap
```

The reported metric is the `STREAM_LOAD ... edge_time_sec=...` line.

## MBFGraph

Prebuilt ini files:

```bash
/home/femu/project/MBFGraph/prebuilt_ini/default_ssd
```

Template:

```bash
DATASET=soc-LiveJournal1

cd /home/femu/project/MBFGraph

cp prebuilt_ini/default_ssd/${DATASET}.ini \
  /tmp/mbfgraph_default_${DATASET}.ini

./bin/snap_stream_load \
  /home/femu/project/datasets/snap/edgelist/${DATASET}.txt \
  /tmp/mbfgraph_default_${DATASET} \
  --no-ini

cd /tmp

/home/femu/project/MBFGraph/bin/benchmark_driver \
  -g /tmp/mbfgraph_default_${DATASET} \
  -b degree_cnt \
  --physical_memory 1073741824 \
  --measure_scatter_gather \
  --heartbeat \
  -a \
  -p 64 \
  --page_size 2048 \
  --bits_per_BF 256 \
  --chip_v_once 4 \
  --max_check_v 1 \
  --kernel CPU \
  --bfs_per_page 256 \
  --bf_mode create

rm -f /tmp/mbfgraph_default_${DATASET} \
  /tmp/mbfgraph_default_${DATASET}.ini \
  /tmp/bf /tmp/vertices /tmp/updates0 /tmp/updates1
```

Use `snap_stream_load`'s `edge_time_sec` as append time. Use `benchmark_driver`'s `CORE::TIME::WALL` as Bloom-filter creation plus scan time. For end-to-end MBFGraph maintenance time, add the two.

