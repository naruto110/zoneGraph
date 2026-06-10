# zoneGraph

zoneGraph is a prototype dynamic graph storage engine for Zoned Namespace (ZNS) SSDs. It is designed around the mismatch between fine-grained graph mutations and the append-only, block/zone-granular persistence model of ZNS devices.

At a high level, zoneGraph turns edge-level insertions and deletions into ZNS-friendly logical block versions:

```text
edge-level graph updates
  -> degree-aware DRAM adjacency buffers
  -> PMA-style logical adjacency blocks
  -> append-only self-describing ZNS records
```

The project accompanies the paper draft "zoneGraph: A ZNS-native Design for Dynamic Graphs".

## Core Ideas

### Hierarchical Buffered Block Array

zoneGraph stores adjacency data in a Hierarchical Buffered Block Array (HBBA):

- Low-degree vertices are kept in compact DRAM buffers to avoid sparse block materialization.
- High-degree vertices are promoted to PMA-style logical adjacency blocks.
- Per-block pending buffers coalesce repeated updates before PMA maintenance.
- Dirty logical blocks are flushed to ZNS as append-only records.

HBBA preserves macro-level ordering across each vertex's logical block sequence, while allowing batched local maintenance within modified blocks.

### ZNS Block Versions

Each persisted adjacency block is written as a self-describing ZNS record containing its logical identity and version metadata. The runtime mapping is:

```text
(vertex id, logical adjacency block id) -> latest physical ZNS record
```

When a logical block is rewritten, the new version is appended to a fresh ZNS location and the old physical version becomes invalid.

### Garbage Collection

zoneGraph reclaims obsolete block versions at zone granularity. Live adjacency blocks are migrated out of a victim zone, logical-to-physical mappings are updated, and the reclaimed zone is reset.

The prototype also includes ATZ-GC, an asymmetric temperature-aware GC policy that uses vertex update temperature to guide block placement and victim-zone selection, reducing repeated movement of cold graph data under skewed updates.

### Graph Reopen

Graph reopen reconstructs DRAM-resident metadata from persisted ZNS records. zoneGraph scans durable block versions, selects the latest record for each logical adjacency block, and rebuilds vertex metadata, PMA block descriptors, and logical-to-physical mappings without reloading the original edge stream.

This is a restart metadata reconstruction path, not a complete power-loss consistency protocol.

## Repository Layout

```text
zoneGraph/
  main.cpp                 Command-line driver
  pma.h                    HBBA/PMA maintenance, dump, GC, reopen logic
  utils.h                  ZNS record I/O and zone management helpers
  global.h                 Global configuration and data structures
  benchmark.h              BFS, PageRank, and WCC query benchmarks
  Makefile                 Build target for main_pma_check_new
  gc_skew_experiment.cpp   GC stress workload driver
  script/                  Plotting scripts for paper figures
  figures/                 Generated paper figures
  result/                  Local result summaries
```

## Requirements

- Linux with a zoned block device visible to `libzbd`
- `g++` with C++17 support
- OpenMP
- `libzbd` development headers and library
- Root permission for resetting and writing zones

On Ubuntu-like systems:

```bash
sudo apt-get install -y g++ make libzbd-dev
```

The default device path is compiled from `global.h`:

```c
#define NVME_DEVICE_PATH "/dev/nvme0n1"
```

Change this macro before building if your ZNS device uses a different path.

## Build

```bash
cd /home/femu/project/zoneGraph
make clean
make
```

The default binary is:

```text
main_pma_check_new
```

## Datasets

The command-line driver uses numeric dataset ids defined in `main.cpp`. Common SNAP datasets in this workspace are:

```text
16  com-DBLP
17  com-Youtube
18  twitch-gamers
19  as-Skitter
20  soc-Pokec
21  soc-LiveJournal1
```

The corresponding paths are under:

```text
/home/femu/project/datasets/snap/edgelist/
```

Input files are interpreted as undirected graphs by default:

```text
1 input edge -> 2 directed edges
```

Deletion workloads follow the same convention:

```text
1 undirected delete -> 2 directed delete attempts
```

## Basic Usage

General command format:

```bash
sudo ./main_pma_check_new <dataset_id> <thread_num> [mode] [options]
```

Modes:

```text
baseline   Insert directly into PMA-style blocks
new        Degree-aware single-threaded path
parallel2  Default optimized parallel path
```

Example: build and query `as-Skitter`:

```bash
cd /home/femu/project/zoneGraph

sudo ./main_pma_check_new 19 64 parallel2 \
  --cache-block-threshold 100000 \
  --pending-merge-threshold 4000000 \
  --run-queries \
  --final-dump
```

Useful options:

```text
--cache-threshold <n>           Low-degree cache threshold, default 16
--cache-block-threshold <n>     Flush dirty PMA blocks after n cached blocks
--pending-merge-threshold <n>   Merge pending PMA edges after n queued edges
--batch-size <n>                Directed edges per maintenance batch
--run-queries                   Run BFS, PageRank, and WCC after maintenance
--query-source <vid>            BFS source vertex, default 1
--pagerank-iters <n>            PageRank iterations, default 20
--final-dump                    Flush dirty PMA blocks after load
--delete-file <path>            Delete undirected edges listed in a text file
--delete-parallel               Execute delete-file with OpenMP workers
--basic-gc                      Use basic victim selection
--atz-gc                        Enable ATZ-GC
--max-data-zones <n>            Limit data zones for GC stress experiments
--gc-proactive-free-blocks <n>  Trigger proactive GC below free-block threshold
--recovery-test                 Clear DRAM metadata and reopen from ZNS records
--recover-zns                   Reopen by scanning existing ZNS records
```

Normal non-recovery runs reset the target ZNS device at startup.

## Graph Reopen Example

Build a graph, force a final dump, and run the reopen validation path:

```bash
sudo ./main_pma_check_new 19 64 parallel2 \
  --cache-block-threshold 100000 \
  --pending-merge-threshold 4000000 \
  --final-dump \
  --recovery-test
```

To reopen from existing ZNS records without resetting zones:

```bash
sudo ./main_pma_check_new 19 64 parallel2 --recover-zns
```

## Edge Deletion Example

```bash
sudo ./main_pma_check_new 19 64 parallel2 \
  --cache-block-threshold 100000 \
  --delete-file /path/to/delete_edges.txt \
  --delete-parallel
```

The delete file should contain one undirected edge per line:

```text
src dst
```

zoneGraph attempts both directed deletions for each input line.

## GC Stress Experiments

The repository includes a dedicated GC stress workload:

```bash
cd /home/femu/project/zoneGraph
./gc_skew_experiment --help
```

