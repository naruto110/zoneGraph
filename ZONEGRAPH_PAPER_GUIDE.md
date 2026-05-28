# zoneGraph Paper Writing Guide

This guide focuses on the three mechanisms that should anchor the paper: graph maintenance, garbage collection, and recovery.

## Core Positioning

zoneGraph is a dynamic graph maintenance system for ZNS SSDs. The key design point is to absorb fine-grained edge updates in DRAM, maintain per-vertex adjacency with PMA-style sorted blocks, and write ZNS in device-friendly block/zone granularity.

The paper should emphasize this contrast:

```text
edge-granular logical updates
  -> DRAM pending/cache maintenance
  -> block-granular ZNS persistence
  -> zone-aware GC and recovery
```

Do not present zoneGraph as issuing one physical SSD write per inserted edge. In the performance path, redo logging is disabled. In the recovery path, logical redo entries are packed into 512B redo pages before being written to file/ZNS redo space.

## Graph Maintenance

### Problem

Dynamic graph updates are naturally edge-granular, but ZNS SSDs require sequential zone writes and strongly prefer coarse writes. Directly writing each inserted edge creates write amplification, random-write pressure, and poor zone utilization.

### Design

Each vertex owns a PMA-like adjacency layout. Neighbor ids are stored in sorted sparse blocks, which preserve query-friendly ordered adjacency while leaving slack for future inserts.

The maintenance path has three levels:

1. Per-edge insertion receives directed edges one at a time.
2. New edges first enter per-vertex pending buffers grouped by target PMA block.
3. Pending values are merged in batches into PMA blocks; if local free slots are insufficient, zoneGraph rebalances a local PMA window or rebuilds/resizes the vertex PMA.

Relevant implementation names:

```text
graph_insert_directed_edge()
pending_insert_by_block
pending_insert_unplaced
maybe_merge_pending_insert_buffers()
merge_pending_insert_buffer()
graph_maintenance_parallel2()
dumpToZoneParallel()
```

### What To Claim

The graph maintenance contribution is not just parallel insertion. It is the separation of update granularity and persistence granularity:

```text
fine-grained edge arrival
batched DRAM merge
block-aligned persistence
```

This lets zoneGraph support streaming updates while avoiding the write amplification of immediate per-edge persistence.

### Evaluation Metrics

Report:

```text
total maintenance time
foreground compute time
dump count
dump time
GC count
GC time
throughput = directed_edges / total_time
```

When `dump count = 0`, clearly state that the run measures DRAM maintenance under the chosen cache threshold, not ZNS dump behavior.

## Garbage Collection

### Problem

ZNS zones cannot be overwritten in place. When graph blocks are updated, old physical block versions become invalid. Over time, zones contain a mix of live and stale graph blocks, so free space must be reclaimed by moving valid blocks and resetting zones.

### Design

zoneGraph tracks the logical mapping:

```text
(vertex id, PMA block id) -> physical zone write pointer
```

On dump, new block versions are appended to zones and the mapping is updated. Old versions become reclaimable. GC selects zones, copies still-live blocks elsewhere, updates mappings, then resets the reclaimed zone.

Relevant implementation names:

```text
mark_zone()
dumpToZoneZoneSerial()
maybe_run_proactive_gc()
garbage_collection_greedy()
requestRelocationBlock()
global_gc_proactive_free_blocks
```

### Writing Angle

The key point is that GC is graph-aware at block granularity. zoneGraph does not need to understand individual edges during GC. It only relocates valid PMA blocks and preserves the logical block mapping.

If discussing ATZ-GC, frame it as an optimization over basic greedy GC:

```text
basic GC: reclaim zones based on valid-block pressure
ATZ-GC: separate hot/warm/cold graph blocks to reduce future relocation
```

Use ATZ-GC only if the experiment section includes clear results for it.

### Evaluation Metrics

Report:

```text
GC count
GC total time
GC moved blocks
write amplification from relocated blocks
foreground slowdown with/without proactive GC
```

## Recovery

### Problem

After a crash, DRAM metadata is lost:

```text
vertex objects
PMA shape metadata
pending/cache state
logical-to-physical block mappings
recent updates not yet dumped
```

ZNS still contains block records, but recovery must identify committed records, rebuild latest mappings, and replay recent edge updates.

### Design

zoneGraph recovery uses two persistent sources:

1. ZNS data records for dumped PMA blocks.
2. Redo log pages for recent edge updates.

The ZNS data record contains enough metadata to reconstruct a PMA block:

```text
src vertex id
PMA block id
valid count
PMA element count
PMA capacity
epoch/sequence/checksum metadata
```

The redo log records logical edge insertions. Redo entries are accumulated in memory and flushed as 512B redo pages, not written as one physical redo write per directed edge.

Relevant implementation names:

```text
appendZNSRedoLog()
flushZNSRedoLog()
decodeZNSRedoLogPageRaw()
recover_graph_from_zns_records()
run_recovery_experiment_from_zns()
replay_redo_logs_after()
clear_dram_graph_metadata_for_recovery()
```

### Recovery Flow

Describe the flow as:

```text
scan metadata zone for committed epochs
scan data zones for valid block records
select latest committed record per logical PMA block
rebuild vertex/PMA metadata and logical-to-physical mapping
materialize missing PMA blocks if needed
read redo pages
replay idempotent edge records
validate recovered graph sample
```

The redo replay should be described as idempotent. Replaying an edge already present should not corrupt the graph because insertion checks/merge logic preserve unique adjacency entries.

### What Not To Claim

Avoid claiming that every directed edge causes an immediate redo write to SSD. The logical API appends a redo entry, but physical redo persistence is page-packed.

Avoid mixing recovery overhead into normal maintenance results unless the run explicitly enables:

```bash
--redo-log
--recovery-test
--recover-zns
```

### Evaluation Metrics

Report:

```text
recovery_time_sec
metadata slots scanned
data slots scanned
redo slots scanned
valid committed epochs
rebuilt vertices
rebuilt blocks
redo records replayed
recovery_vs_load
validation PASS/FAIL
```

## Suggested Paper Structure

1. Introduction
   - Dynamic graph maintenance conflicts with ZNS sequential-write constraints.
   - zoneGraph decouples edge arrival from ZNS persistence.

2. Background and Motivation
   - Dynamic graph update behavior.
   - PMA adjacency layout.
   - ZNS append/reset constraints.

3. System Overview
   - Edge ingestion, pending buffers, PMA blocks, ZNS zones.
   - Logical-to-physical block mapping.

4. Graph Maintenance
   - Per-edge directed insertion.
   - Pending buffer batching.
   - PMA local merge, rebalance, resize.
   - Dump threshold and block-aligned persistence.

5. Garbage Collection
   - Invalid block creation.
   - Zone selection and relocation.
   - Optional hot/warm/cold separation.

6. Recovery
   - ZNS data records.
   - Commit metadata.
   - 512B packed redo pages.
   - Rebuild and replay.

7. Evaluation
   - Load/maintenance throughput against LSMGraph and MBFGraph.
   - Cache threshold sensitivity.
   - Dump/GC overhead.
   - Recovery time and validation.

## Recommended Terminology

Use these terms consistently:

```text
directed edge
pending insert buffer
PMA block
dirty/cache block
block-aligned dump
logical PMA block
physical ZNS record
logical-to-physical mapping
committed epoch
packed redo page
idempotent replay
```

