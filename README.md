# zoneGraph

zoneGraph is a research prototype for maintaining dynamic graphs directly on NVMe Zoned Namespace (ZNS) devices. The system combines an in-memory cache for low-degree vertices with a Packed Memory Array (PMA) stored in zoned flash, and supports parallel ingestion, recovery, and analytics such as PageRank and BFS.

## Highlights
- Zoned storage aware allocator that streams adjacency blocks into sequential zone writes while tracking block metadata in `ZoneStatus`.
- Vertex representation that keeps hot neighbors in an in-DRAM cache (`CACHE_BLK_SIZE`) and spills cold data to PMA-managed blocks (`pm_blk`).
- Parallel ingestion pipelines (`graph_maintenance_parallel*`) that batch edges, flush blocks to ZNS, and run background garbage collection.
- Recovery path (`graph_recovery`) that reloads vertex state from the device without rebuilding the entire graph.
- GAP Benchmark Suite integration (`gapbs/`) providing PageRank, BFS, and supporting utility code.

## Requirements
- Linux with NVMe Zoned Namespace (ZNS) support and a device exposed through `libzbd`.
- C++17-capable compiler with OpenMP support (tested with `g++` 11+).
- Development headers for `libzbd` (e.g., `sudo apt install libzbd-dev`).
- Optional: access to real graph datasets in plain edge-list format.

## Build
Compile the driver with your preferred compiler; the command below is a starting point:

```bash
g++ -std=c++17 -O3 -fopenmp main.cpp -o zoneGraph -lzbd
```

All headers are self-contained; no external build system is required. If you split the code into translation units or add new files, update the compilation command accordingly.

## Running
The executable expects two arguments:

```bash
./zoneGraph <dataset_id> <thread_count>
```

`dataset_id` indexes into the hard-coded path table in `main.cpp::set_path`. Update this table to point to your datasets (edge lists with space-separated `src dst` pairs). `thread_count` controls OpenMP parallelism (`THREAD_NUM`).

When running against a real device:
- Ensure `NVME_DEVICE_PATH` in `global.h` matches the zoned block device (default is `/dev/nullb0` for local testing).
- The program resets all zones before loading data. Use a dedicated test namespace to avoid data loss.

## Data Format
`load_graph()` consumes text files containing integer vertex IDs per line (`src des`). Edges are stored bidirectionally by default; undirected graphs should list each edge once.

## Repository Layout
- `main.cpp` – CLI entry point, device bootstrap, ingestion flow control.
- `global.h` – Global configuration, type definitions, and shared state.
- `utils.h` – Zoned device helpers, garbage collection, and block allocators.
- `pma.h` – PMA implementation, graph ingestion logic, recovery, and analytics kernels.
- `benchmark.h` – PageRank, BFS, and other GAP-inspired routines.
- `gapbs/` – Third-party headers from the GAP Benchmark Suite.

## Notes & Next Steps
- Garbage collection assumes one spare zone (`global_empty_zone`). Adjust thresholds (`ZONE_GC_THRESHOLD`) to tune reclamation behavior.
- The prototype uses `malloc/new` for buffers; consider integrating a custom allocator if you extend the system.
- Thread safety relies on `std::shared_timed_mutex` per vertex. Additional synchronization may be needed for new mutation paths.
- Add your own benchmarks under `gapbs/test/` or extend `benchmark.h` to evaluate new algorithms.

