#include "utils.h"
#include "pma.h"
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

static vector<string> data_path;

static void set_path(vector<string> &paths)
{
    paths.clear();
    paths.push_back("./data.txt");
    paths.push_back("./validation_edges.txt");
    paths.push_back("/home/yingsui/data/livejournal/out.e");
    paths.push_back("/home/yingsui/data/orkut/out.e");
    paths.push_back("/home/yingsui/data/skitter/out.e");
    paths.push_back("/home/yingsui/data/youtube/out.e");
    paths.push_back("/home/yingsui/data/dblp/out.e");
    paths.push_back("/home/yingsui/data/roadCA/out.e");
    paths.push_back("/home/yingsui/data/synthetic/syn19.e");
    paths.push_back("/home/yingsui/data/synthetic/syn20.e");
    paths.push_back("/home/yingsui/data/synthetic/syn21.e");
    paths.push_back("/home/yingsui/data/synthetic/syn22.e");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/ca-GrQc.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/ca-CondMat.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/email-Enron.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/com-Amazon.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/com-DBLP.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/com-Youtube.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/twitch-gamers.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/as-Skitter.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/soc-Pokec.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/soc-LiveJournal1.txt");
    paths.push_back("/home/femu/project/datasets/snap/edgelist/com-Orkut.txt");
}

static bool parse_int_arg(const char *arg, int &value)
{
    errno = 0;
    char *end = nullptr;
    long parsed = strtol(arg, &end, 10);
    if (errno != 0 || end == arg || *end != '\0' || parsed < 0 || parsed > INT_MAX) {
        return false;
    }
    value = static_cast<int>(parsed);
    return true;
}

static void print_usage(const char *program)
{
    cerr << "Usage: " << program
         << " <dataset_id> <thread_num> [hot_vertices] [updates_per_vertex] [options]" << endl;
    cerr << "  --basic-gc" << endl;
    cerr << "  --atz-gc" << endl;
    cerr << "  --max-data-zones <n>" << endl;
    cerr << "  --cache-block-threshold <n>" << endl;
    cerr << "  --pending-merge-threshold <n>" << endl;
    cerr << "  --gc-proactive-free-blocks <n>" << endl;
}

static void reset_stats(void)
{
    total_dump_time = 0.0;
    total_dump_count = 0;
    total_gc_time = 0.0;
    total_gc_count = 0;
    total_gc_moved_blocks = 0;
    atz_hot_zone_gc_count = 0;
    atz_warm_zone_gc_count = 0;
    atz_mixed_zone_gc_count = 0;
    atz_cold_zone_gc_count = 0;
    atz_emergency_gc_count = 0;
    atz_migrated_hot_blocks = 0;
    atz_migrated_warm_blocks = 0;
    atz_migrated_cold_blocks = 0;
    atz_cold_zone_skipped_count = 0;
    atz_cold_zone_protected_count = 0;
    atz_invalid_blocks_reclaimed = 0;
    atz_gc_write_bytes = 0;
    atz_user_write_bytes = 0;
    atz_hot_block_update_writes = 0;
    atz_warm_block_update_writes = 0;
    atz_cold_block_update_writes = 0;
}

static vector<int32_t> pick_hot_vertices(int32_t wanted)
{
    vector<pair<int32_t, int32_t>> candidates;
    candidates.reserve(global_vectex_vec.size());
    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr) {
            continue;
        }
        if (v->degree >= 64 && v->degree < hub_degree_threshold && !v->blk_list.empty()) {
            candidates.push_back({v->degree, v->id});
        }
    }
    sort(candidates.begin(), candidates.end(), greater<pair<int32_t, int32_t>>());
    vector<int32_t> result;
    for (int32_t i = 0; i < wanted && i < static_cast<int32_t>(candidates.size()); i++) {
        result.push_back(candidates[i].second);
    }
    return result;
}

int main(int argc, char *argv[])
{
    set_path(data_path);
    if (argc < 3) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    int dataset_id = 0;
    int thread_num = 0;
    if (!parse_int_arg(argv[1], dataset_id) ||
        dataset_id < 0 ||
        dataset_id >= static_cast<int>(data_path.size()) ||
        !parse_int_arg(argv[2], thread_num) ||
        thread_num <= 0) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    int hot_vertex_count = 2048;
    int updates_per_vertex = 512;
    if (argc >= 4 && argv[3][0] != '-' &&
        (!parse_int_arg(argv[3], hot_vertex_count) || hot_vertex_count <= 0)) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    if (argc >= 5 && argv[4][0] != '-' &&
        (!parse_int_arg(argv[4], updates_per_vertex) || updates_per_vertex <= 0)) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    int argi = 3;
    if (argi < argc && argv[argi][0] != '-') {
        argi++;
    }
    if (argi < argc && argv[argi][0] != '-') {
        argi++;
    }

    THREAD_NUM = thread_num;
    global_max_data_zones_override = 8;
    cache_blk_threshold = 5000;
    pending_merge_threshold = 500000;
    global_gc_proactive_free_blocks = 15000000;

    for (; argi < argc; argi++) {
        string arg = argv[argi];
        if (arg == "--basic-gc") {
            global_basic_gc_enabled = true;
            global_atz_gc_enabled = false;
        } else if (arg == "--atz-gc") {
            global_atz_gc_enabled = true;
            global_basic_gc_enabled = false;
            global_zone_class_enabled = true;
        } else if (arg == "--max-data-zones" ||
                   arg == "--cache-block-threshold" ||
                   arg == "--pending-merge-threshold" ||
                   arg == "--gc-proactive-free-blocks") {
            if (argi + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            int value = 0;
            if (!parse_int_arg(argv[argi + 1], value) || value <= 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            if (arg == "--max-data-zones") {
                global_max_data_zones_override = value;
            } else if (arg == "--cache-block-threshold") {
                cache_blk_threshold = value;
            } else if (arg == "--pending-merge-threshold") {
                pending_merge_threshold = value;
            } else {
                global_gc_proactive_free_blocks = value;
            }
            argi++;
        } else {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    memset(&zns, 0, sizeof(zns));
    zns.path = NVME_DEVICE_PATH;
    zns.dev_fd = -1;
    zns.direct_fd = -1;

    if (!zbd_device_is_zoned(zns.path)) {
        cerr << zns.path << " is not a zoned block device" << endl;
        return EXIT_FAILURE;
    }
    cout << "This is a ZBD" << endl;
    if (zns_open() != 0) {
        cerr << "Failed to open ZNS device: " << zns.path << endl;
        return EXIT_FAILURE;
    }
    if (global_max_data_zones_override > 0 &&
        global_max_data_zones_override < global_num_queues) {
        global_num_queues = global_max_data_zones_override;
    }
    report_device_info(zns.info);
    if (zbd_reset_zones(zns.dev_fd, 0, 0) != 0) {
        perror("zbd_reset_zones");
        return EXIT_FAILURE;
    }
    reset_zns_write_pointers();

    string data_path_for_run = data_path[dataset_id];
    if (load_graph(data_path_for_run) != 0) {
        return EXIT_FAILURE;
    }
    cout << "Dataset id: " << dataset_id
         << " directed edges: " << global_edge_num
         << " max vid: " << global_max_vid
         << " data zones: " << global_num_queues
         << " gc mode: "
         << (global_atz_gc_enabled ? "atz-gc" :
             (global_basic_gc_enabled ? "basic" : "greedy"))
         << " cache_block_threshold: " << cache_blk_threshold
         << " proactive_free_blocks: " << global_gc_proactive_free_blocks << endl;

    reset_stats();
    graph_maintenance_parallel2();
    merge_all_pending_insert_buffers();
    if (global_cache_blk > 0) {
        dumpToZoneParallel();
    }

    global_edge_vec.clear();
    global_edge_vec.shrink_to_fit();
    vector<int32_t> hot_vertices = pick_hot_vertices(hot_vertex_count);
    cout << "Selected hot vertices: " << hot_vertices.size()
         << " requested: " << hot_vertex_count
         << " updates_per_vertex: " << updates_per_vertex << endl;
    if (hot_vertices.empty()) {
        cerr << "No hot candidates found." << endl;
        return EXIT_FAILURE;
    }

    reset_stats();
    double update_start = get_current_time();
    int64_t update_count = 0;
    const int64_t merge_interval = 200000;
    int32_t synthetic_dst_base = 1000000000;

    for (int round = 0; round < updates_per_vertex; round++) {
        for (int32_t src : hot_vertices) {
            int32_t dst = synthetic_dst_base - static_cast<int32_t>(update_count % 500000000);
            graph_insert_directed_edge(src, dst);
            update_count++;
            if (update_count % merge_interval == 0) {
                merge_all_pending_insert_buffers();
                if (global_cache_blk > 0) {
                    dumpToZoneParallel();
                }
            }
        }
    }
    merge_all_pending_insert_buffers();
    if (global_cache_blk > 0) {
        dumpToZoneParallel();
    }
    double update_total = get_current_time() - update_start;

    cout << "Skew update total directed inserts: " << update_count << endl;
    cout << "Skew update total time: " << update_total << endl;
    cout << "Skew update dump count: " << total_dump_count << endl;
    cout << "Skew update dump total time: " << total_dump_time << endl;
    cout << "Skew update GC count: " << total_gc_count << endl;
    cout << "Skew update GC total time: " << total_gc_time << endl;
    cout << "Skew update GC moved blocks: " << total_gc_moved_blocks << endl;
    print_atz_gc_summary();
    double foreground_time = update_total - total_dump_time - total_gc_time;
    if (foreground_time < 0.0) {
        foreground_time = 0.0;
    }
    cout << "Skew update foreground time: " << foreground_time << endl;
    zns_close();
    return EXIT_SUCCESS;
}
