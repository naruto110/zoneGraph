#include "utils.h"
#include <string.h>
#include "pma.h"
#include "benchmark.h"
#include <cerrno>
#include <climits>
#include <iostream>
#include <cstdlib>
// #include <chrono> // test git
// #include "global.h"
using namespace std;

// 数据集路径表。
// 运行方式：./main <dataset_id> <thread_num> [mode] [--print-pma]
// dataset_id 对应下面的下标。
// mode 可选：baseline、new、parallel、parallel2，默认 parallel2。
vector<string> data_path;
void set_path(vector<string>& data_path)  
{
    data_path.clear();
	data_path.push_back("./data.txt");  // 0: 示例数据
	data_path.push_back("./validation_edges.txt");  // 1: 更新/删除验证数据
	data_path.push_back("/home/yingsui/data/livejournal/out.e");  // 2: LiveJournal
    data_path.push_back("/home/yingsui/data/orkut/out.e");  // 3: Orkut
    data_path.push_back("/home/yingsui/data/skitter/out.e");  // 4: Skitter
    data_path.push_back("/home/yingsui/data/youtube/out.e");  // 5: YouTube
    data_path.push_back("/home/yingsui/data/dblp/out.e");  // 6: DBLP
    data_path.push_back("/home/yingsui/data/roadCA/out.e");  // 7: RoadCA

    data_path.push_back("/home/yingsui/data/synthetic/syn19.e");  // 8: synthetic-19
    data_path.push_back("/home/yingsui/data/synthetic/syn20.e");  // 9: synthetic-20
    data_path.push_back("/home/yingsui/data/synthetic/syn21.e");  // 10: synthetic-21
	    data_path.push_back("/home/yingsui/data/synthetic/syn22.e");  // 11: synthetic-22
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/ca-GrQc.txt");  // 12: SNAP ca-GrQc
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/ca-CondMat.txt");  // 13: SNAP ca-CondMat
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/email-Enron.txt");  // 14: SNAP email-Enron
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/com-Amazon.txt");  // 15: SNAP com-Amazon
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/com-DBLP.txt");  // 16: SNAP com-DBLP
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/com-Youtube.txt");  // 17: SNAP com-Youtube
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/twitch-gamers.txt");  // 18: SNAP twitch-gamers
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/as-Skitter.txt");  // 19: SNAP as-Skitter
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/soc-Pokec.txt");  // 20: SNAP soc-Pokec
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/soc-LiveJournal1.txt");  // 21: SNAP soc-LiveJournal1
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/com-Orkut.txt");  // 22: SNAP com-Orkut
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
    cerr << "Usage: " << program << " <dataset_id> <thread_num> [mode] [options]" << endl;
    cerr << "  dataset_id: 0-" << (data_path.size() - 1) << endl;
    cerr << "  mode: baseline | new | parallel | parallel2 (default: parallel2)" << endl;
    cerr << "  --print-pma: print each vertex cache/PMA block after insertion" << endl;
    cerr << "  --shuffle-edges: shuffle loaded edge order with a fixed seed before graph maintenance" << endl;
    cerr << "  --insert <src> <dst>: insert an undirected edge after build" << endl;
    cerr << "  --delete <src> <dst>: delete an undirected edge after insertion" << endl;
    cerr << "  --compact-delete: compact PMA layout after each PMA deletion" << endl;
    cerr << "  --no-cache: disable cache_buf and route updates directly to PMA blocks" << endl;
    cerr << "  --cache-threshold <n>: cache up to n neighbors per vertex before PMA migration (0-"
         << CACHE_BLK_SIZE << ", default " << CACHE_BLK_SIZE << ")" << endl;
    cerr << "  --cache-block-threshold <n>: flush dirty PMA blocks after n cached blocks (default "
         << cache_blk_threshold << ")" << endl;
    cerr << "  --pending-merge-threshold <n>: merge pending PMA edges after n queued edges (default "
         << pending_merge_threshold << ")" << endl;
    cerr << "  --debug-pending-merge: print block values after pending merge/rebalance" << endl;
    cerr << "  --zone-writer: use zone-serial writer when flushing queued/parallel updates" << endl;
    cerr << "  --zone-class: place hot/cold/hub vertices into separate zone pools" << endl;
    cerr << "  --basic-gc: choose the first zone with invalid blocks as victim" << endl;
    cerr << "  --atz-gc: enable Asymmetric Temperature-aware Zone GC" << endl;
    cerr << "  --max-data-zones <n>: limit data zones used for allocation in GC experiments" << endl;
    cerr << "  --gc-proactive-free-blocks <n>: run GC after dump when free data blocks fall below n" << endl;
    cerr << "  --zns-record: write self-describing ZNS records (default; kept for compatibility)" << endl;
    cerr << "  --redo-log: enable packed edge redo logging for recovery validation" << endl;
    cerr << "  --recovery-test: after build, clear DRAM metadata and recover from ZNS+redo" << endl;
    cerr << "  --recover-zns: recover graph by scanning ZNS records; does not reset zones" << endl;
    cerr << "  --final-dump: after load, flush dirty PMA blocks to ZNS (cache-only vertices are not persisted)" << endl;
    cerr << "  --print-dram-metadata: estimate DRAM graph metadata after load, excluding PMA block payload" << endl;
}

int main(int argc, char* argv[]) {
    // cout << "int32_t size: " << sizeof(int32_t) << endl;
    // cout << "int size: " << sizeof(int) << endl;
    
    set_path(data_path);
    if (argc < 3) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    int dataset_id = 0;
    int thread_num = 0;
    if (!parse_int_arg(argv[1], dataset_id) ||
        dataset_id < 0 ||
        dataset_id >= static_cast<int>(data_path.size())) {
        cerr << "Invalid dataset_id: " << argv[1] << endl;
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    if (!parse_int_arg(argv[2], thread_num) || thread_num <= 0) {
        cerr << "Invalid thread_num: " << argv[2] << endl;
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    string mode = "parallel2";
    bool print_pma = false;
    bool shuffle_edges = false;
    bool mode_set = false;
    bool use_zone_writer = false;
    bool recover_zns = false;
    bool recovery_test = false;
    bool skip_final_dump = true;
    bool print_dram_metadata = false;
    vector<pair<int32_t, int32_t>> insert_edges;
    vector<pair<int32_t, int32_t>> delete_edges;
    for (int argi = 3; argi < argc; argi++) {
        string arg = argv[argi];
        if (arg == "--print-pma") {
            print_pma = true;
        } else if (arg == "--shuffle-edges") {
            shuffle_edges = true;
        } else if (arg == "--print-dram-metadata") {
            print_dram_metadata = true;
        } else if (arg == "--debug-pending-merge") {
            global_debug_pending_merge = true;
        } else if (arg == "--compact-delete") {
            global_delete_compaction_enabled = true;
        } else if (arg == "--no-cache") {
            global_cache_buf_threshold = 0;
        } else if (arg == "--cache-threshold") {
            if (argi + 1 >= argc) {
                cerr << "--cache-threshold requires <n>" << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            int threshold = 0;
            if (!parse_int_arg(argv[argi + 1], threshold) ||
                threshold < 0 || threshold > CACHE_BLK_SIZE) {
                cerr << "Invalid cache threshold: " << argv[argi + 1] << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            global_cache_buf_threshold = threshold;
            argi += 1;
        } else if (arg == "--cache-block-threshold") {
            if (argi + 1 >= argc) {
                cerr << "--cache-block-threshold requires <n>" << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            int threshold = 0;
            if (!parse_int_arg(argv[argi + 1], threshold) || threshold <= 0) {
                cerr << "Invalid cache block threshold: " << argv[argi + 1] << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            cache_blk_threshold = threshold;
            argi += 1;
        } else if (arg == "--pending-merge-threshold") {
            if (argi + 1 >= argc) {
                cerr << "--pending-merge-threshold requires <n>" << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            int threshold = 0;
            if (!parse_int_arg(argv[argi + 1], threshold) || threshold <= 0) {
                cerr << "Invalid pending merge threshold: " << argv[argi + 1] << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            pending_merge_threshold = threshold;
            argi += 1;
        } else if (arg == "--zone-writer") {
            use_zone_writer = true;
        } else if (arg == "--zone-class") {
            global_zone_class_enabled = true;
        } else if (arg == "--basic-gc") {
            global_basic_gc_enabled = true;
            global_atz_gc_enabled = false;
        } else if (arg == "--atz-gc") {
            global_atz_gc_enabled = true;
            global_basic_gc_enabled = false;
            global_zone_class_enabled = true;
        } else if (arg == "--max-data-zones") {
            if (argi + 1 >= argc) {
                cerr << "--max-data-zones requires <n>" << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            int max_data_zones = 0;
            if (!parse_int_arg(argv[argi + 1], max_data_zones) || max_data_zones <= 0) {
                cerr << "Invalid max data zones: " << argv[argi + 1] << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            global_max_data_zones_override = max_data_zones;
            argi += 1;
        } else if (arg == "--gc-proactive-free-blocks") {
            if (argi + 1 >= argc) {
                cerr << "--gc-proactive-free-blocks requires <n>" << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            int threshold = 0;
            if (!parse_int_arg(argv[argi + 1], threshold) || threshold < 0) {
                cerr << "Invalid GC proactive free block threshold: " << argv[argi + 1] << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            global_gc_proactive_free_blocks = threshold;
            argi += 1;
        } else if (arg == "--zns-record") {
            global_zns_record_enabled = true;
        } else if (arg == "--redo-log") {
            global_zns_record_enabled = true;
            global_redo_log_enabled = true;
        } else if (arg == "--recovery-test") {
            global_zns_record_enabled = true;
            global_redo_log_enabled = true;
            recovery_test = true;
        } else if (arg == "--recover-zns") {
            global_zns_record_enabled = true;
            global_redo_log_enabled = true;
            recover_zns = true;
        } else if (arg == "--skip-final-dump") {
            skip_final_dump = true;
        } else if (arg == "--final-dump") {
            skip_final_dump = false;
        } else if (arg == "--insert") {
            if (argi + 2 >= argc) {
                cerr << "--insert requires <src> <dst>" << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            int insert_src = 0;
            int insert_dst = 0;
            if (!parse_int_arg(argv[argi + 1], insert_src) ||
                !parse_int_arg(argv[argi + 2], insert_dst)) {
                cerr << "Invalid insert edge: " << argv[argi + 1]
                     << " " << argv[argi + 2] << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            insert_edges.push_back({insert_src, insert_dst});
            argi += 2;
        } else if (arg == "--delete") {
            if (argi + 2 >= argc) {
                cerr << "--delete requires <src> <dst>" << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            int delete_src = 0;
            int delete_dst = 0;
            if (!parse_int_arg(argv[argi + 1], delete_src) ||
                !parse_int_arg(argv[argi + 2], delete_dst)) {
                cerr << "Invalid delete edge: " << argv[argi + 1]
                     << " " << argv[argi + 2] << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            delete_edges.push_back({delete_src, delete_dst});
            argi += 2;
        } else {
            if (mode_set) {
                cerr << "Invalid option: " << arg << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            mode = arg;
            mode_set = true;
        }
    }
    string dataPath = data_path[dataset_id];
    THREAD_NUM = thread_num;
    // 初始化全局 ZNS 设备状态

    memset(&zns, 0, sizeof(zns));
    zns.path = NVME_DEVICE_PATH;
    zns.dev_fd = -1;
    zns.direct_fd = -1;

    // 检测目标设备是否为 zoned block device
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
    // cout << "test 0: " << zns.zones_info[0].blk_info_vec[0].src << endl;
    report_device_info(zns.info);
    // report_zones_info(zns.zones, zns.nr_zones);
    if (!recover_zns && zbd_reset_zones(zns.dev_fd, 0, 0) != 0) {
        perror("zbd_reset_zones");
        return EXIT_FAILURE;
    }
    if (!recover_zns) {
        reset_zns_write_pointers();
        resetRedoLogFile();
    }

    if (recover_zns) {
        recover_graph_from_zns_records();
        if (print_pma) {
            print_all_vertex_pma_blocks();
        }
        return EXIT_SUCCESS;
    }

    bool stream_load = (mode == "parallel2" && dataset_id == 22);
    if (stream_load) {
        if (scan_graph_metadata(dataPath) != 0) {
            return EXIT_FAILURE;
        }
    } else if (load_graph(dataPath) != 0) {
        return EXIT_FAILURE;
    }
    cout << "Total number of edges: " << global_edge_num << " global_max_vid: " << global_max_vid << endl;
    cout << "global_num_queues: " << global_num_queues << endl;

    if (shuffle_edges) {
        if (global_edge_vec.empty()) {
            cerr << "--shuffle-edges is ignored because this mode uses streaming load." << endl;
        } else {
            std::mt19937 g(1);
            shuffle(global_edge_vec.begin(), global_edge_vec.end(), g);
            cout << "Shuffled edge order with seed: 1" << endl;
        }
    }

    // graph_init(global_max_vid);
    // graph_init_new(global_max_vid);

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

    if (mode == "baseline") {
        graph_maintenance_baseline();
    } else if (mode == "new") {
        graph_maintenance_new();
    } else if (mode == "parallel") {
        graph_maintenance_parallel();
    } else if (mode == "parallel2") {
        if (stream_load) {
            graph_maintenance_parallel2_streaming(dataPath);
        } else {
            graph_maintenance_parallel2();
        }
    } else {
        cerr << "Invalid mode: " << mode << endl;
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    for (const auto &edge : insert_edges) {
        bool inserted = graph_insert_undirected_edge(edge.first, edge.second);
        cout << "Insert edge (" << edge.first << ", " << edge.second << "): "
             << (inserted ? "inserted" : "already exists") << endl;
    }
    for (const auto &edge : delete_edges) {
        bool deleted = graph_delete_undirected_edge(edge.first, edge.second);
        cout << "Delete edge (" << edge.first << ", " << edge.second << "): "
             << (deleted ? "deleted" : "not found") << endl;
    }
    merge_all_pending_insert_buffers();
    if (print_dram_metadata) {
        print_dram_metadata_usage();
    }
    if (!skip_final_dump && global_cache_blk > 0) {
        if (use_zone_writer) {
            dumpToZoneZoneSerial();
        } else if (mode == "parallel2") {
            dumpToZoneParallel();
        } else {
            dumpToZone();
        }
        global_cache_blk = 0;
    }
    flushZNSRedoLog();
    closeRedoLogFile();

    if (recovery_test) {
        run_recovery_experiment_from_zns();
    }

    if (print_pma) {
        print_all_vertex_pma_blocks();
    }

    if (!global_zns_record_enabled) {
        graph_recovery();
    }
    // global_edge_vec.clear();
    // global_edge_vec.shrink_to_fit();

    // // 测试垃圾回收
    // // garbage_collection_new();
    // int64_t total_remain_blk = 0;
    // for (int i = 0; i < global_blk_queues.size(); i++) {
    //     total_remain_blk += global_blk_queues[i]->memoryWPs.size();
    // }
    // cout << "global total blks: " << (zns.nr_blocks_perzone * global_blk_queues.size()) << " total_remain_blk: " << total_remain_blk << endl;

    // while (1);

    // graph_maintenance_new();
    // garbage_collection();

    // cout << "test-- print" << endl;
    // print_graph();

    // test_benchmark();

    // int dataToWrite[] = {11, 22, 32, 24, 25};
    // writeBlockToZone(zns.dev_fd, dataToWrite, 5, 201326592);
    // zbd_reset_zones(zns.dev_fd, 201326592, zns.zone_size);

    return EXIT_SUCCESS;
}
