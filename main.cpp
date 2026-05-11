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
	data_path.push_back("/home/yingsui/data/livejournal/out.e");  // 1: LiveJournal
    data_path.push_back("/home/yingsui/data/orkut/out.e");  // 2: Orkut
    data_path.push_back("/home/yingsui/data/skitter/out.e");  // 3: Skitter
    data_path.push_back("/home/yingsui/data/youtube/out.e");  // 4: YouTube
    data_path.push_back("/home/yingsui/data/dblp/out.e");  // 5: DBLP
    data_path.push_back("/home/yingsui/data/roadCA/out.e");  // 6: RoadCA

    data_path.push_back("/home/yingsui/data/synthetic/syn19.e");  // 7: synthetic-19
    data_path.push_back("/home/yingsui/data/synthetic/syn20.e");  // 8: synthetic-20
    data_path.push_back("/home/yingsui/data/synthetic/syn21.e");  // 9: synthetic-21
	    data_path.push_back("/home/yingsui/data/synthetic/syn22.e");  // 10: synthetic-22
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
    cerr << "  --insert <src> <dst>: insert an undirected edge after build" << endl;
    cerr << "  --delete <src> <dst>: delete an undirected edge after insertion" << endl;
    cerr << "  --compact-delete: compact PMA layout after each PMA deletion" << endl;
    cerr << "  --zone-writer: use zone-serial writer when flushing queued/parallel updates" << endl;
    cerr << "  --zone-class: place hot/cold/hub vertices into separate zone pools" << endl;
    cerr << "  --gc-global: use original global GC strategy (default)" << endl;
    cerr << "  --gc-class-aware: choose victim globally, relocate by vertex class" << endl;
    cerr << "  --gc-class-local: GC requested class first, then fallback globally" << endl;
    cerr << "  --zns-record: write self-describing ZNS records for recovery experiments" << endl;
    cerr << "  --recover-zns: recover graph by scanning ZNS records; does not reset zones" << endl;
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
    bool mode_set = false;
    bool use_zone_writer = false;
    bool recover_zns = false;
    vector<pair<int32_t, int32_t>> insert_edges;
    vector<pair<int32_t, int32_t>> delete_edges;
    for (int argi = 3; argi < argc; argi++) {
        string arg = argv[argi];
        if (arg == "--print-pma") {
            print_pma = true;
        } else if (arg == "--compact-delete") {
            global_delete_compaction_enabled = true;
        } else if (arg == "--zone-writer") {
            use_zone_writer = true;
        } else if (arg == "--zone-class") {
            global_zone_class_enabled = true;
        } else if (arg == "--gc-global") {
            global_gc_mode = GC_GLOBAL;
        } else if (arg == "--gc-class-aware") {
            global_gc_mode = GC_CLASS_AWARE;
        } else if (arg == "--gc-class-local") {
            global_gc_mode = GC_CLASS_LOCAL;
        } else if (arg == "--zns-record") {
            global_zns_record_enabled = true;
        } else if (arg == "--recover-zns") {
            global_zns_record_enabled = true;
            recover_zns = true;
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
    // cout << "test 0: " << zns.zones_info[0].blk_info_vec[0].src << endl;
    report_device_info(zns.info);
    // report_zones_info(zns.zones, zns.nr_zones);
    if (!recover_zns && zbd_reset_zones(zns.dev_fd, 0, 0) != 0) {
        perror("zbd_reset_zones");
        return EXIT_FAILURE;
    }

    if (recover_zns) {
        recover_graph_from_zns_records();
        if (print_pma) {
            print_all_vertex_pma_blocks();
        }
        cout << "total_dump_time: " << total_dump_time << endl;
        return EXIT_SUCCESS;
    }

    if (load_graph(dataPath) != 0) {
        return EXIT_FAILURE;
    }
    cout << "Total number of edges: " << global_edge_num << " global_max_vid: " << global_max_vid << endl;
    cout << "global_num_queues: " << global_num_queues << endl;

    // std::random_device rd;  // 获取随机种子
    // std::mt19937 g(rd()); 
    // shuffle(global_edge_vec.begin(), global_edge_vec.end(), g);

    // graph_init(global_max_vid);
    // graph_init_new(global_max_vid);

    if (mode == "baseline") {
        graph_maintenance_baseline();
    } else if (mode == "new") {
        graph_maintenance_new();
    } else if (mode == "parallel") {
        graph_maintenance_parallel();
    } else if (mode == "parallel2") {
        graph_maintenance_parallel2();
    } else {
        cerr << "Invalid mode: " << mode << endl;
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    for (const auto &edge : insert_edges) {
        bool forward_inserted = graph_insert_directed_edge(edge.first, edge.second);
        bool backward_inserted = graph_insert_directed_edge(edge.second, edge.first);
        cout << "Insert edge (" << edge.first << ", " << edge.second << "): "
             << ((forward_inserted || backward_inserted) ? "inserted" : "already exists") << endl;
    }
    for (const auto &edge : delete_edges) {
        bool deleted = graph_delete_undirected_edge(edge.first, edge.second);
        cout << "Delete edge (" << edge.first << ", " << edge.second << "): "
             << (deleted ? "deleted" : "not found") << endl;
    }
    if ((!insert_edges.empty() || !delete_edges.empty()) &&
        global_cache_blk >= cache_blk_threshold) {
        if (use_zone_writer) {
            dumpToZoneZoneSerial();
        } else if (mode == "parallel2") {
            dumpToZoneParallel();
        } else {
            dumpToZone();
        }
        global_cache_blk = 0;
    }
    persist_cache_vertices_as_records();

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

    cout << "total_dump_time: " << total_dump_time << endl;
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
