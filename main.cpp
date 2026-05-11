#include "utils.h"
#include <string.h>
#include "pma.h"
#include "benchmark.h"
#include <iostream>
// #include <chrono> // test git
// #include "global.h"
using namespace std;

// 数据集路径表。
// 运行方式：./main <dataset_id> <thread_num>
// dataset_id 对应下面的下标。
vector<string> data_path;
void set_path(vector<string>& data_path)  
{
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

int main(int argc, char* argv[]) {
    // cout << "int32_t size: " << sizeof(int32_t) << endl;
    // cout << "int size: " << sizeof(int) << endl;
    
    set_path(data_path);
    string dataPath = data_path[atoi(argv[1])];
    THREAD_NUM = atoi(argv[2]);
    // 初始化全局 ZNS 设备状态

    memset(&zns, 0, sizeof(zns));
    zns.path = NVME_DEVICE_PATH;
    zns.dev_fd = -1;

    // 检测目标设备是否为 zoned block device
    if (zbd_device_is_zoned(zns.path)) {
        cout << "This is a ZBD" << endl;
    }

    zns_open();
    // cout << "test 0: " << zns.zones_info[0].blk_info_vec[0].src << endl;
    report_device_info(zns.info);
    // report_zones_info(zns.zones, zns.nr_zones);
    zbd_reset_zones(zns.dev_fd, 0, 0);


    load_graph(dataPath);
    cout << "Total number of edges: " << global_edge_num << " global_max_vid: " << global_max_vid << endl;
    cout << "global_num_queues: " << global_num_queues << endl;

    // std::random_device rd;  // 获取随机种子
    // std::mt19937 g(rd()); 
    // shuffle(global_edge_vec.begin(), global_edge_vec.end(), g);

    // graph_init(global_max_vid);
    // graph_init_new(global_max_vid);

    graph_maintenance_baseline();
    // graph_maintenance_parallel();
    graph_maintenance_parallel2();
    graph_recovery();
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

    return 1;
}
