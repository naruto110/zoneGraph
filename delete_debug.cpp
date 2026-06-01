#include "utils.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "pma.h"

using namespace std;

static void set_debug_path(vector<string> &data_path)
{
    data_path.clear();
    data_path.push_back("./data.txt");
    data_path.push_back("./validation_edges.txt");
    data_path.push_back("/home/yingsui/data/livejournal/out.e");
    data_path.push_back("/home/yingsui/data/orkut/out.e");
    data_path.push_back("/home/yingsui/data/skitter/out.e");
    data_path.push_back("/home/yingsui/data/youtube/out.e");
    data_path.push_back("/home/yingsui/data/dblp/out.e");
    data_path.push_back("/home/yingsui/data/roadCA/out.e");
    data_path.push_back("/home/yingsui/data/synthetic/syn19.e");
    data_path.push_back("/home/yingsui/data/synthetic/syn20.e");
    data_path.push_back("/home/yingsui/data/synthetic/syn21.e");
    data_path.push_back("/home/yingsui/data/synthetic/syn22.e");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/ca-GrQc.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/ca-CondMat.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/email-Enron.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/com-Amazon.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/com-DBLP.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/com-Youtube.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/twitch-gamers.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/as-Skitter.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/soc-Pokec.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/soc-LiveJournal1.txt");
    data_path.push_back("/home/femu/project/datasets/snap/edgelist/com-Orkut.txt");
}

static bool edge_exists_directed(int32_t src, int32_t dst)
{
    if (src < 0 || src >= static_cast<int32_t>(global_vectex_vec.size())) {
        return false;
    }
    Vertex *v = global_vectex_vec[src];
    shared_lock<std::shared_timed_mutex> lock(v->mutex);
    if (v->blk_list.empty()) {
        return cache_contains(v, dst);
    }
    return pending_insert_contains(v, dst) || pma_contains_scan(v, dst);
}

static bool edge_exists_directed_exact(int32_t src, int32_t dst)
{
    if (src < 0 || src >= static_cast<int32_t>(global_vectex_vec.size())) {
        return false;
    }
    Vertex *v = global_vectex_vec[src];
    shared_lock<std::shared_timed_mutex> lock(v->mutex);
    if (v->blk_list.empty()) {
        return cache_contains(v, dst);
    }
    if (pending_insert_contains(v, dst)) {
        return true;
    }
    int32_t tmp_blk_buf[PM_BLK_SIZE];
    for (int32_t blk_idx = 0; blk_idx < static_cast<int32_t>(v->blk_list.size()); blk_idx++) {
        pm_blk &blk = v->blk_list[blk_idx];
        if (blk.valid_count == 0) {
            continue;
        }
        int32_t *blk_ptr = blk.buf;
        if (blk_ptr == nullptr) {
            readBlockFromZoneInto(zns.dev_fd, tmp_blk_buf, PM_BLK_SIZE, blk.global_wp);
            blk_ptr = tmp_blk_buf;
        }
        int32_t local_end = pma_block_end(v->pma, blk_idx) - pma_block_start(blk_idx);
        for (int32_t bit_idx = 0; bit_idx < local_end; bit_idx++) {
            if (blk_ptr[bit_idx] == dst) {
                return true;
            }
        }
    }
    return false;
}

int main(int argc, char **argv)
{
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <dataset_id> <thread_num> <delete_file> <limit> [parallel]" << endl;
        return 2;
    }
    vector<string> data_path;
    set_debug_path(data_path);
    int dataset_id = atoi(argv[1]);
    THREAD_NUM = atoi(argv[2]);
    string delete_file = argv[3];
    int limit = atoi(argv[4]);
    bool parallel = argc >= 6 && string(argv[5]) == "parallel";
    if (dataset_id < 0 || dataset_id >= static_cast<int>(data_path.size()) || THREAD_NUM <= 0 || limit <= 0) {
        return 2;
    }

    memset(&zns, 0, sizeof(zns));
    zns.path = NVME_DEVICE_PATH;
    zns.dev_fd = -1;
    zns.direct_fd = -1;
    if (zns_open() != 0) {
        cerr << "Failed to open ZNS device" << endl;
        return 1;
    }
    if (zbd_reset_zones(zns.dev_fd, 0, 0) != 0) {
        perror("zbd_reset_zones");
        return 1;
    }
    reset_zns_write_pointers();
    resetRedoLogBuffer();

    if (load_graph(data_path[dataset_id]) != 0) {
        return 1;
    }
    graph_maintenance_parallel2();

    ifstream input(delete_file);
    if (!input) {
        cerr << "failed to open " << delete_file << endl;
        return 1;
    }
    vector<pair<int32_t, int32_t>> edges;
    int32_t src = 0;
    int32_t dst = 0;
    while (static_cast<int>(edges.size()) < limit && input >> src >> dst) {
        edges.push_back({src, dst});
    }

    int64_t forward_exists = 0;
    int64_t backward_exists = 0;
    int64_t forward_exists_exact = 0;
    int64_t backward_exists_exact = 0;
    for (const auto &edge : edges) {
        if (edge_exists_directed(edge.first, edge.second)) {
            forward_exists++;
        }
        if (edge_exists_directed(edge.second, edge.first)) {
            backward_exists++;
        }
        if (edge_exists_directed_exact(edge.first, edge.second)) {
            forward_exists_exact++;
        }
        if (edge_exists_directed_exact(edge.second, edge.first)) {
            backward_exists_exact++;
        }
    }

    int64_t forward_deleted = 0;
    int64_t backward_deleted = 0;
    int64_t forward_false_removed = 0;
    int64_t backward_false_removed = 0;
    double start = get_current_time();
    if (parallel) {
        shared_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
#pragma omp parallel for reduction(+ : forward_deleted, backward_deleted) schedule(dynamic, 1024)
        for (int64_t i = 0; i < static_cast<int64_t>(edges.size()); i++) {
            const auto &edge = edges[i];
            if (graph_delete_directed_edge_locked(edge.first, edge.second)) {
                forward_deleted++;
            }
            if (graph_delete_directed_edge_locked(edge.second, edge.first)) {
                backward_deleted++;
            }
        }
    } else {
        for (const auto &edge : edges) {
            bool fd = graph_delete_directed_edge(edge.first, edge.second);
            if (fd) {
                forward_deleted++;
            } else if (!edge_exists_directed_exact(edge.first, edge.second)) {
                forward_false_removed++;
            }
            bool bd = graph_delete_directed_edge(edge.second, edge.first);
            if (bd) {
                backward_deleted++;
            } else if (!edge_exists_directed_exact(edge.second, edge.first)) {
                backward_false_removed++;
            }
        }
    }
    double elapsed = get_current_time() - start;

    int64_t forward_after = 0;
    int64_t backward_after = 0;
    int64_t forward_after_exact = 0;
    int64_t backward_after_exact = 0;
    for (const auto &edge : edges) {
        if (edge_exists_directed(edge.first, edge.second)) {
            forward_after++;
        }
        if (edge_exists_directed(edge.second, edge.first)) {
            backward_after++;
        }
        if (edge_exists_directed_exact(edge.first, edge.second)) {
            forward_after_exact++;
        }
        if (edge_exists_directed_exact(edge.second, edge.first)) {
            backward_after_exact++;
        }
    }

    cout << "DEBUG_INPUT_UNDIRECTED=" << edges.size() << endl;
    cout << "DEBUG_PARALLEL=" << (parallel ? 1 : 0) << endl;
    cout << "DEBUG_FORWARD_EXISTS_BEFORE=" << forward_exists << endl;
    cout << "DEBUG_BACKWARD_EXISTS_BEFORE=" << backward_exists << endl;
    cout << "DEBUG_FORWARD_EXISTS_EXACT_BEFORE=" << forward_exists_exact << endl;
    cout << "DEBUG_BACKWARD_EXISTS_EXACT_BEFORE=" << backward_exists_exact << endl;
    cout << "DEBUG_FORWARD_DELETED=" << forward_deleted << endl;
    cout << "DEBUG_BACKWARD_DELETED=" << backward_deleted << endl;
    cout << "DEBUG_FORWARD_FALSE_REMOVED=" << forward_false_removed << endl;
    cout << "DEBUG_BACKWARD_FALSE_REMOVED=" << backward_false_removed << endl;
    cout << "DEBUG_FORWARD_EXISTS_AFTER=" << forward_after << endl;
    cout << "DEBUG_BACKWARD_EXISTS_AFTER=" << backward_after << endl;
    cout << "DEBUG_FORWARD_EXISTS_EXACT_AFTER=" << forward_after_exact << endl;
    cout << "DEBUG_BACKWARD_EXISTS_EXACT_AFTER=" << backward_after_exact << endl;
    cout << "DEBUG_DELETE_SEC=" << elapsed << endl;
    return 0;
}
