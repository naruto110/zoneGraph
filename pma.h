#ifndef __PMA_H_
#define __PMA_H_

#include "global.h"
// #include "utils.h"

using namespace std;

// 本文件是工程核心：
// 1. 实现 PMA 的查找、插入、重平衡与扩容；
// 2. 实现图加载、构图、刷盘和恢复主流程；
// 3. 负责把逻辑上的邻接表映射到 ZNS 设备块上。

/* Height-based (as opposed to depth-based) thresholds. */
/* Upper density thresholds. */
static const double t_h = 0.75;  /* root. */
static const double t_0 = 1.00;  /* leaves. */
/* Lower density thresholds. */
static const double p_h = 0.50;  /* root. */
static const double p_0 = 0.25;  /* leaves. */
static const uint8_t max_sparseness =  2; // 1 / p_0;
static const int32_t BLOCK_ALIGNED_SEG_LEN = PM_BLK_SIZE / 2;

#define BLK_SCALE 1
// #define PM_BLK_SIZE (32 * BLK_SCALE)  // 每块的元素数可在实验中调整

PMA* pma_creat (void);
static void refresh_loaded_block_meta(Vertex *v);
static void refresh_block_meta_range(Vertex *v, int32_t start_blk, int32_t end_blk);
static int32_t pma_block_start(int32_t blk_idx);
static int32_t pma_block_end(PMA *pma, int32_t blk_idx);
static bool pma_slot_in_range(PMA *pma, int32_t blk_idx, int32_t bit_idx);
static int32_t pma_block_base_height(PMA *pma);
static void ensure_pma_blocks_for_capacity(Vertex *v, PMA *pma);
static void rebuild_pma_from_sorted_values(Vertex *v, vector<int32_t> &values);
static bool pma_contains_scan(Vertex *v, int32_t des);
static int32_t pma_find_index_scan(Vertex *v, int32_t des);
static bool find_target_block_for_insert(Vertex *v, int32_t des, int32_t &target_blk);
static bool pending_insert_contains(Vertex *v, int32_t des);
static bool pending_insert_remove(Vertex *v, int32_t des);
static bool pending_insert_empty(Vertex *v);
static void pending_insert_clear(Vertex *v);
static int64_t pending_insert_count(Vertex *v);
static bool merge_pending_values_into_block(Vertex *v, int32_t blk_idx, vector<int32_t> &values);
static bool find_rebalance_window(Vertex *v, int32_t blk_idx,
                                  unordered_map<int32_t, vector<int32_t>> &pending_by_block,
                                  int32_t &out_window_start, int32_t &out_window_end);
static void resize_with_pending_values(Vertex *v, vector<int32_t> &values);
static void resize_with_pending_blocks(Vertex *v,
                                       unordered_map<int32_t, vector<int32_t>> &pending_by_block);
static void merge_pending_insert_buffer(Vertex *v);
static void merge_all_pending_insert_buffers(void);
static void maybe_merge_pending_insert_buffers(void);
static int32_t find_prev_valid_index(Vertex *src, int32_t from, int32_t index,
                                     unordered_map<int32_t, int32_t *> &local_blk_cache);
static void compute_capacity (PMA* pma);
int32_t binary_pma_find(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to);
int32_t pma_insert(PMA* pma, Vertex *src, int32_t des, int32_t index);
void rebalance(PMA* pma, Vertex *src, int32_t index);
void spread (PMA* pma, Vertex *src, int32_t from, int32_t to, int32_t occupancy);
void resize (PMA* pma, Vertex *src);
static inline int32_t pma_value_at(Vertex *src, int32_t index);
static int32_t count_window_occupancy(Vertex *src, int32_t from, int32_t to);

int load_graph(string path);
int scan_graph_metadata(string path);
void graph_init(int32_t max_vid);
void graph_maintenance_baseline(void);
bool cache_insert_unique(Vertex *v, int32_t des);
bool cache_delete(Vertex *v, int32_t des);
void pma_compact(PMA *pma, Vertex *src);
bool pma_delete(PMA *pma, Vertex *src, int32_t des);
bool graph_insert_directed_edge(int32_t src, int32_t des);
static bool graph_insert_directed_edge_locked(int32_t src, int32_t des);
bool graph_insert_undirected_edge(int32_t src, int32_t des);
bool graph_delete_directed_edge(int32_t src, int32_t des);
static bool graph_delete_directed_edge_locked(int32_t src, int32_t des);
bool graph_delete_undirected_edge(int32_t src, int32_t des);
VertexClass classify_vertex(Vertex *v);
const char *vertex_class_name(VertexClass vclass);
static void note_vertex_update(Vertex *v);
static ATZVertexType classify_atz_vertex(Vertex *v);
static const char *atz_vertex_type_name(ATZVertexType type);
static const char *atz_zone_type_name(ATZZoneType type);
static void update_vertex_temperatures(void);
static void print_atz_gc_summary(void);
bool requestClassifiedBlock(Vertex *v, unsigned long long &wp);
bool garbage_collection_greedy(void);
static void maybe_run_proactive_gc(void);
void print_vertex(int32_t vid, int64_t &total_edges);
void print_graph(void);
void print_vertex_pma_blocks(int32_t vid);
void print_all_vertex_pma_blocks(void);
static void print_pending_merge_debug_snapshot(const char *title);
static void print_pending_merge_block_values(const char *tag, Vertex *v, int32_t blk_idx);

// void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp);
void dumpToZone(void);
void dumpToZoneParallel(void);
void dumpToZoneZoneSerial(void);
void persist_cache_vertices_as_records(void);
void recover_graph_from_zns_records(void);
void run_recovery_experiment_from_zns(void);
void print_dram_metadata_usage(void);

void graph_init_new(int32_t max_vid);
void graph_maintenance_new(void);
void graph_maintenance_parallel(void);
void graph_maintenance_parallel2(void);
void graph_maintenance_parallel2_streaming(const string &path);

void graph_recovery(void);
int32_t binary_pma_find_recovery(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to);


// 创建 PMA 元数据，初始容量很小，后续依靠 resize 扩容。
PMA* pma_creat (void)
{
    PMA* pma = (PMA*)malloc(sizeof(PMA));
    pma->ele_num = 0;
    pma->seg_len = BLOCK_ALIGNED_SEG_LEN;
    pma->arr_cap = pma->seg_len * 2;
    pma->seg_num = pma->arr_cap / pma->seg_len;
    pma->tree_h = floor_lg (pma->seg_num) + 1;
    pma->delta_t = (t_0 - t_h) / pma->tree_h;  // 上界密度按层级线性变化
    pma->delta_p = (p_h - p_0) / pma->tree_h;
    // cout << "tree_h: " << pma->tree_h << " seg_num: " << pma->seg_num << endl;
    // pma->array = (element_type *)memkind_malloc(kind, sizeof(element_type) * pma->arr_cap);  // on NVM
    // pma->row_offset = (vertexOff *)malloc(sizeof(vertexOff) * pma->vertex_cap); // on DRAM
    return (pma);
}

static void refresh_loaded_block_meta(Vertex *v)
{
    for (int32_t blk_idx = 0; blk_idx < static_cast<int32_t>(v->blk_list.size()); blk_idx++) {
        if (v->blk_list[blk_idx].buf != nullptr) {
            refresh_pm_blk_meta(v->blk_list[blk_idx], v->blk_list[blk_idx].buf);
        }
    }
}

static void refresh_block_meta_range(Vertex *v, int32_t start_blk, int32_t end_blk)
{
    if (v->blk_list.empty()) {
        return;
    }
    start_blk = max(start_blk, 0);
    end_blk = min(end_blk, static_cast<int32_t>(v->blk_list.size()) - 1);
    for (int32_t blk_idx = start_blk; blk_idx <= end_blk; blk_idx++) {
        if (v->blk_list[blk_idx].buf != nullptr) {
            refresh_pm_blk_meta(v->blk_list[blk_idx], v->blk_list[blk_idx].buf);
        }
    }
}

static int32_t pma_block_start(int32_t blk_idx)
{
    return blk_idx * PM_BLK_SIZE;
}

static int32_t pma_block_end(PMA *pma, int32_t blk_idx)
{
    return min((blk_idx + 1) * PM_BLK_SIZE, pma->arr_cap);
}

static bool pma_slot_in_range(PMA *pma, int32_t blk_idx, int32_t bit_idx)
{
    int32_t idx = blk_idx * PM_BLK_SIZE + bit_idx;
    return idx >= 0 && idx < pma->arr_cap;
}

static int32_t pma_block_base_height(PMA *pma)
{
    if (pma->seg_len <= 0 || PM_BLK_SIZE < pma->seg_len ||
        PM_BLK_SIZE % pma->seg_len != 0) {
        return 0;
    }
    int32_t ratio = PM_BLK_SIZE / pma->seg_len;
    if ((ratio & (ratio - 1)) != 0) {
        return 0;
    }
    return min(static_cast<int32_t>(floor_lg(ratio)), max(pma->tree_h - 1, 0));
}

static void ensure_pma_blocks_for_capacity(Vertex *v, PMA *pma)
{
    int32_t blk_num = ceil_div(pma->arr_cap, PM_BLK_SIZE);
    while (static_cast<int32_t>(v->blk_list.size()) < blk_num) {
        pm_blk new_blk{};
        new_blk.buf = new int32_t[PM_BLK_SIZE];
        fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
        v->blk_list.push_back(new_blk);
        global_cache_blk++;
    }
}

static int32_t pma_find_index_scan(Vertex *v, int32_t des)
{
    if (v->pma == nullptr || v->blk_list.empty()) {
        return -1;
    }
    int32_t tmp_blk_buf[PM_BLK_SIZE];
    for (int32_t blk_idx = 0; blk_idx < static_cast<int32_t>(v->blk_list.size()); blk_idx++) {
        pm_blk &blk = v->blk_list[blk_idx];
        if (blk.valid_count == 0) {
            continue;
        }
        if (des < blk.min_val || des > blk.max_val) {
            continue;
        }
        int32_t *blk_ptr = blk.buf;
        if (blk_ptr == nullptr) {
            readBlockFromZoneInto(zns.dev_fd, tmp_blk_buf, PM_BLK_SIZE, blk.global_wp);
            refresh_pm_blk_meta(blk, tmp_blk_buf);
            blk_ptr = tmp_blk_buf;
        }
        int32_t local_end = pma_block_end(v->pma, blk_idx) - pma_block_start(blk_idx);
        if (blk.sorted) {
            int32_t lo = 0;
            int32_t hi = local_end - 1;
            while (lo <= hi) {
                int32_t mid = lo + (hi - lo) / 2;
                int32_t mid_idx = mid;
                while (mid_idx >= lo && blk_ptr[mid_idx] < 0) {
                    mid_idx--;
                }
                if (mid_idx < lo) {
                    lo = mid + 1;
                    continue;
                }
                if (blk_ptr[mid_idx] == des) {
                    return blk_idx * PM_BLK_SIZE + mid_idx;
                }
                if (blk_ptr[mid_idx] < des) {
                    lo = mid + 1;
                } else {
                    hi = mid_idx - 1;
                }
            }
            continue;
        }
        for (int32_t bit_idx = 0; bit_idx < local_end; bit_idx++) {
            if (blk_ptr[bit_idx] == des) {
                return blk_idx * PM_BLK_SIZE + bit_idx;
            }
        }
    }
    return -1;
}

static bool pma_contains_scan(Vertex *v, int32_t des)
{
    return pma_find_index_scan(v, des) >= 0;
}

static bool find_target_block_for_insert(Vertex *v, int32_t des, int32_t &target_blk)
{
    if (v->pma == nullptr || v->blk_list.empty()) {
        return false;
    }

    int32_t prev_max = INT32_MIN;
    for (int32_t blk_idx = 0; blk_idx < static_cast<int32_t>(v->blk_list.size()); blk_idx++) {
        int32_t local_end = pma_block_end(v->pma, blk_idx) - pma_block_start(blk_idx);
        if (local_end <= 0) {
            continue;
        }

        int32_t next_min = INT32_MAX;
        for (int32_t next_blk = blk_idx + 1;
             next_blk < static_cast<int32_t>(v->blk_list.size());
             next_blk++) {
            if (v->blk_list[next_blk].valid_count > 0) {
                next_min = v->blk_list[next_blk].min_val;
                break;
            }
        }

        if (des >= prev_max && des <= next_min) {
            target_blk = blk_idx;
            return true;
        }

        if (v->blk_list[blk_idx].max_val != INT32_MIN) {
            prev_max = v->blk_list[blk_idx].max_val;
        }
    }
    return false;
}


static bool pending_insert_contains(Vertex *v, int32_t des)
{
    if (find(v->pending_insert_unplaced.begin(),
             v->pending_insert_unplaced.end(),
             des) != v->pending_insert_unplaced.end()) {
        return true;
    }
    for (auto &entry : v->pending_insert_by_block) {
        if (find(entry.second.begin(), entry.second.end(), des) != entry.second.end()) {
            return true;
        }
    }
    return false;
}

static bool pending_insert_remove(Vertex *v, int32_t des)
{
    auto unplaced_it = find(v->pending_insert_unplaced.begin(),
                            v->pending_insert_unplaced.end(),
                            des);
    if (unplaced_it != v->pending_insert_unplaced.end()) {
        v->pending_insert_unplaced.erase(unplaced_it);
        global_pending_edge_count--;
        return true;
    }
    for (auto it = v->pending_insert_by_block.begin();
         it != v->pending_insert_by_block.end();
         ++it) {
        auto value_it = find(it->second.begin(), it->second.end(), des);
        if (value_it != it->second.end()) {
                it->second.erase(value_it);
                if (it->second.empty()) {
                    v->pending_insert_by_block.erase(it);
                }
                global_pending_edge_count--;
                return true;
            }
    }
    return false;
}

static bool pending_insert_empty(Vertex *v)
{
    return v->pending_insert_unplaced.empty() &&
           v->pending_insert_by_block.empty();
}

static void pending_insert_clear(Vertex *v)
{
    v->pending_insert_unplaced.clear();
    v->pending_insert_by_block.clear();
}

static int64_t pending_insert_count(Vertex *v)
{
    int64_t count = static_cast<int64_t>(v->pending_insert_unplaced.size());
    for (const auto &entry : v->pending_insert_by_block) {
        count += static_cast<int64_t>(entry.second.size());
    }
    return count;
}


static void print_pending_merge_block_values(const char *tag, Vertex *v, int32_t blk_idx)
{
    if (!global_debug_pending_merge || v == nullptr ||
        blk_idx < 0 || blk_idx >= static_cast<int32_t>(v->blk_list.size())) {
        return;
    }

    int32_t tmp_blk_buf[PM_BLK_SIZE];
    int32_t *blk_ptr = v->blk_list[blk_idx].buf;
    if (blk_ptr == nullptr) {
        readBlockFromZoneInto(zns.dev_fd, tmp_blk_buf, PM_BLK_SIZE,
                              v->blk_list[blk_idx].global_wp);
        blk_ptr = tmp_blk_buf;
    }

    cout << tag
         << " vid=" << v->id
         << " blk=" << blk_idx
         << " valid=" << v->blk_list[blk_idx].valid_count
         << " values=[";
    for (int32_t i = 0; i < PM_BLK_SIZE; i++) {
        cout << blk_ptr[i];
        if (i + 1 < PM_BLK_SIZE) {
            cout << ",";
        }
    }
    cout << "]" << endl;
}

static bool merge_pending_values_into_block(Vertex *v, int32_t blk_idx, vector<int32_t> &values)
{
    if (values.empty()) {
        return true;
    }
    if (blk_idx < 0 || blk_idx >= static_cast<int32_t>(v->blk_list.size())) {
        return false;
    }

    int32_t local_end = pma_block_end(v->pma, blk_idx) - pma_block_start(blk_idx);
    if (local_end <= 0) {
        return false;
    }

    pm_blk &blk = v->blk_list[blk_idx];
    int32_t free_slots = local_end - blk.valid_count;
    if (static_cast<int32_t>(values.size()) > free_slots) {
        return false;
    }

    if (blk.buf == nullptr) {
        cacheBlock(v, blk_idx);
    }

    size_t value_idx = 0;
    for (int32_t bit_idx = 0;
         bit_idx < local_end && value_idx < values.size();
         bit_idx++) {
        if (blk.buf[bit_idx] < 0) {
            blk.buf[bit_idx] = values[value_idx++];
        }
    }
    if (value_idx != values.size()) {
        return false;
    }

    v->pma->ele_num += static_cast<int32_t>(values.size());
    refresh_pm_blk_meta(blk, blk.buf);
    if (global_debug_pending_merge) {
        cout << "PENDING_MERGE_FULL_PMA vid=" << v->id << endl;
        print_vertex_pma_blocks(v->id);
    }
    return true;
}

static void collect_block_values(Vertex *v, int32_t blk_idx,
                                 int32_t from, int32_t to,
                                 vector<int32_t> &values)
{
    if (v->blk_list[blk_idx].buf == nullptr) {
        cacheBlock(v, blk_idx);
    }
    int32_t local_from = max(from - pma_block_start(blk_idx), 0);
    int32_t local_to = min(to - pma_block_start(blk_idx), PM_BLK_SIZE) - 1;
    for (int32_t bit_idx = local_from; bit_idx <= local_to; bit_idx++) {
        if (v->blk_list[blk_idx].buf[bit_idx] >= 0) {
            values.push_back(v->blk_list[blk_idx].buf[bit_idx]);
        }
    }
}

static void spread_sorted_values_to_window(Vertex *v, int32_t from, int32_t to,
                                           const vector<int32_t> &values)
{
    int32_t capacity = to - from;
    int32_t start_blk = from / PM_BLK_SIZE;
    int32_t end_blk = (to - 1) / PM_BLK_SIZE;
    for (int32_t blk_idx = start_blk; blk_idx <= end_blk; blk_idx++) {
        if (v->blk_list[blk_idx].buf == nullptr) {
            cacheBlock(v, blk_idx);
        }
        int32_t local_from = blk_idx == start_blk ? from % PM_BLK_SIZE : 0;
        int32_t local_to = blk_idx == end_blk ? (to - 1) % PM_BLK_SIZE : PM_BLK_SIZE - 1;
        for (int32_t bit_idx = local_from; bit_idx <= local_to; bit_idx++) {
            v->blk_list[blk_idx].buf[bit_idx] = -1;
        }
    }

    if (!values.empty()) {
        int32_t frequency = (capacity << 8) / static_cast<int32_t>(values.size());
        int32_t write_index = (to << 8) - frequency;
        for (int32_t ridx = static_cast<int32_t>(values.size()) - 1; ridx >= 0; ridx--) {
            int32_t target = write_index >> 8;
            target = min(max(target, from), to - 1);
            int32_t blk_idx = target / PM_BLK_SIZE;
            int32_t bit_idx = target % PM_BLK_SIZE;
            v->blk_list[blk_idx].buf[bit_idx] = values[ridx];
            write_index -= frequency;
        }
    }
    refresh_block_meta_range(v, start_blk, end_blk);
}

static bool find_rebalance_window(Vertex *v, int32_t blk_idx,
                                  unordered_map<int32_t, vector<int32_t>> &pending_by_block,
                                  int32_t &out_window_start, int32_t &out_window_end)
{
    if (blk_idx < 0 || blk_idx >= static_cast<int32_t>(v->blk_list.size())) {
        return false;
    }

    PMA *pma = v->pma;
    int32_t index = pma_block_start(blk_idx);
    int32_t height = pma_block_base_height(pma);
    double density = 0.0;
    double t_height = 0.0;
    double p_height = 0.0;

    do {
        int32_t window_size = pma->seg_len * (1 << height);
        int32_t window = index / window_size;
        out_window_start = window * window_size;
        out_window_end = min(out_window_start + window_size, pma->arr_cap);
        int32_t start_blk = out_window_start / PM_BLK_SIZE;
        int32_t end_blk = (out_window_end - 1) / PM_BLK_SIZE;

        int32_t inserted_count = 0;
        for (int32_t pb = start_blk; pb <= end_blk; pb++) {
            auto it = pending_by_block.find(pb);
            if (it != pending_by_block.end())
                inserted_count += static_cast<int32_t>(it->second.size());
        }
        int32_t occupancy = count_window_occupancy(v, out_window_start, out_window_end)
                            + inserted_count;
        density = static_cast<double>(occupancy) /
                  static_cast<double>(out_window_end - out_window_start);
        t_height = t_0 - (height * pma->delta_t);
        p_height = p_0 + (height * pma->delta_p);
        height++;
    } while ((density < p_height || density >= t_height) && height < pma->tree_h);

    return !(density < p_height || density >= t_height);
}

static void resize_with_pending_values(Vertex *v, vector<int32_t> &values)
{
    vector<int32_t> all_values;
    all_values.reserve(v->pma->ele_num + values.size());
    int32_t tmp_blk_buf[PM_BLK_SIZE];
    for (int32_t blk_idx = 0; blk_idx < static_cast<int32_t>(v->blk_list.size()); blk_idx++) {
        int32_t *blk_ptr = v->blk_list[blk_idx].buf;
        if (blk_ptr == nullptr) {
            readBlockFromZoneInto(zns.dev_fd, tmp_blk_buf, PM_BLK_SIZE,
                                  v->blk_list[blk_idx].global_wp);
            blk_ptr = tmp_blk_buf;
        }
        int32_t local_end = pma_block_end(v->pma, blk_idx) - pma_block_start(blk_idx);
        for (int32_t bit_idx = 0; bit_idx < local_end; bit_idx++) {
            if (blk_ptr[bit_idx] >= 0) {
                all_values.push_back(blk_ptr[bit_idx]);
            }
        }
    }
    all_values.insert(all_values.end(), values.begin(), values.end());
    sort(all_values.begin(), all_values.end());
    rebuild_pma_from_sorted_values(v, all_values);
    if (global_debug_pending_merge) {
        cout << "PENDING_RESIZE_FULL_PMA vid=" << v->id << endl;
        print_vertex_pma_blocks(v->id);
    }
}

static void resize_with_pending_blocks(Vertex *v,
                                       unordered_map<int32_t, vector<int32_t>> &pending_by_block)
{
    vector<int32_t> all_values;
    all_values.reserve(v->pma->ele_num);
    for (const auto &entry : pending_by_block) {
        all_values.reserve(all_values.capacity() + entry.second.size());
    }

    for (int32_t blk_idx = 0; blk_idx < static_cast<int32_t>(v->blk_list.size()); blk_idx++) {
        vector<int32_t> block_values;
        int32_t block_from = pma_block_start(blk_idx);
        int32_t block_to = pma_block_end(v->pma, blk_idx);
        if (block_to <= block_from) {
            continue;
        }
        collect_block_values(v, blk_idx, block_from, block_to, block_values);
        auto pending_it = pending_by_block.find(blk_idx);
        if (pending_it != pending_by_block.end()) {
            block_values.insert(block_values.end(),
                                pending_it->second.begin(),
                                pending_it->second.end());
        }
        sort(block_values.begin(), block_values.end());
        block_values.erase(unique(block_values.begin(), block_values.end()),
                           block_values.end());
        all_values.insert(all_values.end(), block_values.begin(), block_values.end());
    }

    rebuild_pma_from_sorted_values(v, all_values);
    pending_by_block.clear();
    if (global_debug_pending_merge) {
        cout << "PENDING_RESIZE_FULL_PMA vid=" << v->id << endl;
        print_vertex_pma_blocks(v->id);
    }
}

static void rebuild_pma_from_sorted_values(Vertex *v, vector<int32_t> &values)
{
    values.erase(unique(values.begin(), values.end()), values.end());
    int32_t old_arr_cap = v->pma->arr_cap;
    v->pma->ele_num = static_cast<int32_t>(values.size());
    if (v->pma->arr_cap <= v->pma->ele_num ||
        static_cast<double>(v->pma->ele_num) / v->pma->arr_cap >= t_h ||
        v->pma->seg_len != BLOCK_ALIGNED_SEG_LEN) {
        compute_capacity(v->pma);
        if (v->pma->arr_cap < old_arr_cap) {
            v->pma->seg_num = hyperceil(ceil_div(old_arr_cap, v->pma->seg_len));
            v->pma->arr_cap = v->pma->seg_len * v->pma->seg_num;
        }
    }
    v->pma->tree_h = floor_lg(v->pma->seg_num) + 1;
    v->pma->delta_t = (t_0 - t_h) / v->pma->tree_h;
    v->pma->delta_p = (p_h - p_0) / v->pma->tree_h;
    ensure_pma_blocks_for_capacity(v, v->pma);

    for (int32_t blk_idx = 0; blk_idx < static_cast<int32_t>(v->blk_list.size()); blk_idx++) {
        if (v->blk_list[blk_idx].buf == nullptr) {
            cacheBlock(v, blk_idx);
        }
        fill(v->blk_list[blk_idx].buf, v->blk_list[blk_idx].buf + PM_BLK_SIZE, -1);
    }

    if (!values.empty()) {
        int32_t capacity = v->pma->arr_cap;
        int32_t frequency = (capacity << 8) / static_cast<int32_t>(values.size());
        int32_t write_index = (capacity << 8) - frequency;
        for (int32_t ridx = static_cast<int32_t>(values.size()) - 1; ridx >= 0; ridx--) {
            int32_t target = write_index >> 8;
            if (target >= v->pma->arr_cap) {
                target = v->pma->arr_cap - 1;
            }
            int32_t blk_idx = target / PM_BLK_SIZE;
            int32_t bit_idx = target % PM_BLK_SIZE;
            assert(pma_slot_in_range(v->pma, blk_idx, bit_idx));
            v->blk_list[blk_idx].buf[bit_idx] = values[ridx];
            write_index -= frequency;
        }
    }
    v->degree = static_cast<int32_t>(values.size());
    refresh_loaded_block_meta(v);
}

static void merge_pending_insert_buffer(Vertex *v)
{
    if (pending_insert_empty(v)) return;

    int64_t consumed_pending = pending_insert_count(v);
    unordered_map<int32_t, vector<int32_t>> pending_by_block = v->pending_insert_by_block;
    vector<int32_t> unplaced = v->pending_insert_unplaced;
    pending_insert_clear(v);

    // Route unplaced values into pending_by_block; collect truly unroutable separately
    vector<int32_t> unroutable;
    for (int32_t des : unplaced) {
        int32_t target_blk = -1;
        if (find_target_block_for_insert(v, des, target_blk)) {
            pending_by_block[target_blk].push_back(des);
            pm_blk &blk = v->blk_list[target_blk];
            blk.min_val = min(blk.min_val, des);
            blk.max_val = max(blk.max_val, des);
        } else {
            unroutable.push_back(des);
        }
    }
    if (!unroutable.empty()) {
        // Merge all pending into unroutable and rebuild
        for (auto &[blk_id, vals] : pending_by_block)
            unroutable.insert(unroutable.end(), vals.begin(), vals.end());
        resize_with_pending_values(v, unroutable);
        global_pending_edge_count -= consumed_pending;
        return;
    }

    if (pending_by_block.empty()) {
        global_pending_edge_count -= consumed_pending;
        return;
    }

    // Phase 1: direct merge — insert pending values into free slots
    vector<int32_t> sorted_blk_ids;
    sorted_blk_ids.reserve(pending_by_block.size());
    for (auto &entry : pending_by_block) sorted_blk_ids.push_back(entry.first);
    sort(sorted_blk_ids.begin(), sorted_blk_ids.end());

    vector<int32_t> needs_rebalance;
    for (int32_t blk_idx : sorted_blk_ids) {
        auto it = pending_by_block.find(blk_idx);
        if (it == pending_by_block.end() || it->second.empty()) continue;
        if (merge_pending_values_into_block(v, blk_idx, it->second)) {
            pending_by_block.erase(blk_idx);
        } else {
            needs_rebalance.push_back(blk_idx);
        }
    }

    if (needs_rebalance.empty()) {
        global_pending_edge_count -= consumed_pending;
        return;
    }

    // Phase 2: find rebalance window for each block that needs it
    // All pending counts (including other failing blocks) are included in occupancy
    vector<pair<int32_t, int32_t>> slot_windows;
    slot_windows.reserve(needs_rebalance.size());
    for (int32_t blk_idx : needs_rebalance) {
        int32_t ws = 0, we = 0;
        if (!find_rebalance_window(v, blk_idx, pending_by_block, ws, we)) {
            resize_with_pending_blocks(v, pending_by_block);
            global_pending_edge_count -= consumed_pending;
            return;
        }
        slot_windows.push_back({ws, we});
    }

    // Phase 3: merge overlapping windows into unified intervals
    sort(slot_windows.begin(), slot_windows.end());
    vector<pair<int32_t, int32_t>> merged_windows;
    for (auto &w : slot_windows) {
        if (merged_windows.empty() || w.first >= merged_windows.back().second) {
            merged_windows.push_back(w);
        } else {
            merged_windows.back().second = max(merged_windows.back().second, w.second);
        }
    }

    // Phase 4: spread each merged window
    for (auto &[window_start, window_end] : merged_windows) {
        int32_t start_blk = window_start / PM_BLK_SIZE;
        int32_t end_blk = (window_end - 1) / PM_BLK_SIZE;

        int32_t inserted_count = 0;
        for (int32_t pb = start_blk; pb <= end_blk; pb++) {
            auto it = pending_by_block.find(pb);
            if (it != pending_by_block.end())
                inserted_count += static_cast<int32_t>(it->second.size());
        }

        // Re-check density for the merged window (may be wider than each individual window)
        int32_t occupancy = count_window_occupancy(v, window_start, window_end) + inserted_count;
        double density = static_cast<double>(occupancy) /
                         static_cast<double>(window_end - window_start);
        if (density >= t_h) {
            resize_with_pending_blocks(v, pending_by_block);
            global_pending_edge_count -= consumed_pending;
            return;
        }

        vector<int32_t> window_values;
        window_values.reserve(occupancy);
        for (int32_t cur_blk = start_blk; cur_blk <= end_blk; cur_blk++) {
            vector<int32_t> block_values;
            int32_t block_from = max(window_start, pma_block_start(cur_blk));
            int32_t block_to = min(window_end, pma_block_end(v->pma, cur_blk));
            collect_block_values(v, cur_blk, block_from, block_to, block_values);
            auto pit = pending_by_block.find(cur_blk);
            if (pit != pending_by_block.end()) {
                block_values.insert(block_values.end(),
                                    pit->second.begin(), pit->second.end());
            }
            sort(block_values.begin(), block_values.end());
            block_values.erase(unique(block_values.begin(), block_values.end()),
                               block_values.end());
            window_values.insert(window_values.end(),
                                 block_values.begin(), block_values.end());
        }

        spread_sorted_values_to_window(v, window_start, window_end, window_values);
        v->pma->ele_num += inserted_count;

        if (global_debug_pending_merge) {
            cout << "PENDING_REBALANCE vid=" << v->id
                 << " window=[" << window_start << "," << window_end << ")" << endl;
            print_vertex_pma_blocks(v->id);
        }

        for (int32_t cur_blk = start_blk; cur_blk <= end_blk; cur_blk++) {
            pending_by_block.erase(cur_blk);
        }
    }

    if (!pending_by_block.empty()) {
        resize_with_pending_blocks(v, pending_by_block);
    }

    global_pending_edge_count -= consumed_pending;
}

static void merge_all_pending_insert_buffers(void)
{
    lock_guard<std::mutex> merge_lock(global_pending_merge_mutex);
    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr) {
            continue;
        }
        unique_lock<std::shared_timed_mutex> lock(v->mutex);
        if (pending_insert_empty(v)) {
            continue;
        }
        merge_pending_insert_buffer(v);
    }
}

static bool pma_debug_snapshot(vector<pair<int32_t, int32_t>> &snapshot)
{
    snapshot.clear();
    snapshot.reserve(global_vectex_vec.size());
    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr || v->pma == nullptr) {
            snapshot.push_back({0, 0});
            continue;
        }
        snapshot.push_back({v->pma->ele_num, v->pma->arr_cap});
    }
    return true;
}

static void print_pma_debug_capacity_changes(const vector<pair<int32_t, int32_t>> &before)
{
    cout << "  PMA capacity changes:";
    bool changed = false;
    for (int32_t vidx = 0; vidx < static_cast<int32_t>(global_vectex_vec.size()); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        if (v == nullptr || v->pma == nullptr ||
            vidx >= static_cast<int32_t>(before.size())) {
            continue;
        }
        if (before[vidx].first != v->pma->ele_num ||
            before[vidx].second != v->pma->arr_cap) {
            cout << " v" << vidx
                 << "(ele " << before[vidx].first << "->" << v->pma->ele_num
                 << ", cap " << before[vidx].second << "->" << v->pma->arr_cap
                 << ")";
            changed = true;
        }
    }
    if (!changed) {
        cout << " none";
    }
    cout << endl;
}

static bool pending_merge_needs_rebalance_or_resize(void)
{
    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr || v->pma == nullptr || v->blk_list.empty()) {
            continue;
        }
        if (!v->pending_insert_unplaced.empty()) {
            return true;
        }
        for (const auto &entry : v->pending_insert_by_block) {
            int32_t blk_idx = entry.first;
            if (blk_idx < 0 || blk_idx >= static_cast<int32_t>(v->blk_list.size())) {
                return true;
            }
            int32_t local_end = pma_block_end(v->pma, blk_idx) - pma_block_start(blk_idx);
            int32_t free_slots = local_end - v->blk_list[blk_idx].valid_count;
            if (static_cast<int32_t>(entry.second.size()) > free_slots) {
                return true;
            }
        }
    }
    return false;
}

static void maybe_merge_pending_insert_buffers(void)
{
    if (pending_merge_threshold <= 0 ||
        global_pending_edge_count.load() < pending_merge_threshold) {
        return;
    }
    unique_lock<std::mutex> merge_lock(global_pending_merge_mutex, std::try_to_lock);
    if (!merge_lock.owns_lock() ||
        global_pending_edge_count.load() < pending_merge_threshold) {
        return;
    }
    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr) {
            continue;
        }
        unique_lock<std::shared_timed_mutex> lock(v->mutex);
        if (pending_insert_empty(v)) {
            continue;
        }
        merge_pending_insert_buffer(v);
    }
}

static int32_t find_prev_valid_index(Vertex *src, int32_t from, int32_t index,
                                     unordered_map<int32_t, int32_t *> &local_blk_cache)
{
    while (index >= from) {
        int32_t blk_idx = index / PM_BLK_SIZE;
        int32_t bit_idx = index % PM_BLK_SIZE;
        if (blk_idx < 0 || blk_idx >= static_cast<int32_t>(src->blk_list.size())) {
            return -1;
        }

        pm_blk &blk = src->blk_list[blk_idx];
        if (blk.valid_count == 0) {
            index = blk_idx * PM_BLK_SIZE - 1;
            continue;
        }

        const int32_t *data = blk.buf;
        if (data == nullptr) {
            auto it = local_blk_cache.find(blk_idx);
            if (it == local_blk_cache.end()) {
                int32_t *tmp = nullptr;
                readBlockFromZone(zns.dev_fd, tmp, PM_BLK_SIZE, blk.global_wp);
                local_blk_cache[blk_idx] = tmp;
                data = tmp;
            } else {
                data = it->second;
            }
        }

        int32_t scan_from = 0;
        if (from / PM_BLK_SIZE == blk_idx) {
            scan_from = from % PM_BLK_SIZE;
        }
        for (int32_t bit = bit_idx; bit >= scan_from; bit--) {
            if (data[bit] >= 0) {
                return blk_idx * PM_BLK_SIZE + bit;
            }
        }
        index = blk_idx * PM_BLK_SIZE - 1;
    }
    return -1;
}

// 根据当前元素数重新计算 PMA 容量和段长。
static void compute_capacity (PMA* pma) 
{
    /* Keep leaf segments power-of-two and aligned with PM_BLK_SIZE. */
    pma->seg_len = BLOCK_ALIGNED_SEG_LEN;
    int32_t target_cap = max(pma->seg_len * 2,
                             max_sparseness * max(pma->ele_num, 1));
    pma->seg_num = hyperceil(ceil_div(target_cap, pma->seg_len));
    pma->arr_cap = pma->seg_len * pma->seg_num;

    assert(pma->arr_cap > pma->ele_num);
}

static inline int32_t pma_value_at(Vertex *src, int32_t index) {
    if (index < 0) {
        return -1;
    }

    int32_t blk_idx = index / PM_BLK_SIZE;
    int32_t bit_idx = index % PM_BLK_SIZE;
    if (src->blk_list[blk_idx].buf == nullptr) {
        cacheBlock(src, blk_idx);
    }
    return src->blk_list[blk_idx].buf[bit_idx];
}

static inline int32_t pma_value_at_readonly(Vertex *src, int32_t index) {
    if (index < 0) {
        return -1;
    }

    int32_t blk_idx = index / PM_BLK_SIZE;
    int32_t bit_idx = index % PM_BLK_SIZE;
    if (src->blk_list[blk_idx].buf != nullptr) {
        return src->blk_list[blk_idx].buf[bit_idx];
    }

    int32_t *tmp_blk_ptr = nullptr;
    readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE,
                      src->blk_list[blk_idx].global_wp);
    int32_t value = tmp_blk_ptr[bit_idx];
    delete[] tmp_blk_ptr;
    return value;
}

static int32_t count_window_occupancy(Vertex *src, int32_t from, int32_t to) {
    int32_t occupancy = 0;
    int32_t start_blk = from / PM_BLK_SIZE;
    int32_t end_blk   = (to - 1) / PM_BLK_SIZE;

    for (int32_t b = start_blk; b <= end_blk; b++) {
        int32_t local_from = (b == start_blk) ? (from % PM_BLK_SIZE) : 0;
        int32_t local_to   = (b == end_blk)   ? ((to - 1) % PM_BLK_SIZE) : (PM_BLK_SIZE - 1);

        if (src->blk_list[b].buf != nullptr) {
            for (int32_t bit = local_from; bit <= local_to; bit++) {
                if (src->blk_list[b].buf[bit] >= 0) {
                    occupancy++;
                }
            }
        } else {
            int32_t *tmp = nullptr;
            readBlockFromZone(zns.dev_fd, tmp, PM_BLK_SIZE, src->blk_list[b].global_wp);
            for (int32_t bit = local_from; bit <= local_to; bit++) {
                if (tmp[bit] >= 0) {
                    occupancy++;
                }
            }
            delete[] tmp;
        }
    }
    return occupancy;
}


static inline int32_t get_blk_val(Vertex *src, int32_t global_idx,
                                   unordered_map<int32_t, int32_t *> &local_blk_cache)
{
    int32_t blk_idx = global_idx / PM_BLK_SIZE;
    int32_t bit_idx = global_idx % PM_BLK_SIZE;
    if (src->blk_list[blk_idx].buf != nullptr) {
        return src->blk_list[blk_idx].buf[bit_idx];
    }
    auto it = local_blk_cache.find(blk_idx);
    if (it != local_blk_cache.end()) {
        return it->second[bit_idx];
    }
    return -2; // should not happen after find_prev_valid_index loaded the block
}

static void release_local_blk_cache(unordered_map<int32_t, int32_t *> &cache)
{
    for (auto &kv : cache) {
        delete[] kv.second;
    }
    cache.clear();
}

// 在 PMA 中查找 des 的前驱位置。
// 若访问到未缓存的块，会按需从设备换入。
int32_t binary_pma_find(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to)
{
    if (to < 0) {
        return -1;
    }

    int32_t index;
    unordered_map<int32_t, int32_t *> local_blk_cache;

    while (from <= to) {
        int32_t mid = from + (to - from) / 2;
        int32_t i = find_prev_valid_index(src, from, mid, local_blk_cache);

        if (i < from) {  /* Everything between from and mid (inclusive) is empty. */
            from = mid + 1;
        } else {
            int32_t val = get_blk_val(src, i, local_blk_cache);
            if (val == des) {
                index = i;
                release_local_blk_cache(local_blk_cache);
                return index;
            }
            else if (val < des) {
                from = mid + 1;
            } else {  /* val > des */
                to = i - 1;  // 保证最终得到的是小于 des 的最大位置
            }
        }
    }
    /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
    index = to;
    if (index < 0) {
        release_local_blk_cache(local_blk_cache);
        return -1;
    }

    index = find_prev_valid_index(src, 0, index, local_blk_cache);

    release_local_blk_cache(local_blk_cache);
    return index;
}

// 恢复阶段使用的查找版本；若 des 已存在则返回 -1。
int32_t binary_pma_find_recovery(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to)
{
    if (to < 0) {
        return -1;
    }

    int32_t index;
    unordered_map<int32_t, int32_t *> local_blk_cache;

    while (from <= to) {
        int32_t mid = from + (to - from) / 2;
        int32_t i = find_prev_valid_index(src, from, mid, local_blk_cache);

        if (i < from) {  /* Everything between from and mid (inclusive) is empty. */
            from = mid + 1;
        } else {
            int32_t val = get_blk_val(src, i, local_blk_cache);
            if (val == des) {
                index = i;
                release_local_blk_cache(local_blk_cache);
                return -1; // des already exists
            }
            else if (val < des) {
                from = mid + 1;
            } else {  /* val > des */
                to = i - 1;  // 保证最终得到的是小于 des 的最大位置
            }
        }
    }
    /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
    index = to;
    if (index < 0) {
        release_local_blk_cache(local_blk_cache);
        return -1;
    }

    index = find_prev_valid_index(src, 0, index, local_blk_cache);

    release_local_blk_cache(local_blk_cache);
    return index;
}


// 将 des 插入到 index 之后，插入完成后会触发 rebalance。
int32_t pma_insert(PMA* pma, Vertex *src, int32_t des, int32_t index)
{
    /*
        在 index 之后插入 des。
        插入过程中会按需把相关块换入 DRAM。
    */
    if (index >= 0 && pma_value_at(src, index) == des) {
        return index;
    }

    if (pma->ele_num >= pma->arr_cap) {
        resize(pma, src);
        index = binary_pma_find(pma, src, des, 0, pma->arr_cap - 1);
        if (index >= 0 && pma_value_at(src, index) == des) {
            return index;
        }
    }

    int32_t j = index + 1;
    int32_t blk_idx, bit_idx;
    int32_t inserted_index = index + 1;
    bool inserted = false;
    // int32_t cur_blk = -1;
    // int32_t *tmp_blk_ptr = new int32_t[PM_BLK_SIZE];

    while (j < pma->arr_cap) { // 向右寻找空槽
        blk_idx = j / PM_BLK_SIZE;
        bit_idx = j % PM_BLK_SIZE;

        // cout << "src->blk_list size: " << src->blk_list.size() << " blk_idx: " << blk_idx << endl;

        if (src->blk_list[blk_idx].buf == nullptr) {
            cacheBlock(src, blk_idx); 
        }
     
        if (src->blk_list[blk_idx].buf[bit_idx] < 0)
            break;
        j++;
    }
    
    if (j < pma->arr_cap) {
        while (j > index + 1) /* Push elements to make space for the new element. */
        {
            blk_idx = j / PM_BLK_SIZE;  // 
            bit_idx = j % PM_BLK_SIZE;
            if (src->blk_list[blk_idx].buf == nullptr) {
                cacheBlock(src, blk_idx);  // 把目标块换入 DRAM
            }

            int32_t pre_blk_idx = (j-1) / PM_BLK_SIZE;  // 通过右移元素腾出空位
            int32_t pre_bit_idx = (j-1) % PM_BLK_SIZE;
            if (src->blk_list[pre_blk_idx].buf == nullptr) {
                cacheBlock(src, pre_blk_idx);  // 把源块换入 DRAM
            }

            src->blk_list[blk_idx].buf[bit_idx] = src->blk_list[pre_blk_idx].buf[pre_bit_idx];
            j--;
        }

        blk_idx = (index + 1) / PM_BLK_SIZE;  // 真正的插入位置
        bit_idx = (index + 1) % PM_BLK_SIZE;

        if (src->blk_list[blk_idx].buf == nullptr) {
            cacheBlock(src, blk_idx);  // 把目标块换入 DRAM
        }
        src->blk_list[blk_idx].buf[bit_idx] = des;
        inserted_index = index + 1;  // 更新为元素实际插入位置，供后续 rebalance 使用
        inserted = true;
    } else { /* No empty space to the right of i. Try left. */
        // int32_t be_cur_blk = -1;
        // cout << "test --1 des: " << des << " index: " << index << endl;
        j = index - 1;
        while (j >= 0) {
            blk_idx = j / PM_BLK_SIZE;
            bit_idx = j % PM_BLK_SIZE;

            if (src->blk_list[blk_idx].buf == nullptr) {
                cacheBlock(src, blk_idx);
            }

            if (src->blk_list[blk_idx].buf[bit_idx] < 0)
                break;
            j--;
        }
        if (j >= 0) {
            while (j < index) {
                blk_idx = j / PM_BLK_SIZE;
                bit_idx = j % PM_BLK_SIZE;

                if (src->blk_list[blk_idx].buf == nullptr) {
                    cacheBlock(src, blk_idx);  // 把目标块换入 DRAM
                }

                int32_t be_blk_idx = (j+1) / PM_BLK_SIZE;
                int32_t be_bit_idx = (j+1) % PM_BLK_SIZE;

                if (src->blk_list[be_blk_idx].buf == nullptr) {
                    cacheBlock(src, be_blk_idx);  // 把源块换入 DRAM
                }

                src->blk_list[blk_idx].buf[bit_idx] = src->blk_list[be_blk_idx].buf[be_bit_idx];
                j++;
            }
            blk_idx = (index) / PM_BLK_SIZE; 
            bit_idx = (index) % PM_BLK_SIZE;

            if (src->blk_list[blk_idx].buf == nullptr) {
                cacheBlock(src, blk_idx);  // 把目标块换入 DRAM
            }
            src->blk_list[blk_idx].buf[bit_idx] = des;
            inserted_index = index;
            inserted = true;
        }
    }
    if (!inserted) {
        resize(pma, src);
        index = binary_pma_find(pma, src, des, 0, pma->arr_cap - 1);
        return pma_insert(pma, src, des, index);
    }
    pma->ele_num++;
    src->degree++;
    // print_graph();

    // cout << "test---before rebalance" << endl;
    
    rebalance(pma, src, inserted_index);
    refresh_loaded_block_meta(src);

    return inserted_index;
}

bool pma_delete(PMA *pma, Vertex *src, int32_t des)
{
    if (pma == nullptr || src->blk_list.empty() || pma->ele_num == 0) {
        return false;
    }

    int32_t index = pma_find_index_scan(src, des);
    if (index < 0) {
        return false;
    }

    int32_t blk_idx = index / PM_BLK_SIZE;
    int32_t bit_idx = index % PM_BLK_SIZE;
    if (src->blk_list[blk_idx].buf == nullptr) {
        cacheBlock(src, blk_idx);
    }
    src->blk_list[blk_idx].buf[bit_idx] = -1;
    pma->ele_num--;
    src->degree--;
    if (global_delete_compaction_enabled) {
        pma_compact(pma, src);
    } else {
        refresh_pm_blk_meta(src->blk_list[blk_idx], src->blk_list[blk_idx].buf);
    }
    return true;
}

void pma_compact(PMA *pma, Vertex *src)
{
    if (pma == nullptr || src->blk_list.empty()) {
        return;
    }

    vector<int32_t> values;
    values.reserve(pma->ele_num);
    for (int32_t idx = 0; idx < pma->arr_cap; idx++) {
        int32_t blk_idx = idx / PM_BLK_SIZE;
        int32_t bit_idx = idx % PM_BLK_SIZE;
        if (src->blk_list[blk_idx].buf == nullptr) {
            cacheBlock(src, blk_idx);
        }
        if (src->blk_list[blk_idx].buf[bit_idx] >= 0) {
            values.push_back(src->blk_list[blk_idx].buf[bit_idx]);
        }
    }

    for (int32_t blk_idx = 0; blk_idx < static_cast<int32_t>(src->blk_list.size()); blk_idx++) {
        if (src->blk_list[blk_idx].buf == nullptr) {
            cacheBlock(src, blk_idx);
        }
        fill(src->blk_list[blk_idx].buf, src->blk_list[blk_idx].buf + PM_BLK_SIZE, -1);
    }

    if (values.empty()) {
        refresh_loaded_block_meta(src);
        return;
    }

    int32_t capacity = pma->arr_cap;
    int32_t frequency = (capacity << 8) / static_cast<int32_t>(values.size());
    int32_t write_index = (capacity << 8) - frequency;
    for (int32_t ridx = static_cast<int32_t>(values.size()) - 1; ridx >= 0; ridx--) {
        int32_t target = write_index >> 8;
        int32_t blk_idx = target / PM_BLK_SIZE;
        int32_t bit_idx = target % PM_BLK_SIZE;
        src->blk_list[blk_idx].buf[bit_idx] = values[ridx];
        write_index -= frequency;
    }
    refresh_loaded_block_meta(src);
}


// 逐层扩大窗口，查找满足密度阈值的最小窗口。
void rebalance(PMA* pma, Vertex *src, int32_t index) 
{
    // 通过统计窗口内有效元素密度，决定执行 spread 还是 resize
    int32_t window_start, window_end;
    int32_t height = pma_block_base_height(pma);
    uint32_t occupancy = 0;
    double density, t_height, p_height;
    
    do {
        uint32_t window_size = pma->seg_len * (1 << height);
        uint32_t window = index / window_size;
        window_start = window * window_size;
        window_end = window_start + window_size;
        occupancy = count_window_occupancy(src, window_start, window_end);
        density = (double)occupancy / (double)window_size;
        t_height = t_0 - (height * pma->delta_t);
        p_height = p_0 + (height * pma->delta_p);
        height++;
    } while ((density < p_height ||  // 若窗口密度越界，则继续向上扩大窗口
            density >= t_height) &&
           height < pma->tree_h);
    if (density >= p_height && density < t_height) {
        // 找到了密度合法的窗口
        if (window_end - window_start <= PM_BLK_SIZE) {
            // 窗口就是 block 本身，元素已插入，无需在块内重排
            return;
        }
        spread(pma, src, window_start, window_end, occupancy);
    } else {
        // 所有层级密度均越界（或 tree_h 过小无法扩窗），扩容重建
        resize(pma, src);
    }
}


// 将窗口中的有效元素在 DRAM 中重新均匀铺开。
void spread (PMA* pma, Vertex *src, int32_t from, int32_t to, int32_t occupancy)
{
    /*
        将窗口中的有效元素重新均匀铺开。
        这一步只在 DRAM 中重排，后续再由刷盘逻辑写回设备。
    */
    assert (from < to);
    int32_t capacity = to - from; // 右端点 to 不包含
    int32_t frequency = (capacity << 8) / occupancy;  /* 8-bit fixed point arithmetic. */
    int32_t read_index = to - 1;  // 从窗口末尾开始扫描
    int32_t write_index = (to << 8) - frequency;

    int32_t start_blk = from / PM_BLK_SIZE;
    int32_t end_blk = (to - 1) / PM_BLK_SIZE;

    int32_t span_blks = end_blk - start_blk + 1;

    int32_t *tmp_win_arr = (int32_t*)malloc(capacity * sizeof(int32_t));

    if (tmp_win_arr != NULL) {
    // 先把临时窗口初始化为空洞
        for (size_t i = 0; i < capacity; ++i) {
            tmp_win_arr[i] = -1;
        }
    }

    if (span_blks == 1) { // 当前窗口完全落在同一个块内
        int32_t base_off = start_blk * PM_BLK_SIZE;
        read_index = read_index - base_off;
        write_index = ((to-base_off) << 8) - frequency;  // 转成块内相对索引
        from = from - base_off;

        if (src->blk_list[start_blk].buf == nullptr) {
            cacheBlock(src, start_blk);  // 先把对应块换入 DRAM
        }

        // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[start_blk].global_wp);
        // invalidBlock(src->blk_list[start_blk].global_wp);

        int32_t real_write_idx = (write_index >> 8) - from;
        occupancy -= 1; // read_index 初始位置单独处理
        while (read_index >= from) {
            if (src->blk_list[start_blk].buf[read_index] >= 0) {
                assert(real_write_idx >= 0);
                tmp_win_arr[real_write_idx] = src->blk_list[start_blk].buf[read_index];

                write_index -= frequency;
                read_index--;
                occupancy --;
                real_write_idx = (write_index >> 8) - from;
                if (real_write_idx < occupancy) {
                    real_write_idx = occupancy;
                }

            } else {
                read_index --;
            }
        }
        // 将均匀铺开的结果写回块缓存
        memcpy(&src->blk_list[start_blk].buf[from], tmp_win_arr, sizeof(int32_t) * capacity);
        refresh_pm_blk_meta(src->blk_list[start_blk], src->blk_list[start_blk].buf);
        // 下面是旧版本中“立即写回设备”的思路，当前改为统一刷盘
        // if (graph_blk_que.size() == 0) {
        //     generate_block_baseline(100);
        // }
        // pm_blk new_blk = graph_blk_que.front();
        // graph_blk_que.pop();
        // writeBlockToZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, new_blk.global_wp);
        // new_blk.src = src->id;
        // new_blk.blk_idx = start_blk;
        // // 更新顶点侧的物理地址元数据
        // src->blk_list[start_blk] = new_blk;
        // validBlock(new_blk.global_wp);
        // 旧块在写回后会失效

    } else {  // 当前窗口跨越多个块，需要逐块处理
        occupancy -= 1; // read_index 初始位置单独处理
        int32_t real_write_idx = (write_index >> 8) - from;
        for (int i = end_blk; i >= start_blk; i--) {
            int32_t loc_win_from, loc_win_to; // 当前块在窗口中的局部范围

            if (src->blk_list[i].buf == nullptr) {
                cacheBlock(src, i);  // 先把对应块换入 DRAM
            }

            // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[i].global_wp);

            if (i > start_blk) {
                loc_win_from = 0;
            } else {
                loc_win_from = from - i * PM_BLK_SIZE;
            }
            if (i < end_blk) {
                loc_win_to = PM_BLK_SIZE - 1;
            } else {
                loc_win_to = to - i * PM_BLK_SIZE - 1;
            }
            // cout << "neighbor: ";
            while (loc_win_to >= loc_win_from) {
                if (src->blk_list[i].buf[loc_win_to] >= 0) {
                    assert(occupancy >= 0);
                    tmp_win_arr[real_write_idx] = src->blk_list[i].buf[loc_win_to];
                    write_index -= frequency;
                    real_write_idx = (write_index >> 8) - from;
                    occupancy --;

                    if (real_write_idx < occupancy) {
                        real_write_idx = occupancy;
                    }
                }
                loc_win_to--;
            }
            // assert((write_index >> 8) >= from);
        }

        int32_t cpy_src = 0;
        int32_t cpy_des = from % PM_BLK_SIZE;
        int32_t cpy_size = PM_BLK_SIZE - cpy_des;  // 当前块需要回写的元素数
        // cout << "start blk: " << start_blk << " end blk: " << end_blk << endl; 
        for (int32_t i = start_blk; i <= end_blk; i++) {
            // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[i].global_wp);
            // invalidBlock(src->blk_list[i].global_wp);
            
            memcpy(&src->blk_list[i].buf[cpy_des], &tmp_win_arr[cpy_src], cpy_size * sizeof(int32_t));

            // 下面是旧版本中“立即写回设备”的思路，当前改为统一刷盘
            // if (graph_blk_que.size() == 0) {
            //     generate_block_baseline(100);
            // }
            // pm_blk new_blk = graph_blk_que.front();
            // graph_blk_que.pop();
            // writeBlockToZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, new_blk.global_wp);
            // new_blk.src = src->id;
            // new_blk.blk_idx = i;
            // // 更新顶点侧的物理地址元数据
            // src->blk_list[i] = new_blk;
            // validBlock(new_blk.global_wp);
            
            assert((cpy_des+cpy_size-1)<PM_BLK_SIZE);
            
            cpy_des = 0;
            cpy_src += cpy_size;

            if (i == (end_blk-1)) { // 为最后一个块准备回写长度
                cpy_size = ((to-1) % PM_BLK_SIZE) + 1;  // 右端点是闭区间索引，因此需要 +1
            } else {
                cpy_size = PM_BLK_SIZE;
            }
        }
        refresh_block_meta_range(src, start_blk, end_blk);
    }
    free(tmp_win_arr);
}

// PMA 扩容并按新的稀疏布局重排元素。
void resize (PMA* pma, Vertex *src)
{
    int32_t old_arr_cap = pma->arr_cap;
    compute_capacity(pma);
    pma->tree_h = floor_lg (pma->seg_num) + 1;
    pma->delta_t = (t_0 - t_h) / pma->tree_h;
    pma->delta_p = (p_h - p_0) / pma->tree_h;

    vector<int32_t> values;
    values.reserve(pma->ele_num);
    for (int32_t idx = 0; idx < old_arr_cap; idx++) {
        int32_t blk_idx = idx / PM_BLK_SIZE;
        int32_t bit_idx = idx % PM_BLK_SIZE;
        if (src->blk_list[blk_idx].buf == nullptr) {
            cacheBlock(src, blk_idx);
        }
        if (src->blk_list[blk_idx].buf[bit_idx] >= 0) {
            values.push_back(src->blk_list[blk_idx].buf[bit_idx]);
        }
    }

    // 根据新容量计算所需块数，并补足新增块
    int32_t blk_num = ceil((double)pma->arr_cap / PM_BLK_SIZE);  // 新 PMA 需要的块数
    int32_t ex_blk_num = blk_num - src->blk_list.size();
    

    for (int32_t i=0; i < ex_blk_num; i++) {
        pm_blk new_blk{};
        new_blk.buf = new int32_t[PM_BLK_SIZE];
        fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
        src->blk_list.push_back(new_blk);
        global_cache_blk++;
    }

    for (int32_t blk_idx = 0; blk_idx < src->blk_list.size(); blk_idx++) {
        if (src->blk_list[blk_idx].buf == nullptr) {
            cacheBlock(src, blk_idx);
        }
        fill(src->blk_list[blk_idx].buf, src->blk_list[blk_idx].buf + PM_BLK_SIZE, -1);
    }

    if (values.empty()) {
        refresh_loaded_block_meta(src);
        return;
    }

    int32_t capacity = pma->arr_cap;
    int32_t frequency = (capacity << 8) / static_cast<int32_t>(values.size());
    int32_t write_index = (capacity << 8) - frequency;

    for (int32_t ridx = static_cast<int32_t>(values.size()) - 1; ridx >= 0; ridx--) {
        int32_t target = write_index >> 8;
        int32_t blk_idx = target / PM_BLK_SIZE;
        int32_t bit_idx = target % PM_BLK_SIZE;
        src->blk_list[blk_idx].buf[bit_idx] = values[ridx];
        write_index -= frequency;
    }
    refresh_loaded_block_meta(src);
}

// 读取边文件。当前按无向图处理，因此每条边会被双向展开。
int load_graph(string path)
{
    std::ifstream inputFile(path);

    if (!inputFile) {
        std::cerr << "Failed to open dataset: " << path << std::endl;
        return -1;
    }
    cout << "Load dataset " << path << endl;

    global_edge_vec.clear();
    global_edge_num = 0;
    global_max_vid = 0;

    int32_t src, des;
    while (inputFile >> src >> des) {
        // std::cout << "Read numbers: " << src << " " << des << std::endl;
        new_edge e1, e2;
        e1.src = src;
        e1.des = des;
        global_edge_vec.push_back(e1);
        global_edge_num++;
        global_max_vid = max(src, global_max_vid);

        e2.src = des;
        e2.des = src;
        global_edge_vec.push_back(e2);
        global_edge_num++;
        global_max_vid = max(des, global_max_vid);
    }


    inputFile.close();
    return 0;
}

// 只扫描边文件元数据。parallel2 流式加载用它拿到顶点上界和边数，
// 避免把全部 directed edges 常驻到 global_edge_vec。
int scan_graph_metadata(string path)
{
    std::ifstream inputFile(path);

    if (!inputFile) {
        std::cerr << "Failed to open dataset: " << path << std::endl;
        return -1;
    }
    cout << "Load dataset " << path << endl;

    global_edge_vec.clear();
    global_edge_vec.shrink_to_fit();
    global_edge_num = 0;
    global_max_vid = 0;

    int32_t src, des;
    while (inputFile >> src >> des) {
        global_edge_num += 2;
        global_max_vid = max(global_max_vid, max(src, des));
    }

    inputFile.close();
    return 0;
}

// 基线初始化：每个顶点预分配一个块缓存，所有插入直接进入 PMA。
void graph_init(int32_t max_vid)
{
	generate_block_baseline(9000);
	global_vectex_vec.resize(max_vid + 1);

    // cout << "test ---0" << endl;
	for (int32_t vidx = 0; vidx <= max_vid; vidx ++) {
        // cout << "vidx: " << vidx << " global_blk_vec size:" << global_blk_vec.size() << endl;
		Vertex* v = new Vertex(vidx);
        // 每个顶点先分配一个块缓存，作为初始邻接存储
        pm_blk new_blk{};
        new_blk.buf = new int32_t[PM_BLK_SIZE]; // 在 DRAM 中分配内存
        fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
        v->blk_list.push_back(new_blk);

		v->pma = pma_creat();
        global_vectex_vec[vidx] = v;
        global_cache_blk++;
	}
    cout << "graph initial complete! total cache blk: " << global_cache_blk << endl;
}

// 改进初始化：小度顶点优先使用 cache_buf，不预分配块缓存。
void graph_init_new(int32_t max_vid)
{
    if (global_zone_class_enabled) {
        generate_block_opt();
    } else {
	    generate_block_baseline((zns.nr_blocks_perzone * global_num_queues));
    }
	global_vectex_vec.resize(max_vid + 1);

    // cout << "test ---0" << endl;
	for (int32_t vidx = 0; vidx <= max_vid; vidx ++) {
        // cout << "vidx: " << vidx << " global_blk_vec size:" << global_blk_vec.size() << endl;
		Vertex* v = new Vertex(vidx);
		v->pma = pma_creat();
        global_vectex_vec[vidx] = v;
        // global_cache_blk++;
	}
    cout << "graph initial complete! total cache blk: " << global_cache_blk << endl;
}

// 并行初始化：按 zone 建立独立的可写块队列。
void graph_init_new2(int32_t max_vid)
{
	generate_block_opt();  // 对每个 zone 的块进行独立管理，便于并行刷盘
	global_vectex_vec.resize(max_vid + 1);

    // cout << "test ---0" << endl;
	for (int32_t vidx = 0; vidx <= max_vid; vidx ++) {
        // cout << "vidx: " << vidx << " global_blk_vec size:" << global_blk_vec.size() << endl;
		Vertex* v = new Vertex(vidx);
		v->pma = pma_creat();
        global_vectex_vec[vidx] = v;
        // global_cache_blk++;
	}
    cout << "graph initial complete! total cache blk: " << global_cache_blk << endl;
}

// 基线构图路径：所有边都直接插入 PMA。
void graph_maintenance_baseline(void)
{
    graph_init(global_max_vid);
    cout << "Load graph on NVM baseline begin!" << endl;
	double start_t = get_current_time();

	for (int64_t edgeIdx = 0; edgeIdx < global_edge_num; edgeIdx++)
	{
        new_edge e = global_edge_vec[edgeIdx];  // 从边数组中取出当前边
		int32_t src = e.src;
		int32_t des = e.des;
		Vertex* v = global_vectex_vec[src];
		int32_t index;

        // cout << "insert edge " << edgeIdx << ": (" << src << ", " << des << ")" << endl; 
        // if (v->blk_list[0].buf != nullptr) 
        // {
        //     cout << "error---0" << endl;
        // }
        index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1); // insert after index
        // print_graph();
        // cout << "index: " << index << endl;

        // cout << "test --0" << endl;
        pma_insert(v->pma, v, des, index);
        // cout << "vid: " << v->id << " degree: " << v->degree << endl;
        
        // cout << "testing-- global_cache_blk: " << global_cache_blk << endl;
        if (global_cache_blk >= cache_blk_threshold) {
            // cout << "dump to zone" << endl;
            dumpToZone();
        }
    }
    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);
    
	double end_t = get_current_time();
    double total_t = end_t - start_t;
    double compute_t = total_t - total_dump_time - total_gc_time;
    if (compute_t < 0.0) {
        compute_t = 0.0;
    }
    last_graph_maintenance_total_time = total_t;
    last_graph_maintenance_compute_time = compute_t;
    // print_graph();
    // print_pma(global_vectex_vec[1]);
    cout << "Load graph total time: " << total_t << endl;
    cout << "Load graph dump count: " << total_dump_count << endl;
    cout << "Load graph dump total time: " << total_dump_time << endl;
    cout << "Load graph GC count: " << total_gc_count << endl;
    cout << "Load graph GC total time: " << total_gc_time << endl;
    cout << "Load graph GC moved blocks: " << total_gc_moved_blocks << endl;
    print_atz_gc_summary();
	cout << "Load graph on NVM baseline time: " << compute_t << endl;
    // print_graph();
	}

bool cache_insert_unique(Vertex *v, int32_t des)
{
    for (int32_t i = 0; i < v->degree; i++) {
        if (v->cache_buf[i] == des) {
            return false;
        }
    }
    v->cache_buf[v->degree] = des;
    v->degree++;
    return true;
}

bool cache_delete(Vertex *v, int32_t des)
{
    for (int32_t i = 0; i < v->degree; i++) {
        if (v->cache_buf[i] == des) {
            for (int32_t j = i + 1; j < v->degree; j++) {
                v->cache_buf[j - 1] = v->cache_buf[j];
            }
            v->degree--;
            v->cache_buf[v->degree] = -1;
            return true;
        }
    }
    return false;
}

static bool cache_contains(Vertex *v, int32_t des)
{
    for (int32_t i = 0; i < v->degree; i++) {
        if (v->cache_buf[i] == des) {
            return true;
        }
    }
    return false;
}

static bool graph_insert_directed_edge_locked(int32_t src, int32_t des)
{
    if (src < 0 || src >= static_cast<int32_t>(global_vectex_vec.size())) {
        return false;
    }

    Vertex *v = global_vectex_vec[src];
    unique_lock<std::shared_timed_mutex> lock(v->mutex);
    if (v->blk_list.empty()) {
        if (cache_contains(v, des)) {
            return false;
        }
        if (v->degree < global_cache_buf_threshold) {
            appendZNSRedoLog(EDGE_LOG_INSERT, src, des);
            bool inserted = cache_insert_unique(v, des);
            if (inserted) {
                note_vertex_update(v);
            }
            return inserted;
        }

        int32_t cached_degree = v->degree;
        vector<int32_t> values;
        values.reserve(cached_degree + 1);
        for (int32_t i = 0; i < cached_degree; i++) {
            values.push_back(v->cache_buf[i]);
        }
        fill(v->cache_buf, v->cache_buf + CACHE_BLK_SIZE, -1);
        values.push_back(des);
        sort(values.begin(), values.end());
        appendZNSRedoLog(EDGE_LOG_INSERT, src, des);
        rebuild_pma_from_sorted_values(v, values);
        note_vertex_update(v);
        return true;
    }

    if (pending_insert_contains(v, des)) {
        return false;
    }
    if (pma_contains_scan(v, des)) {
        return false;
    }
    appendZNSRedoLog(EDGE_LOG_INSERT, src, des);
    int32_t target_blk = -1;
    if (find_target_block_for_insert(v, des, target_blk)) {
        v->pending_insert_by_block[target_blk].push_back(des);
        pm_blk &blk = v->blk_list[target_blk];
        blk.min_val = min(blk.min_val, des);
        blk.max_val = max(blk.max_val, des);
    } else {
        v->pending_insert_unplaced.push_back(des);
    }
    global_pending_edge_count++;
    v->degree++;
    note_vertex_update(v);
    return true;
}

bool graph_insert_directed_edge(int32_t src, int32_t des)
{
    shared_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    bool inserted = graph_insert_directed_edge_locked(src, des);
    if (inserted) {
        maybe_merge_pending_insert_buffers();
    }
    return inserted;
}

bool graph_insert_undirected_edge(int32_t src, int32_t des)
{
    shared_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    bool forward_inserted = graph_insert_directed_edge_locked(src, des);
    bool backward_inserted = graph_insert_directed_edge_locked(des, src);
    if (forward_inserted || backward_inserted) {
        maybe_merge_pending_insert_buffers();
    }
    return forward_inserted || backward_inserted;
}

static bool graph_delete_directed_edge_locked(int32_t src, int32_t des)
{
    if (src < 0 || src >= static_cast<int32_t>(global_vectex_vec.size())) {
        return false;
    }
    Vertex *v = global_vectex_vec[src];
    unique_lock<std::shared_timed_mutex> lock(v->mutex);
    bool deleted = false;
    if (v->blk_list.empty()) {
        if (!cache_contains(v, des)) {
            return false;
        }
        appendZNSRedoLog(EDGE_LOG_DELETE, src, des);
        deleted = cache_delete(v, des);
    } else {
        if (pending_insert_remove(v, des)) {
            appendZNSRedoLog(EDGE_LOG_DELETE, src, des);
            v->degree--;
            deleted = true;
        } else if (!pma_contains_scan(v, des)) {
            return false;
        } else {
            appendZNSRedoLog(EDGE_LOG_DELETE, src, des);
            deleted = pma_delete(v->pma, v, des);
        }
    }
    if (deleted) {
        note_vertex_update(v);
    }
    return deleted;
}

bool graph_delete_directed_edge(int32_t src, int32_t des)
{
    shared_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    return graph_delete_directed_edge_locked(src, des);
}

bool graph_delete_undirected_edge(int32_t src, int32_t des)
{
    shared_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    bool forward_deleted = graph_delete_directed_edge_locked(src, des);
    bool backward_deleted = graph_delete_directed_edge_locked(des, src);
    return forward_deleted || backward_deleted;
}

VertexClass classify_vertex(Vertex *v)
{
    if (v->degree >= hub_degree_threshold ||
        static_cast<int32_t>(v->blk_list.size()) >= hub_block_threshold) {
        v->vclass = VERTEX_HUB;
    } else if (v->update_count >= hot_update_threshold) {
        v->vclass = VERTEX_HOT;
    } else {
        v->vclass = VERTEX_COLD;
    }
    return v->vclass;
}

static ATZVertexType classify_atz_vertex(Vertex *v)
{
    if (v == nullptr) {
        return ATZ_VERTEX_COLD;
    }
    if (v->temperature >= atz_hot_threshold) {
        v->atz_type = ATZ_VERTEX_HOT;
    } else if (v->temperature <= atz_cold_threshold) {
        v->atz_type = ATZ_VERTEX_COLD;
    } else {
        v->atz_type = ATZ_VERTEX_WARM;
    }
    return v->atz_type;
}

static void note_vertex_update(Vertex *v)
{
    if (v == nullptr) {
        return;
    }
    v->update_count++;
    v->recent_updates++;
    classify_vertex(v);
}

static const char *atz_vertex_type_name(ATZVertexType type)
{
    switch (type) {
        case ATZ_VERTEX_HOT:
            return "hot";
        case ATZ_VERTEX_WARM:
            return "warm";
        case ATZ_VERTEX_COLD:
        default:
            return "cold";
    }
}

static const char *atz_zone_type_name(ATZZoneType type)
{
    switch (type) {
        case ATZ_ZONE_HOT:
            return "hot";
        case ATZ_ZONE_WARM:
            return "warm";
        case ATZ_ZONE_COLD:
            return "cold";
        case ATZ_ZONE_MIXED:
            return "mixed";
        case ATZ_ZONE_FREE:
        default:
            return "free";
    }
}

static void update_vertex_temperatures(void)
{
    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr) {
            continue;
        }
        v->temperature = atz_temperature_alpha * static_cast<double>(v->recent_updates) +
                         (1.0 - atz_temperature_alpha) * v->temperature;
        v->recent_updates = 0;
        classify_atz_vertex(v);
    }
    atz_temperature_initialized = true;
}

static void print_atz_gc_summary(void)
{
    if (!global_atz_gc_enabled) {
        return;
    }
    double write_amp = atz_user_write_bytes > 0 ?
        static_cast<double>(atz_user_write_bytes + atz_gc_write_bytes) /
        static_cast<double>(atz_user_write_bytes) :
        0.0;
    cout << "ATZ GC zone counts hot/warm/mixed/cold/emergency: "
         << atz_hot_zone_gc_count << "/"
         << atz_warm_zone_gc_count << "/"
         << atz_mixed_zone_gc_count << "/"
         << atz_cold_zone_gc_count << "/"
         << atz_emergency_gc_count << endl;
    cout << "ATZ GC migrated hot/warm/cold blocks: "
         << atz_migrated_hot_blocks << "/"
         << atz_migrated_warm_blocks << "/"
         << atz_migrated_cold_blocks << endl;
    cout << "ATZ GC cold skipped/protected: "
         << atz_cold_zone_skipped_count << "/"
         << atz_cold_zone_protected_count << endl;
    cout << "ATZ GC invalid blocks reclaimed: "
         << atz_invalid_blocks_reclaimed << endl;
    cout << "ATZ write bytes user/gc/write_amp: "
         << atz_user_write_bytes << "/"
         << atz_gc_write_bytes << "/"
         << write_amp << endl;
    cout << "ATZ update writes hot/warm/cold: "
         << atz_hot_block_update_writes << "/"
         << atz_warm_block_update_writes << "/"
         << atz_cold_block_update_writes << endl;
}

const char *vertex_class_name(VertexClass vclass)
{
    switch (vclass) {
        case VERTEX_HOT:
            return "hot";
        case VERTEX_HUB:
            return "hub";
        case VERTEX_COLD:
        default:
            return "cold";
    }
}

static bool requestBlockFromZoneSet(const vector<int32_t> &zone_ids, unsigned long long &wp)
{
    for (int32_t zone_idx : zone_ids) {
        if (zone_idx >= 0 &&
            zone_idx < static_cast<int32_t>(global_blk_queues.size()) &&
            global_blk_queues[zone_idx]->getNextBlock(wp)) {
            return true;
        }
    }
    return false;
}

static bool requestBlockFromZoneSetExcluding(const vector<int32_t> &zone_ids,
                                             int32_t excluded_zone,
                                             unsigned long long &wp)
{
    for (int32_t zone_idx : zone_ids) {
        if (zone_idx == excluded_zone) {
            continue;
        }
        if (zone_idx >= 0 &&
            zone_idx < static_cast<int32_t>(global_blk_queues.size()) &&
            global_blk_queues[zone_idx]->getNextBlock(wp)) {
            return true;
        }
    }
    return false;
}

struct ATZZoneStats {
    int32_t valid_blocks = 0;
    int32_t invalid_blocks = 0;
    int32_t free_blocks = 0;
    int32_t written_blocks = 0;
    int32_t hot_live_blocks = 0;
    int32_t warm_live_blocks = 0;
    int32_t cold_live_blocks = 0;
    double invalid_ratio = 0.0;
};

static void ensure_atz_zone_types(void)
{
    if (global_atz_zone_types.size() != zns.nr_zones) {
        global_atz_zone_types.assign(zns.nr_zones, ATZ_ZONE_FREE);
    }
    for (int32_t zone_idx = 0; zone_idx < static_cast<int32_t>(zns.nr_zones); zone_idx++) {
        if (zone_idx >= global_num_queues) {
            global_atz_zone_types[zone_idx] = ATZ_ZONE_FREE;
        }
    }
}

static ATZZoneStats compute_atz_zone_stats(int32_t zone_idx)
{
    ATZZoneStats stats;
    if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_zone_vec.size())) {
        return stats;
    }
    if (zone_idx < static_cast<int32_t>(global_blk_queues.size())) {
        lock_guard<std::mutex> lock(global_blk_queues[zone_idx]->mtx_);
        stats.free_blocks = static_cast<int32_t>(global_blk_queues[zone_idx]->memoryWPs.size());
    }
    stats.written_blocks = zns.nr_blocks_perzone - stats.free_blocks;
    for (int32_t blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
        const blk_info &info = global_zone_vec[zone_idx]->blk_info_vec[blkidx];
        if (!info.status) {
            continue;
        }
        stats.valid_blocks++;
        Vertex *v = nullptr;
        if (info.src >= 0 && info.src < static_cast<int32_t>(global_vectex_vec.size())) {
            v = global_vectex_vec[info.src];
        }
        ATZVertexType vtype = v == nullptr ? ATZ_VERTEX_COLD : classify_atz_vertex(v);
        if (vtype == ATZ_VERTEX_HOT) {
            stats.hot_live_blocks++;
        } else if (vtype == ATZ_VERTEX_WARM) {
            stats.warm_live_blocks++;
        } else {
            stats.cold_live_blocks++;
        }
    }
    stats.invalid_blocks = max(stats.written_blocks - stats.valid_blocks, 0);
    if (stats.written_blocks > 0) {
        stats.invalid_ratio =
            static_cast<double>(stats.invalid_blocks) /
            static_cast<double>(stats.written_blocks);
    }
    return stats;
}

static vector<int32_t> atz_zones_of_type(ATZZoneType type)
{
    vector<int32_t> zones;
    ensure_atz_zone_types();
    for (int32_t zone_idx = 0; zone_idx < global_num_queues; zone_idx++) {
        if (global_atz_zone_types[zone_idx] == type) {
            zones.push_back(zone_idx);
        }
    }
    return zones;
}

static bool requestATZBlockFromType(ATZZoneType type,
                                    int32_t excluded_zone,
                                    unsigned long long &wp)
{
    ensure_atz_zone_types();
    for (int32_t zone_idx = 0; zone_idx < global_num_queues; zone_idx++) {
        if (zone_idx == excluded_zone) {
            continue;
        }
        if (global_atz_zone_types[zone_idx] != type &&
            global_atz_zone_types[zone_idx] != ATZ_ZONE_FREE) {
            continue;
        }
        if (zone_idx >= static_cast<int32_t>(global_blk_queues.size())) {
            continue;
        }
        if (global_blk_queues[zone_idx]->getNextBlock(wp)) {
            if (global_atz_zone_types[zone_idx] == ATZ_ZONE_FREE) {
                global_atz_zone_types[zone_idx] = type;
            }
            return true;
        }
    }
    return false;
}

static bool requestATZBlock(Vertex *v,
                            int32_t excluded_zone,
                            bool emergency,
                            unsigned long long &wp)
{
    ATZVertexType vtype = classify_atz_vertex(v);
    ATZZoneType preferred = ATZ_ZONE_MIXED;
    if (!atz_temperature_initialized) {
        preferred = ATZ_ZONE_MIXED;
    } else if (vtype == ATZ_VERTEX_HOT) {
        preferred = ATZ_ZONE_HOT;
    } else if (vtype == ATZ_VERTEX_WARM) {
        preferred = ATZ_ZONE_WARM;
    } else {
        preferred = ATZ_ZONE_COLD;
    }

    if (requestATZBlockFromType(preferred, excluded_zone, wp)) {
        return true;
    }
    if (vtype == ATZ_VERTEX_HOT && !emergency) {
        return requestATZBlockFromType(ATZ_ZONE_WARM, excluded_zone, wp) ||
               requestATZBlockFromType(ATZ_ZONE_MIXED, excluded_zone, wp);
    }
    if (vtype == ATZ_VERTEX_COLD && !emergency) {
        return requestATZBlockFromType(ATZ_ZONE_MIXED, excluded_zone, wp) ||
               requestATZBlockFromType(ATZ_ZONE_WARM, excluded_zone, wp);
    }
    return requestATZBlockFromType(ATZ_ZONE_MIXED, excluded_zone, wp) ||
           requestATZBlockFromType(ATZ_ZONE_WARM, excluded_zone, wp) ||
           requestATZBlockFromType(ATZ_ZONE_HOT, excluded_zone, wp) ||
           requestATZBlockFromType(ATZ_ZONE_COLD, excluded_zone, wp);
}

static int64_t countFreeBlocksInZoneSet(const vector<int32_t> &zone_ids)
{
    int64_t free_blocks = 0;
    for (int32_t zone_idx : zone_ids) {
        if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_blk_queues.size())) {
            continue;
        }
        lock_guard<std::mutex> lock(global_blk_queues[zone_idx]->mtx_);
        free_blocks += static_cast<int64_t>(global_blk_queues[zone_idx]->memoryWPs.size());
    }
    return free_blocks;
}

static int64_t countFreeBlocksAllDataZones(void)
{
    int64_t free_blocks = 0;
    for (int32_t zone_idx = 0; zone_idx < global_num_queues; zone_idx++) {
        if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_blk_queues.size())) {
            continue;
        }
        lock_guard<std::mutex> lock(global_blk_queues[zone_idx]->mtx_);
        free_blocks += static_cast<int64_t>(global_blk_queues[zone_idx]->memoryWPs.size());
    }
    return free_blocks;
}

static void maybe_run_proactive_gc(void)
{
    if (global_gc_proactive_free_blocks <= 0 ||
        global_blk_queues.empty() ||
        global_num_queues <= 0) {
        return;
    }
    int64_t free_blocks = countFreeBlocksAllDataZones();
    if (free_blocks > global_gc_proactive_free_blocks) {
        return;
    }
    cout << "Proactive GC trigger: free_blocks=" << free_blocks
         << " threshold=" << global_gc_proactive_free_blocks << endl;
    garbage_collection_greedy();
}

static bool requestRelocationBlock(Vertex *v, int32_t excluded_zone, unsigned long long &wp)
{
    if (global_atz_gc_enabled) {
        bool emergency = false;
        return requestATZBlock(v, excluded_zone, emergency, wp);
    }

    if (!global_zone_class_enabled) {
        for (int32_t zone_idx = 0; zone_idx < global_num_queues; zone_idx++) {
            if (zone_idx == excluded_zone) {
                continue;
            }
            if (zone_idx >= 0 &&
                zone_idx < static_cast<int32_t>(global_blk_queues.size()) &&
                global_blk_queues[zone_idx]->getNextBlock(wp)) {
                return true;
            }
        }
        return false;
    }

    VertexClass vclass = classify_vertex(v);
    if (vclass == VERTEX_HUB &&
        requestBlockFromZoneSetExcluding(hub_zone_ids, excluded_zone, wp)) {
        return true;
    }
    if (vclass == VERTEX_HOT &&
        requestBlockFromZoneSetExcluding(hot_zone_ids, excluded_zone, wp)) {
        return true;
    }
    if (vclass == VERTEX_COLD &&
        requestBlockFromZoneSetExcluding(cold_zone_ids, excluded_zone, wp)) {
        return true;
    }
    if (requestBlockFromZoneSetExcluding(cold_zone_ids, excluded_zone, wp)) {
        return true;
    }
    if (requestBlockFromZoneSetExcluding(hot_zone_ids, excluded_zone, wp)) {
        return true;
    }
    if (requestBlockFromZoneSetExcluding(hub_zone_ids, excluded_zone, wp)) {
        return true;
    }
    return false;
}

static int32_t countValidBlocksInZone(int32_t zone_idx)
{
    if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_zone_vec.size())) {
        return 0;
    }
    int32_t valid = 0;
    for (int32_t blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
        if (global_zone_vec[zone_idx]->blk_info_vec[blkidx].status) {
            valid++;
        }
    }
    return valid;
}

static int32_t chooseVictimZone(const vector<int32_t> &zone_ids)
{
    int32_t best_zone = -1;
    int32_t best_invalid = -1;
    for (int32_t zone_idx : zone_ids) {
        if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_zone_vec.size())) {
            continue;
        }
        int32_t valid = countValidBlocksInZone(zone_idx);
        int32_t free_blocks = 0;
        if (zone_idx < static_cast<int32_t>(global_blk_queues.size())) {
            lock_guard<std::mutex> lock(global_blk_queues[zone_idx]->mtx_);
            free_blocks = static_cast<int32_t>(global_blk_queues[zone_idx]->memoryWPs.size());
        }
        int32_t written = zns.nr_blocks_perzone - free_blocks;
        int32_t invalid = max(written - valid, 0);
        if (invalid == 0) {
            continue;
        }
        if (invalid > best_invalid) {
            best_invalid = invalid;
            best_zone = zone_idx;
        }
    }
    return best_zone;
}

static int32_t chooseBasicVictimZone(const vector<int32_t> &zone_ids)
{
    for (int32_t zone_idx : zone_ids) {
        if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_zone_vec.size())) {
            continue;
        }
        int32_t valid = countValidBlocksInZone(zone_idx);
        int32_t free_blocks = 0;
        if (zone_idx < static_cast<int32_t>(global_blk_queues.size())) {
            lock_guard<std::mutex> lock(global_blk_queues[zone_idx]->mtx_);
            free_blocks = static_cast<int32_t>(global_blk_queues[zone_idx]->memoryWPs.size());
        }
        int32_t written = zns.nr_blocks_perzone - free_blocks;
        int32_t invalid = max(written - valid, 0);
        if (invalid > 0) {
            return zone_idx;
        }
    }
    return -1;
}

static int32_t select_victim_zone_atz_gc(bool &emergency_out,
                                         ATZZoneStats &victim_stats_out,
                                         bool &cold_protected_out)
{
    ensure_atz_zone_types();
    emergency_out = countFreeBlocksAllDataZones() <= atz_free_zone_emergency_watermark;
    cold_protected_out = false;
    int32_t best_zone = -1;
    double best_score = -1.0;
    const double epsilon = 1.0;

    for (int32_t zone_idx = 0; zone_idx < global_num_queues; zone_idx++) {
        ATZZoneType ztype = global_atz_zone_types[zone_idx];
        if (ztype == ATZ_ZONE_FREE) {
            continue;
        }
        ATZZoneStats stats = compute_atz_zone_stats(zone_idx);
        if (stats.invalid_blocks == 0 && !emergency_out) {
            continue;
        }
        if (ztype == ATZ_ZONE_COLD &&
            !emergency_out &&
            stats.invalid_ratio < atz_cold_gc_threshold) {
            atz_cold_zone_skipped_count++;
            atz_cold_zone_protected_count++;
            cold_protected_out = true;
            continue;
        }

        double denominator =
            static_cast<double>(stats.hot_live_blocks) +
            static_cast<double>(stats.warm_live_blocks) +
            atz_cold_live_beta * static_cast<double>(stats.cold_live_blocks) +
            epsilon;
        double score = static_cast<double>(stats.invalid_blocks) / denominator;
        if (score > best_score) {
            best_score = score;
            best_zone = zone_idx;
            victim_stats_out = stats;
        }
    }
    return best_zone;
}

static void resetZoneMetadataAndQueue(int32_t zone_idx, bool lock_queue = true)
{
    if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_blk_queues.size())) {
        return;
    }
    zbd_reset_zones(zns.dev_fd,
                    static_cast<unsigned long long>(zone_idx) * zns.zone_size,
                    zns.zone_size);
    reset_zns_write_pointer_for_zone(zone_idx);

    ZoneStatus *zone = global_zone_vec[zone_idx];
    zone->val_blocks = 0;
    for (int32_t blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
        zone->blk_info_vec[blkidx] = blk_info();
    }

    unique_lock<std::mutex> queue_lock(global_blk_queues[zone_idx]->mtx_, defer_lock);
    if (lock_queue) {
        queue_lock.lock();
    }
    while (!global_blk_queues[zone_idx]->memoryWPs.empty()) {
        global_blk_queues[zone_idx]->memoryWPs.pop();
    }
    unsigned long long wp = static_cast<unsigned long long>(zone_idx) * zns.zone_size;
    for (int32_t blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
        global_blk_queues[zone_idx]->memoryWPs.push(wp);
        wp += zns_alloc_unit_size();
    }
    if (global_atz_gc_enabled &&
        zone_idx >= 0 &&
        zone_idx < static_cast<int32_t>(global_atz_zone_types.size())) {
        global_atz_zone_types[zone_idx] = ATZ_ZONE_FREE;
    }
}

static bool garbageCollectVictimZone(int32_t victim_zone,
                                     bool emergency_gc = false,
                                     bool cold_protected = false,
                                     const ATZZoneStats *precomputed_stats = nullptr)
{
    if (victim_zone < 0 || victim_zone >= static_cast<int32_t>(global_zone_vec.size())) {
        return false;
    }
    if (victim_zone >= static_cast<int32_t>(global_blk_queues.size())) {
        return false;
    }

    double start_t = get_current_time();
    auto finish_gc = [&](bool result) {
        total_gc_time += (get_current_time() - start_t);
        total_gc_count++;
        return result;
    };
    lock_guard<std::mutex> gc_lock(global_gc_mutex);
    ATZZoneStats victim_stats =
        precomputed_stats == nullptr ?
        compute_atz_zone_stats(victim_zone) :
        *precomputed_stats;
    ATZZoneType victim_type =
        (global_atz_gc_enabled &&
         victim_zone >= 0 &&
         victim_zone < static_cast<int32_t>(global_atz_zone_types.size())) ?
        global_atz_zone_types[victim_zone] :
        ATZ_ZONE_MIXED;
    cout << "GC victim zone: " << victim_zone;
    if (global_atz_gc_enabled) {
        cout << " type=" << atz_zone_type_name(victim_type)
             << " invalid=" << victim_stats.invalid_blocks
             << " valid=" << victim_stats.valid_blocks
             << " hot_live=" << victim_stats.hot_live_blocks
             << " warm_live=" << victim_stats.warm_live_blocks
             << " cold_live=" << victim_stats.cold_live_blocks
             << " emergency=" << (emergency_gc ? 1 : 0)
             << " cold_protected=" << (cold_protected ? 1 : 0);
    }
    cout << endl;
    unique_lock<std::mutex> victim_queue_lock(global_blk_queues[victim_zone]->mtx_);
    int32_t moved = 0;
    int32_t moved_hot = 0;
    int32_t moved_warm = 0;
    int32_t moved_cold = 0;
    int32_t *array = nullptr;
    unsigned long long rp = static_cast<unsigned long long>(victim_zone) * zns.zone_size;
    beginZNSEpoch();
    for (int32_t blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
        blk_info info = global_zone_vec[victim_zone]->blk_info_vec[blkidx];
        if (info.status) {
            Vertex *v = global_vectex_vec[info.src];
            unique_lock<std::shared_timed_mutex> vertex_lock(v->mutex, defer_lock);
            if (!vertex_lock.try_lock()) {
                cout << "GC abort: busy vertex " << info.src << endl;
                commitZNSEpoch();
                delete[] array;
                return finish_gc(false);
            }
            blk_info cur_info = global_zone_vec[victim_zone]->blk_info_vec[blkidx];
            if (!cur_info.status || cur_info.src != info.src ||
                cur_info.blkIndex != info.blkIndex) {
                rp += zns_alloc_unit_size();
                continue;
            }
            unsigned long long new_wp = 0;
            if (!requestRelocationBlock(v, victim_zone, new_wp)) {
                commitZNSEpoch();
                delete[] array;
                return finish_gc(false);
            }

            ATZVertexType vtype = classify_atz_vertex(v);
            readBlockFromZone(zns.dev_fd, array, PM_BLK_SIZE, rp);
            refresh_pm_blk_meta(v->blk_list[info.blkIndex], array);
            if (global_zns_record_enabled) {
                writeZNSRecord(zns.dev_fd, info.src, info.blkIndex, v->pma, array, PM_BLK_SIZE, new_wp);
            } else {
                writeBlockToZone(zns.dev_fd, array, PM_BLK_SIZE, new_wp);
            }
            v->blk_list[info.blkIndex].global_wp = new_wp;
            v->blk_list[info.blkIndex].buf = nullptr;
            global_zone_vec[victim_zone]->blk_info_vec[blkidx].status = false;
            mark_zone(info.src, info.blkIndex, new_wp);
            moved++;
            if (global_atz_gc_enabled) {
                if (vtype == ATZ_VERTEX_HOT) {
                    moved_hot++;
                } else if (vtype == ATZ_VERTEX_WARM) {
                    moved_warm++;
                } else {
                    moved_cold++;
                }
            }
        }
        rp += zns_alloc_unit_size();
    }
    commitZNSEpoch();
    delete[] array;

    resetZoneMetadataAndQueue(victim_zone, false);
    cout << "GC moved blocks: " << moved << endl;
    total_gc_moved_blocks += moved;
    if (global_atz_gc_enabled) {
        atz_gc_write_bytes += static_cast<int64_t>(moved) *
                              static_cast<int64_t>(zns_alloc_unit_size());
        atz_migrated_hot_blocks += moved_hot;
        atz_migrated_warm_blocks += moved_warm;
        atz_migrated_cold_blocks += moved_cold;
        atz_invalid_blocks_reclaimed += victim_stats.invalid_blocks;
        if (emergency_gc) {
            atz_emergency_gc_count++;
        }
        if (victim_type == ATZ_ZONE_HOT) {
            atz_hot_zone_gc_count++;
        } else if (victim_type == ATZ_ZONE_WARM) {
            atz_warm_zone_gc_count++;
        } else if (victim_type == ATZ_ZONE_COLD) {
            atz_cold_zone_gc_count++;
        } else {
            atz_mixed_zone_gc_count++;
        }
        cout << "ATZ GC migrated hot=" << moved_hot
             << " warm=" << moved_warm
             << " cold=" << moved_cold
             << " reclaimed_invalid=" << victim_stats.invalid_blocks << endl;
    }
    return finish_gc(true);
}

bool garbage_collection_greedy(void)
{
    if (global_blk_queues.empty()) {
        return garbage_collection();
    }

    if (global_atz_gc_enabled) {
        update_vertex_temperatures();
        bool emergency = false;
        bool cold_protected = false;
        ATZZoneStats stats;
        int32_t victim = select_victim_zone_atz_gc(emergency, stats, cold_protected);
        return victim >= 0 && garbageCollectVictimZone(victim,
                                                       emergency,
                                                       cold_protected,
                                                       &stats);
    }

    vector<int32_t> candidate_zones;
    if (global_zone_class_enabled) {
        candidate_zones.insert(candidate_zones.end(), hot_zone_ids.begin(), hot_zone_ids.end());
        candidate_zones.insert(candidate_zones.end(), cold_zone_ids.begin(), cold_zone_ids.end());
        candidate_zones.insert(candidate_zones.end(), hub_zone_ids.begin(), hub_zone_ids.end());
    } else {
        for (int32_t zone_idx = 0; zone_idx < global_num_queues; zone_idx++) {
            candidate_zones.push_back(zone_idx);
        }
    }

    int32_t victim = global_basic_gc_enabled ?
        chooseBasicVictimZone(candidate_zones) :
        chooseVictimZone(candidate_zones);
    return victim >= 0 && garbageCollectVictimZone(victim);
}

bool requestClassifiedBlock(Vertex *v, unsigned long long &wp)
{
    if (global_atz_gc_enabled) {
        bool emergency = countFreeBlocksAllDataZones() <= atz_free_zone_emergency_watermark;
        bool ok = requestATZBlock(v, -1, emergency, wp);
        if (ok) {
            ATZVertexType type = classify_atz_vertex(v);
            if (type == ATZ_VERTEX_HOT) {
                atz_hot_block_update_writes++;
            } else if (type == ATZ_VERTEX_WARM) {
                atz_warm_block_update_writes++;
            } else {
                atz_cold_block_update_writes++;
            }
            atz_user_write_bytes += static_cast<int64_t>(zns_alloc_unit_size());
        }
        return ok;
    }

    VertexClass vclass = classify_vertex(v);
    if (vclass == VERTEX_HUB && requestBlockFromZoneSet(hub_zone_ids, wp)) {
        return true;
    }
    if (vclass == VERTEX_HOT && requestBlockFromZoneSet(hot_zone_ids, wp)) {
        return true;
    }
    if (vclass == VERTEX_COLD && requestBlockFromZoneSet(cold_zone_ids, wp)) {
        return true;
    }

    // Fallback across all pools keeps experiments running when one class fills.
    if (requestBlockFromZoneSet(cold_zone_ids, wp)) {
        return true;
    }
    if (requestBlockFromZoneSet(hot_zone_ids, wp)) {
        return true;
    }
    if (requestBlockFromZoneSet(hub_zone_ids, wp)) {
        return true;
    }

    bool resume_epoch = global_zns_record_enabled && global_zns_epoch_active;
    if (resume_epoch) {
        commitZNSEpoch();
    }
    bool gc_ok = garbage_collection_greedy();
    if (resume_epoch) {
        beginZNSEpoch();
    }
    if (gc_ok) {
        if (vclass == VERTEX_HUB && requestBlockFromZoneSet(hub_zone_ids, wp)) {
            return true;
        }
        if (vclass == VERTEX_HOT && requestBlockFromZoneSet(hot_zone_ids, wp)) {
            return true;
        }
        if (vclass == VERTEX_COLD && requestBlockFromZoneSet(cold_zone_ids, wp)) {
            return true;
        }
        if (requestBlockFromZoneSet(cold_zone_ids, wp)) {
            return true;
        }
        if (requestBlockFromZoneSet(hot_zone_ids, wp)) {
            return true;
        }
        if (requestBlockFromZoneSet(hub_zone_ids, wp)) {
            return true;
        }
    }
    return false;
}

// 改进构图路径：小度顶点先进入 cache_buf，再按需迁入 PMA。
void graph_maintenance_new(void)
{
    /*
        与 baseline 的区别在于：
        小度顶点先写入 cache_buf，超过阈值后再迁入 PMA。
    */
    graph_init_new(global_max_vid);
    cout << "Load graph on NVM new begin!" << endl;
	double start_t = get_current_time();

	for (int64_t edgeIdx = 0; edgeIdx < global_edge_num; edgeIdx++)
	{
        new_edge e = global_edge_vec[edgeIdx];  // 从边数组中取出当前边
		int32_t src = e.src;
		int32_t des = e.des;
        // cout << "insert edge " << edgeIdx << ": (" << src << ", " << des << ")" << endl; 

        graph_insert_directed_edge(src, des);
        if (global_cache_blk >= cache_blk_threshold) {
            if (!global_debug_pending_merge) {
                cout << "dump to zone" << endl;
            }
            dumpToZone();
        }
    }

    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);
    
	double end_t = get_current_time();
    double total_t = end_t - start_t;
    double compute_t = total_t - total_dump_time - total_gc_time;
    if (compute_t < 0.0) {
        compute_t = 0.0;
    }
    last_graph_maintenance_total_time = total_t;
    last_graph_maintenance_compute_time = compute_t;
    // print_graph();
    // print_pma(global_vectex_vec[1]);
    cout << "Load graph total time: " << total_t << endl;
    cout << "Load graph dump count: " << total_dump_count << endl;
    cout << "Load graph dump total time: " << total_dump_time << endl;
    cout << "Load graph GC count: " << total_gc_count << endl;
    cout << "Load graph GC total time: " << total_gc_time << endl;
    cout << "Load graph GC moved blocks: " << total_gc_moved_blocks << endl;
    print_atz_gc_summary();
	cout << "Load graph on NVM new time: " << compute_t << endl;
    // print_graph();
}

// 并行插入版本：顶点级锁保证同一顶点更新串行化。
void graph_maintenance_parallel(void)
{
    /*
        并行插入版本：
        小度顶点先写入 cache_buf，batch 结束后统一串行刷盘。
    */
    graph_init_new(global_max_vid);
    cout << "Load graph on NVM new begin!" << endl;
	double start_t = get_current_time();

    // 每个 batch 完成后检查 DRAM 中缓存块数量，达到阈值则刷盘
    omp_set_num_threads(THREAD_NUM);
    int64_t start_e = 0;
    int64_t end_e = min(start_e + global_batch_size, global_edge_num);
    // global_batch_size = global_edge_num;
    while (end_e <= global_edge_num) {
        #pragma omp parallel for
        for (int64_t edgeIdx = start_e; edgeIdx < end_e; edgeIdx++)
        {
            new_edge e = global_edge_vec[edgeIdx];  // 从边数组中取出当前边
            graph_insert_directed_edge(e.src, e.des);
        }

        if (global_cache_blk >= cache_blk_threshold) {
            if (!global_debug_pending_merge) {
                cout << "dump to zone" << endl;
            }
            dumpToZone();
        }

        if (end_e == global_edge_num) {
            break;
        }
        start_e = end_e;
        end_e = min(start_e + global_batch_size, global_edge_num);
    }

	double end_t = get_current_time();
    double total_t = end_t - start_t;
    double compute_t = total_t - total_dump_time - total_gc_time;
    if (compute_t < 0.0) {
        compute_t = 0.0;
    }
    last_graph_maintenance_total_time = total_t;
    last_graph_maintenance_compute_time = compute_t;
    cout << "Load graph total time: " << total_t << endl;
    cout << "Load graph dump count: " << total_dump_count << endl;
    cout << "Load graph dump total time: " << total_dump_time << endl;
    cout << "Load graph GC count: " << total_gc_count << endl;
    cout << "Load graph GC total time: " << total_gc_time << endl;
    cout << "Load graph GC moved blocks: " << total_gc_moved_blocks << endl;
    print_atz_gc_summary();
	cout << "Load graph on NVM new time: " << compute_t << endl;
    // print_graph();
}

// 并行插入 + 并行刷盘版本，是当前 main 默认使用的构图流程。
void graph_maintenance_parallel2(void)
{
    /*
        并行插入 + 并行刷盘版本：
        每个 batch 结束后检查脏 PMA block 数量，达到阈值才刷盘。
    */
    graph_init_new2(global_max_vid);
    cout << "Load graph on NVM new begin!" << endl;
	double start_t = get_current_time();

    // batch 结束后检查缓存块数量，达到阈值则并行刷盘
    omp_set_num_threads(THREAD_NUM);
    int64_t start_e = 0;
    int64_t end_e = min(start_e + global_batch_size, global_edge_num);
    
    while (end_e <= global_edge_num) {
        // cout << "end_e: " << end_e << endl;
        #pragma omp parallel for
        for (int64_t edgeIdx = start_e; edgeIdx < end_e; edgeIdx++)
        {
            new_edge e = global_edge_vec[edgeIdx];  // 从边数组中取出当前边
            graph_insert_directed_edge(e.src, e.des);
        }

        if (global_cache_blk >= cache_blk_threshold) {
            // cout << "dump to zone" << endl;
            dumpToZoneParallel();
        }

        if (end_e == global_edge_num) {
            break;
        }
        start_e = end_e;
        end_e = min(start_e + global_batch_size, global_edge_num);
    }
	double end_t = get_current_time();
    double total_t = end_t - start_t;
    double compute_t = total_t - total_dump_time - total_gc_time;
    if (compute_t < 0.0) {
        compute_t = 0.0;
    }
    last_graph_maintenance_total_time = total_t;
    last_graph_maintenance_compute_time = compute_t;
    cout << "Load graph total time: " << total_t << endl;
    cout << "Load graph dump count: " << total_dump_count << endl;
    cout << "Load graph dump total time: " << total_dump_time << endl;
    cout << "Load graph GC count: " << total_gc_count << endl;
    cout << "Load graph GC total time: " << total_gc_time << endl;
    cout << "Load graph GC moved blocks: " << total_gc_moved_blocks << endl;
    print_atz_gc_summary();
	cout << "Load graph on NVM new time: " << compute_t << endl;
    // print_graph();
}

void graph_maintenance_parallel2_streaming(const string &path)
{
    /*
        parallel2 的流式构图路径：
        文件预扫描只确定 global_max_vid/global_edge_num，正式加载时按 batch
        读入少量 directed edges，插入完成后立即释放 batch。
    */
    graph_init_new2(global_max_vid);
    cout << "Load graph on NVM new begin!" << endl;
	double start_t = get_current_time();

    std::ifstream inputFile(path);
    if (!inputFile) {
        std::cerr << "Failed to open dataset: " << path << std::endl;
        exit(EXIT_FAILURE);
    }

    omp_set_num_threads(THREAD_NUM);
    vector<new_edge> edge_batch;
    edge_batch.reserve(global_batch_size);
    int64_t processed_edges = 0;

    auto flush_batch = [&]() {
        if (edge_batch.empty()) {
            return;
        }

        int64_t batch_size = static_cast<int64_t>(edge_batch.size());
        #pragma omp parallel for
        for (int64_t edgeIdx = 0; edgeIdx < batch_size; edgeIdx++) {
            new_edge e = edge_batch[edgeIdx];
            graph_insert_directed_edge(e.src, e.des);
        }

        processed_edges += batch_size;
        edge_batch.clear();

        if (global_cache_blk >= cache_blk_threshold) {
            dumpToZoneParallel();
        }
    };

    int32_t src, des;
    while (inputFile >> src >> des) {
        new_edge e1;
        e1.src = src;
        e1.des = des;
        edge_batch.push_back(e1);
        if (static_cast<int64_t>(edge_batch.size()) >= global_batch_size) {
            flush_batch();
        }

        new_edge e2;
        e2.src = des;
        e2.des = src;
        edge_batch.push_back(e2);
        if (static_cast<int64_t>(edge_batch.size()) >= global_batch_size) {
            flush_batch();
        }
    }
    flush_batch();
    inputFile.close();

    if (processed_edges != global_edge_num) {
        cerr << "Warning: streamed edge count " << processed_edges
             << " differs from scanned edge count " << global_edge_num << endl;
    }

	double end_t = get_current_time();
    double total_t = end_t - start_t;
    double compute_t = total_t - total_dump_time - total_gc_time;
    if (compute_t < 0.0) {
        compute_t = 0.0;
    }
    last_graph_maintenance_total_time = total_t;
    last_graph_maintenance_compute_time = compute_t;
    cout << "Load graph total time: " << total_t << endl;
    cout << "Load graph dump count: " << total_dump_count << endl;
    cout << "Load graph dump total time: " << total_dump_time << endl;
    cout << "Load graph GC count: " << total_gc_count << endl;
    cout << "Load graph GC total time: " << total_gc_time << endl;
    cout << "Load graph GC moved blocks: " << total_gc_moved_blocks << endl;
    print_atz_gc_summary();
	cout << "Load graph on NVM new time: " << compute_t << endl;
}

void print_vertex(int32_t vid, int64_t &total_edges)
{
    Vertex *v = global_vectex_vec[vid];
    int32_t *blk_ptr = new int32_t[PM_BLK_SIZE];
    if (v->degree == 0) {
        delete[] blk_ptr;
        return;
    }
    // cout << "vid: " << vid << " blknum: " << v->blk_list.size() << " {";

    if (v->blk_list.empty()) {
        total_edges += v->degree;
        delete[] blk_ptr;
        return;
    }
    
    for (int i = 0; i < v->blk_list.size(); i++) {
        if (v->blk_list[i].buf == nullptr) {
            // cout << "test ==" << endl;
            // cout << "read wp: " << v->blk_list[i].global_wp << " blk: " << (v->blk_list[i].global_wp/ZONE_BLK_SIZE)
            // << " status: " << test_block(0, (v->blk_list[i].global_wp/ZONE_BLK_SIZE)) << endl;
            readBlockFromZone(zns.dev_fd, blk_ptr, PM_BLK_SIZE, v->blk_list[i].global_wp); // 从设备读取一个块
            for (int j = 0; j < PM_BLK_SIZE; j++) {
                if (blk_ptr[j] >= 0) {
                    // cout << blk_ptr[j] << ", ";
                    total_edges++;
                }
            }
        } else {
            for (int j = 0; j < PM_BLK_SIZE; j++) {
                if (v->blk_list[i].buf[j] >= 0) {
                    // cout << v->blk_list[i].buf[j] << ", ";
                    total_edges++;
                }
            }
        }
    }
    delete[] blk_ptr;
    // cout << "}" << endl;
}

void print_graph(void)
{
    int64_t total_edges = 0;
    for (int32_t vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
        print_vertex(vidx, total_edges);
    }
    cout << "Total number of edges: " << total_edges << endl;
}

void print_vertex_pma_blocks(int32_t vid)
{
    Vertex *v = global_vectex_vec[vid];
    cout << "Vertex " << vid
         << " degree=" << v->degree
         << " pma_ele=" << v->pma->ele_num
         << " pma_cap=" << v->pma->arr_cap
         << " blk_count=" << v->blk_list.size()
         << " class=" << vertex_class_name(classify_vertex(v))
         << " updates=" << v->update_count << endl;

    if (v->degree == 0) {
        cout << "  empty" << endl;
        return;
    }

    if (v->blk_list.empty()) {
        cout << "  cache_buf: [";
        for (int32_t i = 0; i < v->degree; i++) {
            cout << v->cache_buf[i];
            if (i + 1 < v->degree) {
                cout << ", ";
            }
        }
        cout << "]" << endl;
        return;
    }

    int32_t *tmp_blk_ptr = nullptr;
    for (int32_t blkidx = 0; blkidx < static_cast<int32_t>(v->blk_list.size()); blkidx++) {
        int32_t *blk_ptr = v->blk_list[blkidx].buf;
        bool loaded_from_zone = false;
        if (blk_ptr == nullptr) {
            readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, v->blk_list[blkidx].global_wp);
            blk_ptr = tmp_blk_ptr;
            loaded_from_zone = true;
        }

        cout << "  blk[" << blkidx << "] wp=" << v->blk_list[blkidx].global_wp
             << " source=" << (loaded_from_zone ? "zone" : "dram")
             << " values=[";
        for (int32_t i = 0; i < PM_BLK_SIZE; i++) {
            cout << blk_ptr[i];
            if (i + 1 < PM_BLK_SIZE) {
                cout << ", ";
            }
        }
        cout << "]" << endl;

        if (loaded_from_zone) {
            delete[] tmp_blk_ptr;
            tmp_blk_ptr = nullptr;
        }
    }
}

static bool vertex_has_dirty_block(Vertex *v)
{
    if (v == nullptr) {
        return false;
    }
    for (const pm_blk &blk : v->blk_list) {
        if (blk.buf != nullptr) {
            return true;
        }
    }
    return false;
}

static void collect_vertex_neighbors_for_debug(Vertex *v, vector<int32_t> &neighbors)
{
    neighbors.clear();
    if (v == nullptr || v->degree == 0) {
        return;
    }
    if (v->blk_list.empty()) {
        for (int32_t i = 0; i < v->degree; i++) {
            if (v->cache_buf[i] >= 0) {
                neighbors.push_back(v->cache_buf[i]);
            }
        }
    } else {
        int32_t tmp_blk_buf[PM_BLK_SIZE];
        for (int32_t blkidx = 0; blkidx < static_cast<int32_t>(v->blk_list.size()); blkidx++) {
            int32_t *blk_ptr = v->blk_list[blkidx].buf;
            if (blk_ptr == nullptr) {
                readBlockFromZoneInto(zns.dev_fd, tmp_blk_buf, PM_BLK_SIZE,
                                      v->blk_list[blkidx].global_wp);
                blk_ptr = tmp_blk_buf;
            }
            int32_t local_end = pma_block_end(v->pma, blkidx) - pma_block_start(blkidx);
            for (int32_t i = 0; i < local_end; i++) {
                if (blk_ptr[i] >= 0) {
                    neighbors.push_back(blk_ptr[i]);
                }
            }
        }
        for (const auto &entry : v->pending_insert_by_block) {
            neighbors.insert(neighbors.end(), entry.second.begin(), entry.second.end());
        }
        neighbors.insert(neighbors.end(),
                         v->pending_insert_unplaced.begin(),
                         v->pending_insert_unplaced.end());
    }
    sort(neighbors.begin(), neighbors.end());
    neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
}

static void print_vertex_pending_debug(int32_t vid)
{
    Vertex *v = global_vectex_vec[vid];
    cout << "Vertex " << vid
         << " pending_count=" << pending_insert_count(v)
         << " dirty=" << (vertex_has_dirty_block(v) ? "yes" : "no")
         << endl;

    if (!v->pending_insert_by_block.empty()) {
        cout << "  pending_by_block:";
        vector<pair<int32_t, vector<int32_t>>> ordered(v->pending_insert_by_block.begin(),
                                                       v->pending_insert_by_block.end());
        sort(ordered.begin(), ordered.end(),
             [](const auto &a, const auto &b) {
                 return a.first < b.first;
             });
        for (const auto &entry : ordered) {
            cout << " blk[" << entry.first << "]=[";
            for (size_t i = 0; i < entry.second.size(); i++) {
                cout << entry.second[i];
                if (i + 1 < entry.second.size()) {
                    cout << ", ";
                }
            }
            cout << "]";
        }
        cout << endl;
    }
    if (!v->pending_insert_unplaced.empty()) {
        cout << "  pending_unplaced=[";
        for (size_t i = 0; i < v->pending_insert_unplaced.size(); i++) {
            cout << v->pending_insert_unplaced[i];
            if (i + 1 < v->pending_insert_unplaced.size()) {
                cout << ", ";
            }
        }
        cout << "]" << endl;
    }

    vector<int32_t> neighbors;
    collect_vertex_neighbors_for_debug(v, neighbors);
    cout << "  effective_neighbors=[";
    for (size_t i = 0; i < neighbors.size(); i++) {
        cout << neighbors[i];
        if (i + 1 < neighbors.size()) {
            cout << ", ";
        }
    }
    cout << "]" << endl;
    print_vertex_pma_blocks(vid);
}

static void print_pending_merge_debug_snapshot(const char *title)
{
    cout << "==== " << title
         << " pending_edges=" << global_pending_edge_count.load()
         << " dirty_blocks=" << global_cache_blk.load()
         << " ====" << endl;
    for (int32_t vidx = 0; vidx < static_cast<int32_t>(global_vectex_vec.size()); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        if (v == nullptr) {
            continue;
        }
        if (pending_insert_empty(v) && !vertex_has_dirty_block(v) &&
            (v->degree == 0 || v->blk_list.empty())) {
            continue;
        }
        print_vertex_pending_debug(vidx);
    }
    cout << "==== end " << title << " ====" << endl;
}

void print_all_vertex_pma_blocks(void)
{
    cout << "==== Vertex PMA Blocks Begin ====" << endl;
    for (int32_t vidx = 0; vidx < static_cast<int32_t>(global_vectex_vec.size()); vidx++) {
        print_vertex_pma_blocks(vidx);
    }
    cout << "==== Vertex PMA Blocks End ====" << endl;
}

// // 记录物理块到顶点逻辑块的映射，供恢复和 GC 使用
// void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp)
// {
//     // cout << "zone size: " << zns.info.zone_size << endl;
//     int32_t zone_idx = wp / zns.info.zone_size;
//     int32_t info_idx = (wp % zns.info.zone_size) / ZONE_BLK_SIZE;
//     global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status = 1;
//     global_zone_vec[zone_idx]->blk_info_vec[blk_idx].src = vid;
//     global_zone_vec[zone_idx]->blk_info_vec[blk_idx].blkIndex = blk_idx;
// }

// 当 DRAM 中缓存块达到阈值后，将其写回设备
// 串行刷盘：把 DRAM 中的块缓存写回设备，并更新物理地址映射。
void dumpToZone(void)
{
    // int total_blks = 0;
    shared_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    double start_t = get_current_time();
    int64_t dumped_blocks = 0;
    beginZNSEpoch();
    for (int32_t vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        unique_lock<std::shared_timed_mutex> lock(v->mutex);
        // cout << "dump vid: " << v->id << " degree: " << v->degree << " blks: " << v->blk_list.size() << endl;
        if ((v->degree > 0) && (v->blk_list.size() > 0)) {
            for (int blkidx = 0; blkidx < v->blk_list.size(); blkidx++) {
                if (v->blk_list[blkidx].buf == nullptr) { // 该块当前不在 DRAM 中
                    continue;
                }
                // total_blks++;
                // cout << "dump vid: " << v->id << endl;

                pm_blk new_blk{};
                new_blk.buf = nullptr;
                if (global_zone_class_enabled && !global_blk_queues.empty()) {
                    if (!requestClassifiedBlock(v, new_blk.global_wp)) {
                        std::cerr << "No free classified blocks for dump." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                } else {
                    if (graph_blk_que.size() == 0) {
                        commitZNSEpoch();
                        lock.unlock();
                        bool gc_ok = garbage_collection();
                        lock.lock();
                        beginZNSEpoch();
                        if (!gc_ok || graph_blk_que.empty()) {
                            std::cerr << "No free baseline blocks after GC." << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        if (v->blk_list[blkidx].buf == nullptr) {
                            continue;
                        }
                    }
                    new_blk = graph_blk_que.front();
                    graph_blk_que.pop();
                }

                refresh_pm_blk_meta(v->blk_list[blkidx], v->blk_list[blkidx].buf);
                new_blk.min_val = v->blk_list[blkidx].min_val;
                new_blk.max_val = v->blk_list[blkidx].max_val;
                new_blk.valid_count = v->blk_list[blkidx].valid_count;
                writeGraphBlockToZone(v, blkidx, new_blk.global_wp);
                delete[] v->blk_list[blkidx].buf;
                mark_zone(v->id, blkidx, new_blk.global_wp);  // 更新块映射，供 GC 使用
                // 更新顶点侧的物理地址元数据
                v->blk_list[blkidx] = new_blk;
                // v->blk_list[blkidx].buf = nullptr;
                // v->blk_list[blkidx].global_wp = new_blk.global_wp;

                global_cache_blk--;
                dumped_blocks++;
            }
        }
    }
    commitZNSEpoch();
    double end_t = get_current_time();
    if (!global_debug_pending_merge) {
        cout << "Dump to zone time: " << end_t - start_t << endl;
    }
    if (dumped_blocks > 0) {
        total_dump_time += (end_t - start_t);
        total_dump_count++;
    }
    maybe_run_proactive_gc();
    // cout << "total_blks: " << total_blks << endl;
}

// 并行刷盘：每个线程从可写块队列中申请物理地址并刷回设备。
void dumpToZoneParallel(void)
{
    // ZNS requires each zone to be written strictly at its current write
    // pointer. The old parallel writer could issue writes to the same zone out
    // of order, so the correctness-preserving path uses the zone-serial writer.
    dumpToZoneZoneSerial();
}

struct ZoneWriteBatch {
    bool active = false;
    unsigned long long start_wp = 0;
    vector<char> data;
};

static void appendGraphBlockToZoneBatch(vector<ZoneWriteBatch> &batches,
                                        Vertex *v,
                                        int32_t blkidx,
                                        unsigned long long wp)
{
    int32_t zone_idx = static_cast<int32_t>(wp / static_cast<unsigned long long>(zns.zone_size));
    if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(batches.size())) {
        std::cerr << "Invalid zone for batched dump wp=" << wp << std::endl;
        exit(EXIT_FAILURE);
    }

    ZoneWriteBatch &batch = batches[zone_idx];
    if (!batch.active) {
        batch.active = true;
        batch.start_wp = wp;
    }
    unsigned long long expected_wp =
        batch.start_wp + static_cast<unsigned long long>(batch.data.size());
    if (wp != expected_wp) {
        std::cerr << "Non-contiguous batched dump write: zone=" << zone_idx
                  << " expected_wp=" << expected_wp
                  << " actual_wp=" << wp << std::endl;
        exit(EXIT_FAILURE);
    }

    size_t old_size = batch.data.size();
    batch.data.resize(old_size + zns_alloc_unit_size(), 0);
    if (global_zns_record_enabled) {
        buildZNSRecord(batch.data.data() + old_size,
                       v->id,
                       blkidx,
                       v->pma,
                       v->blk_list[blkidx].buf,
                       PM_BLK_SIZE);
    } else {
        memcpy(batch.data.data() + old_size,
               v->blk_list[blkidx].buf,
               PM_BLK_SIZE * sizeof(int32_t));
    }
}

static void flushZoneWriteBatches(vector<ZoneWriteBatch> &batches)
{
    for (int32_t zone_idx = 0; zone_idx < static_cast<int32_t>(batches.size()); zone_idx++) {
        ZoneWriteBatch &batch = batches[zone_idx];
        if (!batch.active || batch.data.empty()) {
            continue;
        }
        checkedZoneWriteBatch(zns.dev_fd,
                              batch.data.data(),
                              batch.data.size(),
                              batch.start_wp,
                              "Failed to write batched ZNS records");
    }
}

void dumpToZoneZoneSerial(void)
{
    if (global_blk_queues.empty() || global_num_queues <= 0) {
        dumpToZone();
        return;
    }

    shared_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    double start_t = get_current_time();
    int64_t dumped_blocks = 0;
    beginZNSEpoch();
    vector<ZoneWriteBatch> zone_batches(zns.nr_zones);
    int next_zone = 0;
    for (int32_t vidx = 0; vidx < static_cast<int32_t>(global_vectex_vec.size()); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        unique_lock<std::shared_timed_mutex> lock(v->mutex);
        if ((v->degree > 0) && (v->blk_list.size() > 0)) {
            for (int blkidx = 0; blkidx < static_cast<int32_t>(v->blk_list.size()); blkidx++) {
                if (v->blk_list[blkidx].buf == nullptr) {
                    continue;
                }

                unsigned long long wp = 0;
                bool acquired = false;
                if (global_zone_class_enabled) {
                    acquired = requestClassifiedBlock(v, wp);
                } else {
                    for (int tries = 0; tries < global_num_queues; tries++) {
                        int zone_idx = (next_zone + tries) % global_num_queues;
                        if (global_blk_queues[zone_idx]->getNextBlock(wp)) {
                            next_zone = (zone_idx + 1) % global_num_queues;
                            acquired = true;
                            break;
                        }
                    }
                }
                if (!acquired) {
                    std::cerr << "No free blocks for zone-serial dump." << std::endl;
                    exit(EXIT_FAILURE);
                }

                refresh_pm_blk_meta(v->blk_list[blkidx], v->blk_list[blkidx].buf);
                appendGraphBlockToZoneBatch(zone_batches, v, blkidx, wp);
                delete[] v->blk_list[blkidx].buf;
                mark_zone(v->id, blkidx, wp);
                v->blk_list[blkidx].buf = nullptr;
                v->blk_list[blkidx].global_wp = wp;
                global_cache_blk--;
                dumped_blocks++;
            }
        }
    }
    flushZoneWriteBatches(zone_batches);
    commitZNSEpoch();
    double end_t = get_current_time();
    if (!global_debug_pending_merge) {
        cout << "Dump to zone-serial writer time: " << end_t - start_t << endl;
    }
    if (dumped_blocks > 0) {
        total_dump_time += (end_t - start_t);
        total_dump_count++;
    }
    maybe_run_proactive_gc();
}

static bool request_any_record_block(Vertex *v, unsigned long long &wp)
{
    if (global_zone_class_enabled && !global_blk_queues.empty()) {
        return requestClassifiedBlock(v, wp);
    }
    if (!global_blk_queues.empty()) {
        for (int zone_idx = 0; zone_idx < global_num_queues; zone_idx++) {
            if (global_blk_queues[zone_idx]->getNextBlock(wp)) {
                return true;
            }
        }
    }
    if (!graph_blk_que.empty()) {
        pm_blk blk = graph_blk_que.front();
        graph_blk_que.pop();
        wp = blk.global_wp;
        return true;
    }
    bool resume_epoch = global_zns_record_enabled && global_zns_epoch_active;
    if (resume_epoch) {
        commitZNSEpoch();
    }
    garbage_collection_greedy();
    if (resume_epoch) {
        beginZNSEpoch();
    }
    if (!graph_blk_que.empty()) {
        pm_blk blk = graph_blk_que.front();
        graph_blk_que.pop();
        wp = blk.global_wp;
        return true;
    }
    return false;
}

void persist_cache_vertices_as_records(void)
{
    if (!global_zns_record_enabled) {
        return;
    }

    shared_lock<std::shared_timed_mutex> graph_lock(global_query_update_mutex);
    beginZNSEpoch();
    vector<ZoneWriteBatch> zone_batches(zns.nr_zones);
    for (int32_t vidx = 0; vidx < static_cast<int32_t>(global_vectex_vec.size()); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        unique_lock<std::shared_timed_mutex> lock(v->mutex);
        if (v->degree == 0 || !v->blk_list.empty()) {
            continue;
        }

        pm_blk blk{};
        blk.buf = new int32_t[PM_BLK_SIZE];
        fill(blk.buf, blk.buf + PM_BLK_SIZE, -1);
        for (int32_t i = 0; i < v->degree; i++) {
            blk.buf[i] = v->cache_buf[i];
        }
        refresh_pm_blk_meta(blk, blk.buf);
        v->blk_list.push_back(blk);
        v->pma->ele_num = v->degree;
        v->pma->arr_cap = max(v->pma->arr_cap, CACHE_BLK_SIZE);

        unsigned long long wp = 0;
        if (!request_any_record_block(v, wp)) {
            std::cerr << "No free block to persist cache vertex record." << std::endl;
            exit(EXIT_FAILURE);
        }
        appendGraphBlockToZoneBatch(zone_batches, v, 0, wp);
        delete[] v->blk_list[0].buf;
        v->blk_list[0].buf = nullptr;
        v->blk_list[0].global_wp = wp;
        mark_zone(v->id, 0, wp);
    }
    flushZoneWriteBatches(zone_batches);
    commitZNSEpoch();
}

static void print_size_line(const char *name, size_t bytes)
{
    double mib = static_cast<double>(bytes) / (1024.0 * 1024.0);
    cout << name << ": " << bytes << " bytes (" << mib << " MiB)" << endl;
}

void print_dram_metadata_usage(void)
{
    size_t vertex_ptr_array_bytes = global_vectex_vec.capacity() * sizeof(Vertex *);
    size_t vertex_object_bytes = 0;
    size_t pma_object_bytes = 0;
    size_t blk_list_meta_bytes = 0;
    size_t pending_dynamic_bytes = 0;
    size_t cached_block_payload_bytes = 0;
    size_t cache_buf_bytes = 0;
    size_t nonnull_vertices = 0;
    size_t pma_block_meta_count = 0;
    size_t cached_block_count = 0;
    size_t pending_edge_capacity = 0;

    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr) {
            continue;
        }
        nonnull_vertices++;
        vertex_object_bytes += sizeof(Vertex);
        cache_buf_bytes += sizeof(v->cache_buf);
        if (v->pma != nullptr) {
            pma_object_bytes += sizeof(PMA);
        }

        blk_list_meta_bytes += v->blk_list.capacity() * sizeof(pm_blk);
        pma_block_meta_count += v->blk_list.size();
        for (const pm_blk &blk : v->blk_list) {
            if (blk.buf != nullptr) {
                cached_block_count++;
                cached_block_payload_bytes += PM_BLK_SIZE * sizeof(int32_t);
            }
        }

        pending_dynamic_bytes += v->pending_insert_unplaced.capacity() * sizeof(int32_t);
        pending_edge_capacity += v->pending_insert_unplaced.capacity();
        pending_dynamic_bytes += v->pending_insert_by_block.bucket_count() * sizeof(void *);
        for (const auto &entry : v->pending_insert_by_block) {
            pending_dynamic_bytes += sizeof(entry);
            pending_dynamic_bytes += entry.second.capacity() * sizeof(int32_t);
            pending_edge_capacity += entry.second.capacity();
        }
    }

    size_t zone_metadata_bytes = global_zone_vec.capacity() * sizeof(ZoneStatus *);
    for (ZoneStatus *zone : global_zone_vec) {
        if (zone == nullptr) {
            continue;
        }
        zone_metadata_bytes += sizeof(ZoneStatus);
        zone_metadata_bytes += zone->blk_info_vec.capacity() * sizeof(blk_info);
    }
    zone_metadata_bytes += global_zone_write_wp.capacity() * sizeof(unsigned long long);
    zone_metadata_bytes += global_zone_write_limit.capacity() * sizeof(unsigned long long);
    zone_metadata_bytes += hot_zone_ids.capacity() * sizeof(int32_t);
    zone_metadata_bytes += cold_zone_ids.capacity() * sizeof(int32_t);
    zone_metadata_bytes += hub_zone_ids.capacity() * sizeof(int32_t);

    size_t free_wp_queue_payload_bytes = global_blk_queues.capacity() * sizeof(MemoryWPQueue *);
    for (MemoryWPQueue *queue : global_blk_queues) {
        if (queue == nullptr) {
            continue;
        }
        free_wp_queue_payload_bytes += sizeof(MemoryWPQueue);
        free_wp_queue_payload_bytes += queue->memoryWPs.size() * sizeof(unsigned long long);
    }

    size_t input_edge_vec_bytes = global_edge_vec.capacity() * sizeof(new_edge);
    size_t graph_metadata_bytes = vertex_ptr_array_bytes +
                                  vertex_object_bytes +
                                  pma_object_bytes +
                                  blk_list_meta_bytes +
                                  pending_dynamic_bytes +
                                  zone_metadata_bytes +
                                  free_wp_queue_payload_bytes;

    size_t graph_metadata_excluding_cache_buf =
        graph_metadata_bytes >= cache_buf_bytes ?
        graph_metadata_bytes - cache_buf_bytes :
        graph_metadata_bytes;

    cout << "==== DRAM metadata usage estimate ====" << endl;
    cout << "nonnull vertices: " << nonnull_vertices << endl;
    cout << "PMA block metadata entries: " << pma_block_meta_count << endl;
    cout << "cached PMA block payload count excluded: " << cached_block_count << endl;
    cout << "pending edge capacity: " << pending_edge_capacity << endl;
    print_size_line("vertex pointer array", vertex_ptr_array_bytes);
    print_size_line("Vertex objects", vertex_object_bytes);
    print_size_line("  cache_buf inside Vertex", cache_buf_bytes);
    print_size_line("PMA objects", pma_object_bytes);
    print_size_line("blk_list metadata arrays", blk_list_meta_bytes);
    print_size_line("pending dynamic containers", pending_dynamic_bytes);
    print_size_line("zone/free-block metadata", zone_metadata_bytes + free_wp_queue_payload_bytes);
    print_size_line("graph metadata total, excluding PMA block payload", graph_metadata_bytes);
    print_size_line("graph metadata total, also excluding cache_buf", graph_metadata_excluding_cache_buf);
    print_size_line("excluded cached PMA block payload", cached_block_payload_bytes);
    print_size_line("input global_edge_vec, not graph metadata", input_edge_vec_bytes);
    cout << "==== end DRAM metadata usage estimate ====" << endl;
}

struct RecoveredRecord {
    ZNSBlockHeader header;
    unsigned long long wp;
};

struct RecoveredEpochStats {
    uint64_t first_seq = 0;
    uint64_t last_seq = 0;
    uint32_t data_record_count = 0;
    uint32_t aggregate_checksum = 0;
};

static bool epoch_stats_match(const ZNSBlockHeader &commit,
                              const RecoveredEpochStats &stats)
{
    return commit.first_seq == stats.first_seq &&
           commit.last_seq == stats.last_seq &&
           commit.data_record_count == stats.data_record_count &&
           commit.aggregate_checksum == stats.aggregate_checksum;
}

static uint64_t recovered_record_key(int32_t vid, int32_t blk_idx)
{
    return (static_cast<uint64_t>(static_cast<uint32_t>(vid)) << 32) |
           static_cast<uint32_t>(blk_idx);
}

static void refresh_recovered_pma_shape(Vertex *v)
{
    if (v == nullptr || v->pma == nullptr) {
        return;
    }
    v->pma->seg_len = BLOCK_ALIGNED_SEG_LEN;
    int32_t min_cap = max(v->pma->seg_len * 2, v->pma->arr_cap);
    v->pma->seg_num = hyperceil(ceil_div(min_cap, v->pma->seg_len));
    v->pma->arr_cap = v->pma->seg_len * v->pma->seg_num;
    v->pma->tree_h = floor_lg(v->pma->seg_num) + 1;
    v->pma->delta_t = (t_0 - t_h) / v->pma->tree_h;
    v->pma->delta_p = (p_h - p_0) / v->pma->tree_h;
}

static uint64_t materialize_missing_recovered_blocks(const unordered_set<uint64_t> &recovered_keys)
{
    uint64_t missing_blocks = 0;
    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr || v->pma == nullptr || v->blk_list.empty()) {
            continue;
        }
        for (int32_t blk_idx = 0; blk_idx < static_cast<int32_t>(v->blk_list.size()); blk_idx++) {
            if (recovered_keys.find(recovered_record_key(v->id, blk_idx)) != recovered_keys.end()) {
                continue;
            }
            if (v->blk_list[blk_idx].buf != nullptr) {
                continue;
            }
            v->blk_list[blk_idx].buf = new int32_t[PM_BLK_SIZE];
            fill(v->blk_list[blk_idx].buf,
                 v->blk_list[blk_idx].buf + PM_BLK_SIZE,
                 -1);
            reset_pm_blk_meta(v->blk_list[blk_idx]);
            global_cache_blk++;
            missing_blocks++;
        }
    }
    return missing_blocks;
}

static void ensure_recovered_vertex_capacity(int32_t vid)
{
    if (vid < 0) {
        return;
    }
    if (static_cast<int32_t>(global_vectex_vec.size()) <= vid) {
        int32_t old_size = static_cast<int32_t>(global_vectex_vec.size());
        global_vectex_vec.resize(vid + 1, nullptr);
        for (int32_t i = old_size; i <= vid; i++) {
            Vertex *v = new Vertex(i);
            v->pma = pma_creat();
            global_vectex_vec[i] = v;
        }
    } else if (global_vectex_vec[vid] == nullptr) {
        Vertex *v = new Vertex(vid);
        v->pma = pma_creat();
        global_vectex_vec[vid] = v;
    }
    global_max_vid = max(global_max_vid, vid);
}

static uint64_t replay_redo_logs_after(uint64_t committed_op_id,
                                       const vector<EdgeLogRecord> &records)
{
    bool prev_replay_state = global_edge_log_replay_active;
    global_edge_log_replay_active = true;
    uint64_t replayed = 0;
    for (const EdgeLogRecord &record : records) {
        if (record.op_id <= committed_op_id) {
            continue;
        }
        ensure_recovered_vertex_capacity(record.src);
        ensure_recovered_vertex_capacity(record.dst);
        if (record.op_type == EDGE_LOG_INSERT) {
            graph_insert_directed_edge_locked(record.src, record.dst);
        } else if (record.op_type == EDGE_LOG_DELETE) {
            graph_delete_directed_edge_locked(record.src, record.dst);
        }
        global_last_applied_op_id = max(global_last_applied_op_id, record.op_id);
        replayed++;
    }
    global_edge_log_replay_active = prev_replay_state;
    cout << "Redo log replayed records: " << replayed
         << " after committed_op_id: " << committed_op_id << endl;
    return replayed;
}

static bool read_redo_log_file_records(vector<EdgeLogRecord> &records,
                                       uint64_t &records_valid,
                                       uint64_t &max_op_id,
                                       int32_t &max_vid)
{
    records_valid = 0;
    if (!global_redo_file_enabled || global_redo_log_file_path.empty()) {
        return false;
    }
    ifstream input(global_redo_log_file_path, ios::binary);
    if (!input.good()) {
        return false;
    }

    input.seekg(0, ios::end);
    streamoff file_size = input.tellg();
    input.seekg(0, ios::beg);
    if (file_size <= 0) {
        return false;
    }

    if (file_size % ZNS_RECORD_SIZE == 0) {
        bool page_file_ok = true;
        while (input.good() && input.tellg() < file_size) {
            char raw[ZNS_RECORD_SIZE];
            input.read(raw, ZNS_RECORD_SIZE);
            if (input.gcount() != ZNS_RECORD_SIZE) {
                page_file_ok = false;
                break;
            }
            vector<EdgeLogRecord> page_records;
            if (!decodeZNSRedoLogPageRaw(raw, page_records)) {
                page_file_ok = false;
                break;
            }
            for (const EdgeLogRecord &record : page_records) {
                records.push_back(record);
                records_valid++;
                max_op_id = max(max_op_id, record.op_id);
                max_vid = max(max_vid, max(record.src, record.dst));
            }
        }
        if (page_file_ok) {
            sort(records.begin(), records.end(),
                 [](const EdgeLogRecord &a, const EdgeLogRecord &b) {
                     return a.op_id < b.op_id;
                 });
            return records_valid > 0;
        }
        records.clear();
        records_valid = 0;
        input.clear();
        input.seekg(0, ios::beg);
    }

    if (file_size % static_cast<streamoff>(sizeof(EdgeLogEntry)) != 0) {
        return false;
    }

    EdgeLogEntry entry;
    while (input.read(reinterpret_cast<char *>(&entry), sizeof(entry))) {
        if (entry.op_type != EDGE_LOG_INSERT &&
            entry.op_type != EDGE_LOG_DELETE) {
            records.clear();
            return false;
        }
        EdgeLogRecord record;
        memset(&record, 0, sizeof(record));
        record.magic = ZNS_REDO_MAGIC;
        record.version = ZNS_BLOCK_VERSION;
        record.op_type = entry.op_type;
        record.op_id = entry.op_id;
        record.src = entry.src;
        record.dst = entry.dst;
        record.tail_magic = ZNS_REDO_TAIL_MAGIC;
        record.checksum = zns_redo_checksum(record);
        records.push_back(record);
        records_valid++;
        max_op_id = max(max_op_id, record.op_id);
        max_vid = max(max_vid, max(record.src, record.dst));
    }
    if (!input.eof()) {
        records.clear();
        records_valid = 0;
        return false;
    }
    sort(records.begin(), records.end(),
         [](const EdgeLogRecord &a, const EdgeLogRecord &b) {
             return a.op_id < b.op_id;
         });
    return records_valid > 0;
}

static unsigned long long recovery_zone_scan_limit(int32_t zone_idx)
{
    unsigned long long zone_start =
        static_cast<unsigned long long>(zone_idx) *
        static_cast<unsigned long long>(zns.zone_size);
    unsigned long long zone_end = zone_start + zns.zone_size;
    if (zone_idx >= 0 &&
        zone_idx < static_cast<int32_t>(global_zone_write_wp.size()) &&
        global_zone_write_wp[zone_idx] >= zone_start) {
        return min(global_zone_write_wp[zone_idx], zone_end);
    }
    return zone_end;
}

struct RecoverySampleVertex {
    int32_t vid = -1;
    vector<int32_t> neighbors;
};

static vector<RecoverySampleVertex> build_recovery_validation_sample(size_t sample_limit)
{
    vector<pair<int32_t, int32_t>> candidates;
    candidates.reserve(global_vectex_vec.size());
    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr || v->degree <= 0) {
            continue;
        }
        candidates.push_back({v->degree, v->id});
    }
    sort(candidates.begin(), candidates.end(), greater<pair<int32_t, int32_t>>());

    vector<int32_t> selected;
    unordered_set<int32_t> selected_set;
    size_t top_count = min(sample_limit / 2, candidates.size());
    for (size_t i = 0; i < top_count; i++) {
        selected.push_back(candidates[i].second);
        selected_set.insert(candidates[i].second);
    }

    if (!candidates.empty()) {
        size_t remaining = sample_limit > selected.size() ?
            sample_limit - selected.size() : 0;
        size_t step = max<size_t>(1, candidates.size() / max<size_t>(remaining, 1));
        for (size_t i = 0; i < candidates.size() && selected.size() < sample_limit; i += step) {
            int32_t vid = candidates[i].second;
            if (selected_set.insert(vid).second) {
                selected.push_back(vid);
            }
        }
    }

    vector<RecoverySampleVertex> sample;
    sample.reserve(selected.size());
    for (int32_t vid : selected) {
        if (vid < 0 || vid >= static_cast<int32_t>(global_vectex_vec.size())) {
            continue;
        }
        Vertex *v = global_vectex_vec[vid];
        if (v == nullptr) {
            continue;
        }
        RecoverySampleVertex entry;
        entry.vid = vid;
        collect_vertex_neighbors_for_debug(v, entry.neighbors);
        sample.push_back(std::move(entry));
    }
    return sample;
}

static bool validate_recovery_sample(const vector<RecoverySampleVertex> &sample)
{
    uint64_t checked_edges = 0;
    uint64_t mismatched_vertices = 0;
    uint64_t missing_vertices = 0;
    uint64_t degree_mismatches = 0;
    uint64_t edge_mismatches = 0;
    int printed = 0;

    for (const RecoverySampleVertex &expected : sample) {
        checked_edges += expected.neighbors.size();
        if (expected.vid < 0 ||
            expected.vid >= static_cast<int32_t>(global_vectex_vec.size()) ||
            global_vectex_vec[expected.vid] == nullptr) {
            missing_vertices++;
            mismatched_vertices++;
            continue;
        }

        vector<int32_t> actual;
        collect_vertex_neighbors_for_debug(global_vectex_vec[expected.vid], actual);
        bool ok = (actual == expected.neighbors);
        if (!ok) {
            mismatched_vertices++;
            if (actual.size() != expected.neighbors.size()) {
                degree_mismatches++;
            }
            vector<int32_t> expected_minus_actual;
            vector<int32_t> actual_minus_expected;
            set_difference(expected.neighbors.begin(), expected.neighbors.end(),
                           actual.begin(), actual.end(),
                           back_inserter(expected_minus_actual));
            set_difference(actual.begin(), actual.end(),
                           expected.neighbors.begin(), expected.neighbors.end(),
                           back_inserter(actual_minus_expected));
            edge_mismatches += expected_minus_actual.size() + actual_minus_expected.size();
            if (printed < 5) {
                cout << "Recovery sample mismatch vid=" << expected.vid
                     << " expected_degree=" << expected.neighbors.size()
                     << " actual_degree=" << actual.size()
                     << " missing_edges=" << expected_minus_actual.size()
                     << " extra_edges=" << actual_minus_expected.size()
                     << endl;
                printed++;
            }
        }
    }

    cout << "Recovery sample validation: checked_vertices=" << sample.size()
         << " checked_edges=" << checked_edges
         << " mismatched_vertices=" << mismatched_vertices
         << " missing_vertices=" << missing_vertices
         << " degree_mismatches=" << degree_mismatches
         << " edge_mismatches=" << edge_mismatches
         << " result=" << (mismatched_vertices == 0 ? "PASS" : "FAIL")
         << endl;
    return mismatched_vertices == 0;
}

static void clear_dram_graph_metadata_for_recovery(void)
{
    for (Vertex *v : global_vectex_vec) {
        if (v == nullptr) {
            continue;
        }
        for (pm_blk &blk : v->blk_list) {
            delete[] blk.buf;
            blk.buf = nullptr;
        }
        free(v->pma);
        v->pma = nullptr;
        delete v;
    }
    global_vectex_vec.clear();
    global_edge_vec.clear();
    global_edge_num = 0;
    global_max_vid = 0;
    global_cache_blk = 0;
    global_pending_edge_count = 0;
    global_committed_op_id = 0;
    global_last_applied_op_id = 0;
    global_edge_log_replay_active = false;

    for (ZoneStatus *zone : global_zone_vec) {
        if (zone == nullptr) {
            continue;
        }
        zone->val_blocks = 0;
        for (blk_info &info : zone->blk_info_vec) {
            info = blk_info();
        }
    }
}

void recover_graph_from_zns_records(void)
{
    if (!global_zns_record_enabled) {
        std::cerr << "--recover-zns requires ZNS record mode." << std::endl;
        exit(EXIT_FAILURE);
    }

    unordered_map<uint64_t, ZNSBlockHeader> committed_epochs;
    unordered_map<uint64_t, RecoveredEpochStats> observed_epoch_stats;
    unordered_set<uint64_t> valid_epochs;
    unordered_map<uint64_t, RecoveredRecord> latest_records;
    unordered_set<uint64_t> recovered_block_keys;
    vector<EdgeLogRecord> redo_records;
    int32_t max_vid = -1;
    uint64_t max_seq = 0;
    uint64_t max_epoch = 0;
    uint64_t max_op_id = 0;
    uint64_t meta_slots_scanned = 0;
    uint64_t data_slots_scanned = 0;
    uint64_t redo_slots_scanned = 0;
    uint64_t commit_records_valid = 0;
    uint64_t data_records_valid = 0;
    uint64_t data_records_committed_epoch = 0;
    uint64_t redo_records_valid = 0;
    bool redo_file_used = false;
    uint64_t rebuilt_vertices = 0;
    uint64_t rebuilt_blocks = 0;
    uint64_t missing_blocks_materialized = 0;
    uint64_t rebuilt_edges = 0;
    double start_t = get_current_time();

    if (global_meta_zone < 0) {
        std::cerr << "ZNS-native recovery requires a metadata zone." << std::endl;
        exit(EXIT_FAILURE);
    }

    unsigned long long meta_wp =
        static_cast<unsigned long long>(global_meta_zone) *
        static_cast<unsigned long long>(zns.zone_size);
    unsigned long long meta_scan_limit =
        global_meta_wp >= meta_wp ? global_meta_wp : (meta_wp + zns.zone_size);
    while (meta_wp < meta_scan_limit) {
        meta_slots_scanned++;
        ZNSBlockHeader header;
        if (readZNSRecordChecked(zns.dev_fd, header, nullptr, 0, meta_wp) &&
            header.record_type == ZNS_RECORD_COMMIT) {
            commit_records_valid++;
            auto it = committed_epochs.find(header.epoch);
            if (it == committed_epochs.end() || header.seq > it->second.seq) {
                committed_epochs[header.epoch] = header;
            }
            max_seq = max(max_seq, header.seq);
            max_epoch = max(max_epoch, header.epoch);
            max_op_id = max(max_op_id, header.committed_op_id);
            global_meta_wp = meta_wp + zns_alloc_unit_size();
        }
        meta_wp += zns_alloc_unit_size();
    }

    for (int32_t zidx = 0; zidx < static_cast<int32_t>(zns.nr_zones); zidx++) {
        if (zidx == global_meta_zone || zidx == global_redo_zone) {
            continue;
        }
        unsigned long long wp = static_cast<unsigned long long>(zidx) * zns.zone_size;
        unsigned long long scan_limit = recovery_zone_scan_limit(zidx);
        while (wp < scan_limit) {
            data_slots_scanned++;
            ZNSBlockHeader header;
            int32_t payload[PM_BLK_SIZE];
            if (!readZNSRecordChecked(zns.dev_fd, header, payload, PM_BLK_SIZE, wp) ||
                header.record_type != ZNS_RECORD_DATA) {
                wp += zns_alloc_unit_size();
                continue;
            }

            data_records_valid++;
            max_seq = max(max_seq, header.seq);
            max_epoch = max(max_epoch, header.epoch);
            if (committed_epochs.find(header.epoch) != committed_epochs.end()) {
                data_records_committed_epoch++;
                RecoveredEpochStats &stats = observed_epoch_stats[header.epoch];
                if (stats.data_record_count == 0) {
                    stats.first_seq = header.seq;
                    stats.last_seq = header.seq;
                } else {
                    stats.first_seq = min(stats.first_seq, header.seq);
                    stats.last_seq = max(stats.last_seq, header.seq);
                }
                stats.data_record_count++;
                stats.aggregate_checksum =
                    combine_epoch_checksum(stats.aggregate_checksum,
                                           header.seq,
                                           header.checksum);
            }
            wp += zns_alloc_unit_size();
        }
    }

    for (const auto &entry : committed_epochs) {
        auto stats_it = observed_epoch_stats.find(entry.first);
        if (stats_it != observed_epoch_stats.end() &&
            epoch_stats_match(entry.second, stats_it->second)) {
            valid_epochs.insert(entry.first);
            global_committed_op_id =
                max(global_committed_op_id, entry.second.committed_op_id);
        }
    }

    redo_file_used = read_redo_log_file_records(redo_records,
                                                redo_records_valid,
                                                max_op_id,
                                                max_vid);

    if (!redo_file_used && global_redo_zone >= 0) {
        unsigned long long redo_wp =
            static_cast<unsigned long long>(global_redo_zone) *
            static_cast<unsigned long long>(zns.zone_size);
        unsigned long long redo_scan_limit =
            global_redo_wp >= redo_wp ? global_redo_wp : (redo_wp + zns.zone_size);
        while (redo_wp < redo_scan_limit) {
            redo_slots_scanned++;
            vector<EdgeLogRecord> page_records;
            if (!readZNSRedoLogPageChecked(zns.dev_fd, page_records, redo_wp)) {
                redo_wp += zns_alloc_unit_size();
                continue;
            }
            redo_records_valid += page_records.size();
            for (const EdgeLogRecord &record : page_records) {
                redo_records.push_back(record);
                max_op_id = max(max_op_id, record.op_id);
                max_vid = max(max_vid, max(record.src, record.dst));
            }
            global_redo_wp = redo_wp + zns_alloc_unit_size();
            redo_wp += zns_alloc_unit_size();
        }
        sort(redo_records.begin(), redo_records.end(),
             [](const EdgeLogRecord &a, const EdgeLogRecord &b) {
                 return a.op_id < b.op_id;
             });
    }

    for (int32_t zidx = 0; zidx < static_cast<int32_t>(zns.nr_zones); zidx++) {
        if (zidx == global_meta_zone || zidx == global_redo_zone) {
            continue;
        }
        unsigned long long wp = static_cast<unsigned long long>(zidx) * zns.zone_size;
        unsigned long long scan_limit = recovery_zone_scan_limit(zidx);
        while (wp < scan_limit) {
            ZNSBlockHeader header;
            int32_t payload[PM_BLK_SIZE];
            if (!readZNSRecordChecked(zns.dev_fd, header, payload, PM_BLK_SIZE, wp) ||
                header.record_type != ZNS_RECORD_DATA ||
                valid_epochs.find(header.epoch) == valid_epochs.end()) {
                wp += zns_alloc_unit_size();
                continue;
            }

            uint64_t key = recovered_record_key(header.src_vid, header.pma_blk_idx);
            auto it = latest_records.find(key);
            if (it == latest_records.end() ||
                header.epoch > it->second.header.epoch ||
                (header.epoch == it->second.header.epoch &&
                 header.seq > it->second.header.seq)) {
                latest_records[key] = {header, wp};
            }
            max_vid = max(max_vid, header.src_vid);
            wp += zns_alloc_unit_size();
        }
    }

    global_zns_record_seq = max_seq + 1;
    global_current_epoch = max_epoch + 1;
    global_next_op_id = max_op_id + 1;
    global_last_applied_op_id = global_committed_op_id;
    global_max_vid = max_vid;
    global_vectex_vec.clear();
    if (max_vid >= 0) {
        global_vectex_vec.resize(max_vid + 1);
    }
    for (int32_t vid = 0; vid <= max_vid; vid++) {
        Vertex *v = new Vertex(vid);
        v->pma = pma_creat();
        global_vectex_vec[vid] = v;
        rebuilt_vertices++;
    }

    for (const auto &entry : latest_records) {
        const RecoveredRecord &record = entry.second;
        int32_t vid = record.header.src_vid;
        int32_t blk_idx = record.header.pma_blk_idx;
        if (vid < 0 || vid >= static_cast<int32_t>(global_vectex_vec.size()) || blk_idx < 0) {
            continue;
        }
        Vertex *v = global_vectex_vec[vid];
        if (static_cast<int32_t>(v->blk_list.size()) <= blk_idx) {
            v->blk_list.resize(blk_idx + 1);
        }
        v->blk_list[blk_idx].buf = nullptr;
        v->blk_list[blk_idx].global_wp = record.wp;
        recovered_block_keys.insert(recovered_record_key(vid, blk_idx));
        ZNSBlockHeader checked_header;
        int32_t payload[PM_BLK_SIZE];
        if (readZNSRecordChecked(zns.dev_fd, checked_header, payload, PM_BLK_SIZE, record.wp) &&
            checked_header.record_type == ZNS_RECORD_DATA) {
            refresh_pm_blk_meta(v->blk_list[blk_idx], payload);
        } else {
            reset_pm_blk_meta(v->blk_list[blk_idx]);
        }
        v->degree += record.header.valid_count;
        rebuilt_edges += record.header.valid_count;
        v->pma->ele_num = max(v->pma->ele_num, static_cast<int32_t>(record.header.pma_ele));
        v->pma->arr_cap = max(v->pma->arr_cap,
                              max(static_cast<int32_t>(record.header.pma_cap),
                                  (blk_idx + 1) * PM_BLK_SIZE));
        mark_zone(vid, blk_idx, record.wp);
        rebuilt_blocks++;
    }

    for (Vertex *v : global_vectex_vec) {
        refresh_recovered_pma_shape(v);
        if (v != nullptr && v->pma != nullptr) {
            ensure_pma_blocks_for_capacity(v, v->pma);
        }
    }
    missing_blocks_materialized =
        materialize_missing_recovered_blocks(recovered_block_keys);

    // A COMMIT high-water mark is only sufficient when the checkpoint also
    // covers cache-only vertices. Replaying the whole idempotent edge log is
    // conservative and restores DRAM-only small-degree updates after crashes.
    uint64_t redo_records_replayed = replay_redo_logs_after(0, redo_records);

    double end_t = get_current_time();
    cout << "ZNS-native recovery records: " << latest_records.size()
         << " committed_epochs: " << committed_epochs.size()
         << " valid_epochs: " << valid_epochs.size()
         << " committed_op_id: " << global_committed_op_id
         << " redo_records: " << redo_records.size()
         << " max_vid: " << global_max_vid
         << " time: " << (end_t - start_t) << endl;
    cout << "ZNS recovery stats: "
         << "meta_slots_scanned=" << meta_slots_scanned
         << " commit_records_valid=" << commit_records_valid
         << " committed_epochs=" << committed_epochs.size()
         << " valid_epochs=" << valid_epochs.size()
         << " data_slots_scanned=" << data_slots_scanned
         << " data_records_valid=" << data_records_valid
         << " data_records_committed_epoch=" << data_records_committed_epoch
         << " latest_records=" << latest_records.size()
         << " redo_slots_scanned=" << redo_slots_scanned
         << " redo_file_used=" << (redo_file_used ? 1 : 0)
         << " redo_records_valid=" << redo_records_valid
         << " redo_records_replayed=" << redo_records_replayed
         << " rebuilt_vertices=" << rebuilt_vertices
         << " rebuilt_blocks=" << rebuilt_blocks
         << " missing_blocks_materialized=" << missing_blocks_materialized
         << " rebuilt_edges_before_redo=" << rebuilt_edges
         << " committed_op_id=" << global_committed_op_id
         << " max_seq=" << max_seq
         << " max_epoch=" << max_epoch
         << " max_op_id=" << max_op_id
         << " recovery_time_sec=" << (end_t - start_t)
         << endl;
}

void run_recovery_experiment_from_zns(void)
{
    if (!global_zns_record_enabled) {
        std::cerr << "Recovery experiment requires ZNS record mode." << std::endl;
        exit(EXIT_FAILURE);
    }

    flushZNSRedoLog();
    vector<RecoverySampleVertex> validation_sample =
        build_recovery_validation_sample(128);
    double start_t = get_current_time();
    clear_dram_graph_metadata_for_recovery();
    double after_clear_t = get_current_time();
    recover_graph_from_zns_records();
    bool validation_ok = validate_recovery_sample(validation_sample);
    double end_t = get_current_time();
    cout << "Recovery experiment time: total=" << (end_t - start_t)
         << " clear_dram=" << (after_clear_t - start_t)
         << " zns_and_redo=" << (end_t - after_clear_t)
         << " graph_maintenance_total=" << last_graph_maintenance_total_time
         << " graph_maintenance_compute=" << last_graph_maintenance_compute_time
         << " recovery_vs_load="
         << (last_graph_maintenance_total_time > 0.0 ?
             (end_t - after_clear_t) / last_graph_maintenance_total_time : 0.0)
         << " validation=" << (validation_ok ? "PASS" : "FAIL")
         << endl;
}

// 恢复流程：先根据 zone 元数据恢复物理地址与度数，再尝试重建部分 PMA 状态。
void graph_recovery(void)
{
    for (int vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
        global_vectex_vec[vidx]->degree = 0;
    }
    double start_t = get_current_time();
    int32_t *tmp_blk_ptr = nullptr;
    unsigned long long wp = 0;
    for (int zidx = 0; zidx < zns.nr_zones; zidx++) {
        for (int blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
            if (global_zone_vec[zidx]->blk_info_vec[blkidx].status) {
                Vertex *v = global_vectex_vec[global_zone_vec[zidx]->blk_info_vec[blkidx].src];
                // if (v->blk_list[global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex].global_wp != wp) {
                //     cout << "error test: " << global_zone_vec[zidx]->blk_info_vec[blkidx].status << endl;
                // }
                // cout << "ori wp: " << v->blk_list[global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex].global_wp << " cur wp: " << wp << endl;
                int32_t logical_blk_idx = global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex;
                v->blk_list[logical_blk_idx].global_wp = wp;
                readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, wp);
                refresh_pm_blk_meta(v->blk_list[logical_blk_idx], tmp_blk_ptr);
                for (int j=0; j < PM_BLK_SIZE; j++) {
                    if (tmp_blk_ptr[j] >= 0){
                        v->degree++;
                    }
                }
                delete[] tmp_blk_ptr;
                tmp_blk_ptr = nullptr;
            }
            wp += zns_alloc_unit_size();
        }
    }
    cout << "test --0" << endl;
    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for
    for (int64_t edgeIdx = (global_edge_num-1); edgeIdx >= 0; edgeIdx--) {
        new_edge e = global_edge_vec[edgeIdx];  // 从边数组中取出当前边
        int32_t src = e.src;
        int32_t des = e.des;
        Vertex* v = global_vectex_vec[src];
        unique_lock<std::shared_timed_mutex> lock(v->mutex); // 保证同一顶点同一时刻只处理一条边
        int32_t index;
        if (v->blk_list.empty()) { // 小度顶点当前不做复杂恢复
            // v->cache_buf[v->degree] = des;
            // v->degree++;
            continue;
        } 

        // index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
        if (!v->recovery_flag) {
            index = binary_pma_find_recovery(v->pma, v, des, 0, v->pma->arr_cap-1);
            v->recovery_flag = true;
            if (index == -1) {
                continue;
            } else {
                cout << "need to inert" << endl;
                pma_insert(v->pma, v, des, index);
            }
        }
        
    }

    double end_t = get_current_time();
	cout << "Graph recovery time: " << end_t - start_t << endl;
}

#endif
