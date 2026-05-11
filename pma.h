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

#define BLK_SCALE 1
// #define PM_BLK_SIZE (32 * BLK_SCALE)  // 每块的元素数可在实验中调整

PMA* pma_creat (void);
static void compute_capacity (PMA* pma);
int32_t binary_pma_find(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to);
int32_t pma_insert(PMA* pma, Vertex *src, int32_t des, int32_t index);
void rebalance(PMA* pma, Vertex *src, int32_t index);
void spread (PMA* pma, Vertex *src, int32_t from, int32_t to, int32_t occupancy);
void resize (PMA* pma, Vertex *src);
static inline int32_t pma_value_at(Vertex *src, int32_t index);
static int32_t count_window_occupancy(Vertex *src, int32_t from, int32_t to);

int load_graph(string path);
void graph_init(int32_t max_vid);
void graph_maintenance_baseline(void);
bool cache_insert_unique(Vertex *v, int32_t des);
bool cache_delete(Vertex *v, int32_t des);
void pma_compact(PMA *pma, Vertex *src);
bool pma_delete(PMA *pma, Vertex *src, int32_t des);
bool graph_insert_directed_edge(int32_t src, int32_t des);
bool graph_delete_directed_edge(int32_t src, int32_t des);
bool graph_delete_undirected_edge(int32_t src, int32_t des);
VertexClass classify_vertex(Vertex *v);
const char *vertex_class_name(VertexClass vclass);
bool requestClassifiedBlock(Vertex *v, unsigned long long &wp);
bool garbage_collection_by_mode(VertexClass requested_class);
void print_vertex(int32_t vid, int64_t &total_edges);
void print_graph(void);
void print_vertex_pma_blocks(int32_t vid);
void print_all_vertex_pma_blocks(void);

// void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp);
void dumpToZone(void);
void dumpToZoneParallel(void);
void dumpToZoneZoneSerial(void);
void persist_cache_vertices_as_records(void);
void recover_graph_from_zns_records(void);

void graph_init_new(int32_t max_vid);
void graph_maintenance_new(void);
void graph_maintenance_parallel(void);
void graph_maintenance_parallel2(void);

void graph_recovery(void);
int32_t binary_pma_find_recovery(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to);


// 创建 PMA 元数据，初始容量很小，后续依靠 resize 扩容。
PMA* pma_creat (void)
{
    PMA* pma = (PMA*)malloc(sizeof(PMA));
    pma->ele_num = 0;
    pma->arr_cap = 4;
    pma->seg_len = 2;
    pma->seg_num = pma->arr_cap / pma->seg_len;
    pma->tree_h = floor_lg (pma->seg_num) + 1;
    pma->delta_t = (t_0 - t_h) / pma->tree_h;  // 上界密度按层级线性变化
    pma->delta_p = (p_h - p_0) / pma->tree_h;
    // cout << "tree_h: " << pma->tree_h << " seg_num: " << pma->seg_num << endl;
    // pma->array = (element_type *)memkind_malloc(kind, sizeof(element_type) * pma->arr_cap);  // on NVM
    // pma->row_offset = (vertexOff *)malloc(sizeof(vertexOff) * pma->vertex_cap); // on DRAM
    return (pma);
}

// 根据当前元素数重新计算 PMA 容量和段长。
static void compute_capacity (PMA* pma) 
{
    /* 根据已有元素数重新计算 seg_len 和 arr_cap */
    pma->seg_len = ceil_lg (pma->ele_num);  /* Ideal segment size. */
    pma->seg_num = ceil_div (pma->ele_num, pma->seg_len);  /* Ideal number of segments. */
    /* The number of segments has to be a power of 2. */
    pma->seg_num = hyperceil (pma->seg_num);
    /* Update the segment size accordingly. */
    pma->seg_len = ceil_div (pma->ele_num, pma->seg_num);
    pma->arr_cap = pma->seg_len * pma->seg_num;
    /* Scale up as much as possible. */
    pma->arr_cap = max_sparseness * pma->arr_cap;
    pma->seg_len = max_sparseness * pma->seg_len;

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

static int32_t count_window_occupancy(Vertex *src, int32_t from, int32_t to) {
    int32_t occupancy = 0;
    for (int32_t idx = from; idx < to; idx++) {
        int32_t blk_idx = idx / PM_BLK_SIZE;
        int32_t bit_idx = idx % PM_BLK_SIZE;
        if (src->blk_list[blk_idx].buf != nullptr) {
            occupancy += src->blk_list[blk_idx].buf[bit_idx] >= 0 ? 1 : 0;
        } else {
            int32_t *tmp_blk_ptr = nullptr;
            readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp);
            occupancy += tmp_blk_ptr[bit_idx] >= 0 ? 1 : 0;
            delete[] tmp_blk_ptr;
        }
    }
    return occupancy;
}


// 在 PMA 中查找 des 的前驱位置。
// 若访问到未缓存的块，会按需从设备换入。
int32_t binary_pma_find(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to) 
{
    if (to < 0) {
        return -1;
    }

    int32_t index;
    int32_t blk_idx, bit_idx;
    // int32_t cur_blk;
    // int32_t *tmp_blk_ptr = new int32_t[PM_BLK_SIZE];  // 临时块缓冲
    vector<int32_t> release_blk_vec;

    while (from <= to) {
        int32_t mid = from + (to - from) / 2;
        int32_t i = mid;
        /* Start scanning left until we find a non-empty slot or we reach past the
        * beginning of the subarray. */
        blk_idx = i / PM_BLK_SIZE;
        bit_idx = i % PM_BLK_SIZE;

        if (src->blk_list[blk_idx].buf == nullptr) {
            readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从设备读入一个块
            release_blk_vec.push_back(blk_idx);
        }
        
        // cur_blk = blk_idx;
        // while (i >= from && (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)){
        while (i >= from && (src->blk_list[blk_idx].buf[bit_idx] < 0)){
            i--;
            if (i < from) {
                break;
            }
            blk_idx = i / PM_BLK_SIZE;
            bit_idx = i % PM_BLK_SIZE;
            // 跨块后按需把对应块换入内存
            if (src->blk_list[blk_idx].buf == nullptr) {
                readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从设备读入一个块
                release_blk_vec.push_back(blk_idx);
            }
        }

        if (i < from) {  /* Everything between from and mid (inclusive) is empty. */
            from = mid + 1;
        } else {
            // if (src->blk_list[blk_idx].blk_ptr[bit_idx] == des) {
            if (src->blk_list[blk_idx].buf[bit_idx] == des) {
                index = i;
                for (int jj = 0; jj < release_blk_vec.size(); jj++) {
                    delete[] src->blk_list[release_blk_vec[jj]].buf;
                    src->blk_list[release_blk_vec[jj]].buf = nullptr; // 查找结束后释放临时换入的块
                }
                // cout << "error: insert the same neighbor. src: " << src << " des: " << des << endl;
                return index;
            }
            else if (src->blk_list[blk_idx].buf[bit_idx] < des) {
                from = mid + 1;
            } else {  /* pma->array [i].key > key */
                to = i - 1;  // 保证最终得到的是小于 des 的最大位置
            }
        }
    }
    /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
    index = to;
    if (index < 0) {
        for (int jj = 0; jj < release_blk_vec.size(); jj++) {
            delete[] src->blk_list[release_blk_vec[jj]].buf;
            src->blk_list[release_blk_vec[jj]].buf = nullptr;
        }
        return -1;
    }

    blk_idx = index / PM_BLK_SIZE;
    bit_idx = index % PM_BLK_SIZE;

    if (src->blk_list[blk_idx].buf == nullptr) {
        readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从设备读入一个块
        release_blk_vec.push_back(blk_idx);
    }

    // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp);

    // cur_blk = blk_idx;
    while (index >= 0 && (src->blk_list[blk_idx].buf[bit_idx] < 0)) {
        index--;  // 跳过空洞，确保找到真实前驱
        if (index < 0) {
            break;
        }
        blk_idx = index / PM_BLK_SIZE;
        bit_idx = index % PM_BLK_SIZE;
        // 跨块后按需把对应块换入内存

        if (src->blk_list[blk_idx].buf == nullptr) {
            readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从设备读入一个块
            release_blk_vec.push_back(blk_idx);
        }
    }

    for (int jj = 0; jj < release_blk_vec.size(); jj++) {
        delete[] src->blk_list[release_blk_vec[jj]].buf;
        src->blk_list[release_blk_vec[jj]].buf = nullptr; // 查找结束后释放临时换入的块
    }
    return index;
}

// 恢复阶段使用的查找版本；若 des 已存在则返回 -1。
int32_t binary_pma_find_recovery(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to) 
{
    if (to < 0) {
        return -1;
    }

    int32_t index;
    int32_t blk_idx, bit_idx;
    // int32_t cur_blk;
    // int32_t *tmp_blk_ptr = new int32_t[PM_BLK_SIZE];  // 临时块缓冲
    vector<int32_t> release_blk_vec;

    while (from <= to) {
        int32_t mid = from + (to - from) / 2;
        int32_t i = mid;
        /* Start scanning left until we find a non-empty slot or we reach past the
        * beginning of the subarray. */
        blk_idx = i / PM_BLK_SIZE;
        bit_idx = i % PM_BLK_SIZE;

        if (src->blk_list[blk_idx].buf == nullptr) {
            readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从设备读入一个块
            release_blk_vec.push_back(blk_idx);
        }
        
        // cur_blk = blk_idx;
        // while (i >= from && (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)){
        while (i >= from && (src->blk_list[blk_idx].buf[bit_idx] < 0)){
            i--;
            if (i < from) {
                break;
            }
            blk_idx = i / PM_BLK_SIZE;
            bit_idx = i % PM_BLK_SIZE;
            // 跨块后按需把对应块换入内存
            if (src->blk_list[blk_idx].buf == nullptr) {
                readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从设备读入一个块
                release_blk_vec.push_back(blk_idx);
            }
        }

        if (i < from) {  /* Everything between from and mid (inclusive) is empty. */
            from = mid + 1;
        } else {
            // if (src->blk_list[blk_idx].blk_ptr[bit_idx] == des) {
            if (src->blk_list[blk_idx].buf[bit_idx] == des) {
                index = i;
                for (int jj = 0; jj < release_blk_vec.size(); jj++) {
                    delete[] src->blk_list[release_blk_vec[jj]].buf;
                    src->blk_list[release_blk_vec[jj]].buf = nullptr; // 查找结束后释放临时换入的块
                }
                // cout << "error: insert the same neighbor. src: " << src << " des: " << des << endl;
                return -1;
            }
            else if (src->blk_list[blk_idx].buf[bit_idx] < des) {
                from = mid + 1;
            } else {  /* pma->array [i].key > key */
                to = i - 1;  // 保证最终得到的是小于 des 的最大位置
            }
        }
    }
    /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
    index = to;
    if (index < 0) {
        for (int jj = 0; jj < release_blk_vec.size(); jj++) {
            delete[] src->blk_list[release_blk_vec[jj]].buf;
            src->blk_list[release_blk_vec[jj]].buf = nullptr;
        }
        return -1;
    }

    blk_idx = index / PM_BLK_SIZE;
    bit_idx = index % PM_BLK_SIZE;

    if (src->blk_list[blk_idx].buf == nullptr) {
        readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从设备读入一个块
        release_blk_vec.push_back(blk_idx);
    }

    // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp);

    // cur_blk = blk_idx;
    while (index >= 0 && (src->blk_list[blk_idx].buf[bit_idx] < 0)) {
        index--;  // 跳过空洞，确保找到真实前驱
        if (index < 0) {
            break;
        }
        blk_idx = index / PM_BLK_SIZE;
        bit_idx = index % PM_BLK_SIZE;
        // 跨块后按需把对应块换入内存

        if (src->blk_list[blk_idx].buf == nullptr) {
            readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从设备读入一个块
            release_blk_vec.push_back(blk_idx);
        }
    }

    for (int jj = 0; jj < release_blk_vec.size(); jj++) {
        delete[] src->blk_list[release_blk_vec[jj]].buf;
        src->blk_list[release_blk_vec[jj]].buf = nullptr; // 查找结束后释放临时换入的块
    }
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

    return inserted_index;
}

bool pma_delete(PMA *pma, Vertex *src, int32_t des)
{
    if (pma == nullptr || src->blk_list.empty() || pma->ele_num == 0) {
        return false;
    }

    int32_t index = binary_pma_find(pma, src, des, 0, pma->arr_cap - 1);
    if (index < 0 || pma_value_at(src, index) != des) {
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
}


// 逐层扩大窗口，查找满足密度阈值的最小窗口。
void rebalance(PMA* pma, Vertex *src, int32_t index) 
{
    // 通过统计窗口内有效元素密度，决定执行 spread 还是 resize
    int32_t window_start, window_end;
    int32_t height = 0;
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
      /* Found a window within threshold. */
    if (height == 1)
        return;
    
    if (density >= p_height && density < t_height) {
        // cout << "spread process ---" << endl;
        spread (pma, src, window_start, window_end, occupancy);
    } else {
        // print_graph();
        // cout << "resize process ---" << endl;
        resize (pma, src);
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

    int32_t cur_blk = -1;
    int32_t *tmp_blk_ptr = new int32_t[PM_BLK_SIZE];

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

        // tmpResBitmapIdx = 0;
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
    }
    free(tmp_win_arr);
    delete[] tmp_blk_ptr;
    // cout << "spread done" << endl;
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
        pm_blk new_blk;
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
        pm_blk new_blk;
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
            global_cache_blk = 0;
        }
    }
    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);
    
	double end_t = get_current_time();
    // print_graph();
    // print_pma(global_vectex_vec[1]);
	cout << "Load graph on NVM baseline time: " << end_t - start_t << endl;
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

bool graph_insert_directed_edge(int32_t src, int32_t des)
{
    if (src < 0 || src >= static_cast<int32_t>(global_vectex_vec.size())) {
        return false;
    }

    Vertex *v = global_vectex_vec[src];
    unique_lock<std::shared_timed_mutex> lock(v->mutex);
    if (v->blk_list.empty()) {
        for (int32_t i = 0; i < v->degree; i++) {
            if (v->cache_buf[i] == des) {
                return false;
            }
        }
        if (v->degree < CACHE_BLK_SIZE) {
            bool inserted = cache_insert_unique(v, des);
            if (inserted) {
                v->update_count++;
                classify_vertex(v);
            }
            return inserted;
        }

        pm_blk new_blk;
        new_blk.buf = new int32_t[PM_BLK_SIZE];
        fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
        v->blk_list.push_back(new_blk);
        global_cache_blk++;
        int32_t cached_degree = v->degree;
        v->degree = 0;
        for (int32_t i = 0; i < cached_degree; i++) {
            int32_t cached_des = v->cache_buf[i];
            int32_t index = binary_pma_find(v->pma, v, cached_des, 0, v->pma->arr_cap - 1);
            pma_insert(v->pma, v, cached_des, index);
        }
        fill(v->cache_buf, v->cache_buf + CACHE_BLK_SIZE, -1);
    }

    int32_t index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap - 1);
    int32_t old_degree = v->degree;
    pma_insert(v->pma, v, des, index);
    bool inserted = v->degree != old_degree;
    if (inserted) {
        v->update_count++;
        classify_vertex(v);
    }
    return inserted;
}

bool graph_delete_directed_edge(int32_t src, int32_t des)
{
    if (src < 0 || src >= static_cast<int32_t>(global_vectex_vec.size())) {
        return false;
    }
    Vertex *v = global_vectex_vec[src];
    unique_lock<std::shared_timed_mutex> lock(v->mutex);
    bool deleted = false;
    if (v->blk_list.empty()) {
        deleted = cache_delete(v, des);
    } else {
        deleted = pma_delete(v->pma, v, des);
    }
    if (deleted) {
        v->update_count++;
        classify_vertex(v);
    }
    return deleted;
}

bool graph_delete_undirected_edge(int32_t src, int32_t des)
{
    bool forward_deleted = graph_delete_directed_edge(src, des);
    bool backward_deleted = graph_delete_directed_edge(des, src);
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

static const vector<int32_t> &zoneIdsForClass(VertexClass vclass)
{
    if (vclass == VERTEX_HOT) {
        return hot_zone_ids;
    }
    if (vclass == VERTEX_HUB) {
        return hub_zone_ids;
    }
    return cold_zone_ids;
}

static bool requestRelocationBlock(Vertex *v, int32_t excluded_zone, unsigned long long &wp)
{
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
        int32_t invalid = zns.nr_blocks_perzone - valid;
        if (invalid > best_invalid) {
            best_invalid = invalid;
            best_zone = zone_idx;
        }
    }
    return best_zone;
}

static void resetZoneMetadataAndQueue(int32_t zone_idx)
{
    if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_blk_queues.size())) {
        return;
    }
    zbd_reset_zones(zns.dev_fd,
                    static_cast<unsigned long long>(zone_idx) * zns.zone_size,
                    zns.zone_size);

    ZoneStatus *zone = global_zone_vec[zone_idx];
    zone->val_blocks = 0;
    for (int32_t blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
        zone->blk_info_vec[blkidx] = blk_info();
    }

    lock_guard<std::mutex> lock(global_blk_queues[zone_idx]->mtx_);
    while (!global_blk_queues[zone_idx]->memoryWPs.empty()) {
        global_blk_queues[zone_idx]->memoryWPs.pop();
    }
    unsigned long long wp = static_cast<unsigned long long>(zone_idx) * zns.zone_size;
    for (int32_t blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
        global_blk_queues[zone_idx]->memoryWPs.push(wp);
        wp += zns_alloc_unit_size();
    }
}

static bool garbageCollectVictimZone(int32_t victim_zone)
{
    if (victim_zone < 0 || victim_zone >= static_cast<int32_t>(global_zone_vec.size())) {
        return false;
    }

    cout << "GC victim zone: " << victim_zone << endl;
    int32_t moved = 0;
    int32_t *array = nullptr;
    unsigned long long rp = static_cast<unsigned long long>(victim_zone) * zns.zone_size;
    for (int32_t blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
        blk_info info = global_zone_vec[victim_zone]->blk_info_vec[blkidx];
        if (info.status) {
            Vertex *v = global_vectex_vec[info.src];
            unsigned long long new_wp = 0;
            if (!requestRelocationBlock(v, victim_zone, new_wp)) {
                delete[] array;
                return false;
            }

            readBlockFromZone(zns.dev_fd, array, PM_BLK_SIZE, rp);
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
        }
        rp += zns_alloc_unit_size();
    }
    delete[] array;

    resetZoneMetadataAndQueue(victim_zone);
    cout << "GC moved blocks: " << moved << endl;
    return true;
}

bool garbage_collection_by_mode(VertexClass requested_class)
{
    if (!global_zone_class_enabled || global_blk_queues.empty()) {
        garbage_collection();
        return true;
    }

    int32_t victim = -1;
    vector<int32_t> all_class_zones;
    all_class_zones.insert(all_class_zones.end(), hot_zone_ids.begin(), hot_zone_ids.end());
    all_class_zones.insert(all_class_zones.end(), cold_zone_ids.begin(), cold_zone_ids.end());
    all_class_zones.insert(all_class_zones.end(), hub_zone_ids.begin(), hub_zone_ids.end());

    if (global_gc_mode == GC_GLOBAL) {
        victim = chooseVictimZone(all_class_zones);
        return victim >= 0 && garbageCollectVictimZone(victim);
    }

    if (global_gc_mode == GC_CLASS_LOCAL) {
        victim = chooseVictimZone(zoneIdsForClass(requested_class));
        if (victim >= 0 && garbageCollectVictimZone(victim)) {
            return true;
        }
    }

    if (global_gc_mode == GC_CLASS_AWARE || global_gc_mode == GC_CLASS_LOCAL) {
        victim = chooseVictimZone(all_class_zones);
        if (victim >= 0 && garbageCollectVictimZone(victim)) {
            return true;
        }
    }

    return false;
}

bool requestClassifiedBlock(Vertex *v, unsigned long long &wp)
{
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

    if (garbage_collection_by_mode(vclass)) {
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
            cout << "dump to zone" << endl;
            dumpToZone();
            global_cache_blk = 0;
        }
    }

    if (global_cache_blk > 0) {
        dumpToZone();
        global_cache_blk = 0;
    }
    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);
    
	double end_t = get_current_time();
    // print_graph();
    // print_pma(global_vectex_vec[1]);
	cout << "Load graph on NVM new time: " << end_t - start_t << endl;
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
            int32_t src = e.src;
            int32_t des = e.des;
            Vertex* v = global_vectex_vec[src];
            unique_lock<std::shared_timed_mutex> lock(v->mutex); // 保证同一顶点同一时刻只处理一条边
            int32_t index;
            // cout << "insert edge " << edgeIdx << ": (" << src << ", " << des << ")" << endl; 

            if (v->blk_list.empty() && v->degree < CACHE_BLK_SIZE) { // 小度顶点直接进入缓存
                // cout << "test-- cache insert" << endl;
                cache_insert_unique(v, des);
            } else if (v->blk_list.empty() && v->degree == CACHE_BLK_SIZE) { // 达到阈值后，把缓存中的边整体迁入 PMA
                // cout << "test-- boardline insert" << endl;
                pm_blk new_blk;
                new_blk.buf = new int32_t[PM_BLK_SIZE]; // 在 DRAM 中分配内存
                fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
                v->blk_list.push_back(new_blk);
                global_cache_blk++;
                v->degree = 0;
                for (int i = 0; i < CACHE_BLK_SIZE; i++) {
                    des = v->cache_buf[i];
                    index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
                    pma_insert(v->pma, v, des, index);
                }
                des = e.des;  // 再插入当前这条新边
                index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
                pma_insert(v->pma, v, des, index);
            } else {
                // cout << "test-- PMA insert" << endl;
                index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
                pma_insert(v->pma, v, des, index);

            }
        }

        if (global_cache_blk >= cache_blk_threshold) {
            cout << "dump to zone" << endl;
            dumpToZone();
            global_cache_blk = 0;
        }

        if (end_e == global_edge_num) {
            break;
        }
        start_e = end_e;
        end_e = min(start_e + global_batch_size, global_edge_num);
    }

	double end_t = get_current_time();
	cout << "Load graph on NVM new time: " << end_t - start_t << endl;
    // print_graph();
}

// 并行插入 + 并行刷盘版本，是当前 main 默认使用的构图流程。
void graph_maintenance_parallel2(void)
{
    /*
        并行插入 + 并行刷盘版本：
        当前默认将整个数据集作为一个 batch，在 batch 结束后统一并行刷盘。
    */
    graph_init_new2(global_max_vid);
    cout << "Load graph on NVM new begin!" << endl;
	double start_t = get_current_time();

    // batch 结束后检查缓存块数量，达到阈值则并行刷盘
    omp_set_num_threads(THREAD_NUM);
    int64_t start_e = 0;
    global_batch_size = global_edge_num;
    int64_t end_e = min(start_e + global_batch_size, global_edge_num);
    
    while (end_e <= global_edge_num) {
        // cout << "end_e: " << end_e << endl;
        #pragma omp parallel for
        for (int64_t edgeIdx = start_e; edgeIdx < end_e; edgeIdx++)
        {
            new_edge e = global_edge_vec[edgeIdx];  // 从边数组中取出当前边
            int32_t src = e.src;
            int32_t des = e.des;
            Vertex* v = global_vectex_vec[src];
            unique_lock<std::shared_timed_mutex> lock(v->mutex); // 保证同一顶点同一时刻只处理一条边
            int32_t index;
            // cout << "insert edge " << edgeIdx << ": (" << src << ", " << des << ")" << endl; 

            if (v->blk_list.empty() && v->degree < CACHE_BLK_SIZE) { // 小度顶点直接进入缓存
                // cout << "test-- cache insert" << endl;
                cache_insert_unique(v, des);
            } else if (v->blk_list.empty() && v->degree == CACHE_BLK_SIZE) { // 达到阈值后，把缓存中的边整体迁入 PMA
                // cout << "test-- boardline insert" << endl;
                pm_blk new_blk;
                new_blk.buf = new int32_t[PM_BLK_SIZE]; // 在 DRAM 中分配内存
                fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
                v->blk_list.push_back(new_blk);
                global_cache_blk++;
                v->degree = 0;
                for (int i = 0; i < CACHE_BLK_SIZE; i++) {
                    des = v->cache_buf[i];
                    index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
                    pma_insert(v->pma, v, des, index);
                }
                des = e.des;  // 再插入当前这条新边
                index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
                pma_insert(v->pma, v, des, index);
            } else {
                // cout << "test-- PMA insert" << endl;
                index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
                pma_insert(v->pma, v, des, index);

            }
        }

        if (global_cache_blk >= cache_blk_threshold) {
            // cout << "dump to zone" << endl;
            dumpToZoneParallel();
            global_cache_blk = 0;
        }

        // if ((end_e % (global_batch_size * 10)) == 0)
        // {
        //     garbage_collection_new();
        // }
        

        if (end_e == global_edge_num) {
            break;
        }
        start_e = end_e;
        end_e = min(start_e + global_batch_size, global_edge_num);
    }


	double end_t = get_current_time();
	cout << "Load graph on NVM new time: " << end_t - start_t << endl;
    // print_graph();
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

        cout << "  blk[" << blkidx << "] wp=";
        if (loaded_from_zone) {
            cout << v->blk_list[blkidx].global_wp;
        } else {
            cout << "unassigned";
        }
        cout << " source=" << (loaded_from_zone ? "zone" : "dram")
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
    double start_t = get_current_time();
    for (int32_t vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        // cout << "dump vid: " << v->id << " degree: " << v->degree << " blks: " << v->blk_list.size() << endl;
        if ((v->degree > 0) && (v->blk_list.size() > 0)) {
            for (int blkidx = 0; blkidx < v->blk_list.size(); blkidx++) {
                if (v->blk_list[blkidx].buf == nullptr) { // 该块当前不在 DRAM 中
                    continue;
                }
                // total_blks++;
                // cout << "dump vid: " << v->id << endl;

                pm_blk new_blk;
                if (global_zone_class_enabled && !global_blk_queues.empty()) {
                    if (!requestClassifiedBlock(v, new_blk.global_wp)) {
                        std::cerr << "No free classified blocks for dump." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                } else {
                    if (graph_blk_que.size() == 0) {
                        garbage_collection();
                    }
                    new_blk = graph_blk_que.front();
                    graph_blk_que.pop();
                }

                writeGraphBlockToZone(v, blkidx, new_blk.global_wp);
                delete[] v->blk_list[blkidx].buf;
                mark_zone(v->id, blkidx, new_blk.global_wp);  // 更新块映射，供 GC 使用
                // 更新顶点侧的物理地址元数据
                v->blk_list[blkidx] = new_blk;
                // v->blk_list[blkidx].buf = nullptr;
                // v->blk_list[blkidx].global_wp = new_blk.global_wp;

                global_cache_blk--;
            }
        }
    }
    double end_t = get_current_time();
    cout << "Dump to zone time: " << end_t - start_t << endl;
    // cout << "total_blks: " << total_blks << endl;
}

// 并行刷盘：每个线程从可写块队列中申请物理地址并刷回设备。
void dumpToZoneParallel(void)
{
    omp_set_num_threads(THREAD_NUM);
    double start_t = get_current_time();
    #pragma omp parallel for
    for (int32_t vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
        int threadId = omp_get_thread_num();
        Vertex *v = global_vectex_vec[vidx];
        // cout << "dump vid: " << v->id << " degree: " << v->degree << " blks: " << v->blk_list.size() << endl;
        if ((v->degree > 0) && (v->blk_list.size() > 0)) {
            for (int blkidx = 0; blkidx < v->blk_list.size(); blkidx++) {
                if (v->blk_list[blkidx].buf == nullptr) { // 该块当前不在 DRAM 中
                    continue;
                }

                // if (graph_blk_que.size() == 0) {
                //     garbage_collection();
                // }

                unsigned long long wp = 0;
                if (global_zone_class_enabled) {
                    if (!requestClassifiedBlock(v, wp)) {
                        std::cerr << "No free classified blocks for parallel dump." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                } else {
                    wp = requestMemoryBlocks(threadId);
                }
                // cout << "dump wp: " << wp << endl;
                
                writeGraphBlockToZone(v, blkidx, wp);
                // cout << " write to zone done." << endl;
                delete[] v->blk_list[blkidx].buf;
                mark_zone(v->id, blkidx, wp);  // 更新块映射，供 GC 和恢复使用
                // 更新顶点侧的物理地址元数据
                // v->blk_list[blkidx] = new_blk;
                v->blk_list[blkidx].buf = nullptr;
                v->blk_list[blkidx].global_wp = wp;

                global_cache_blk--;
            }
        }
    }
    double end_t = get_current_time();
    cout << "Dump to zone parallel time: " << end_t - start_t << endl;
    total_dump_time += (end_t - start_t);
    // cout << "total_blks: " << total_blks << endl;
}

void dumpToZoneZoneSerial(void)
{
    if (global_blk_queues.empty() || global_num_queues <= 0) {
        dumpToZone();
        return;
    }

    double start_t = get_current_time();
    int next_zone = 0;
    for (int32_t vidx = 0; vidx < static_cast<int32_t>(global_vectex_vec.size()); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
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

                writeGraphBlockToZone(v, blkidx, wp);
                delete[] v->blk_list[blkidx].buf;
                mark_zone(v->id, blkidx, wp);
                v->blk_list[blkidx].buf = nullptr;
                v->blk_list[blkidx].global_wp = wp;
                global_cache_blk--;
            }
        }
    }
    double end_t = get_current_time();
    cout << "Dump to zone-serial writer time: " << end_t - start_t << endl;
    total_dump_time += (end_t - start_t);
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
    garbage_collection_by_mode(VERTEX_COLD);
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

    for (int32_t vidx = 0; vidx < static_cast<int32_t>(global_vectex_vec.size()); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        if (v->degree == 0 || !v->blk_list.empty()) {
            continue;
        }

        pm_blk blk;
        blk.buf = new int32_t[PM_BLK_SIZE];
        fill(blk.buf, blk.buf + PM_BLK_SIZE, -1);
        for (int32_t i = 0; i < v->degree; i++) {
            blk.buf[i] = v->cache_buf[i];
        }
        v->blk_list.push_back(blk);
        v->pma->ele_num = v->degree;
        v->pma->arr_cap = max(v->pma->arr_cap, CACHE_BLK_SIZE);

        unsigned long long wp = 0;
        if (!request_any_record_block(v, wp)) {
            std::cerr << "No free block to persist cache vertex record." << std::endl;
            exit(EXIT_FAILURE);
        }
        writeGraphBlockToZone(v, 0, wp);
        delete[] v->blk_list[0].buf;
        v->blk_list[0].buf = nullptr;
        v->blk_list[0].global_wp = wp;
        mark_zone(v->id, 0, wp);
    }
}

struct RecoveredRecord {
    ZNSBlockHeader header;
    unsigned long long wp;
};

void recover_graph_from_zns_records(void)
{
    if (!global_zns_record_enabled) {
        std::cerr << "--recover-zns requires --zns-record." << std::endl;
        exit(EXIT_FAILURE);
    }

    unordered_map<uint64_t, RecoveredRecord> latest_records;
    int32_t max_vid = -1;
    uint64_t max_seq = 0;
    double start_t = get_current_time();

    for (int32_t zidx = 0; zidx < static_cast<int32_t>(zns.nr_zones); zidx++) {
        unsigned long long wp = static_cast<unsigned long long>(zidx) * zns.zone_size;
        for (int32_t blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
            ZNSBlockHeader header;
            if (!readZNSRecordHeader(zns.dev_fd, header, wp)) {
                wp += zns_alloc_unit_size();
                continue;
            }
            int32_t *payload = nullptr;
            readZNSRecordPayload(zns.dev_fd, payload, PM_BLK_SIZE, wp);
            delete[] payload;

            uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(header.src_vid)) << 32) |
                           static_cast<uint32_t>(header.pma_blk_idx);
            auto it = latest_records.find(key);
            if (it == latest_records.end() || header.seq > it->second.header.seq) {
                latest_records[key] = {header, wp};
            }
            max_vid = max(max_vid, header.src_vid);
            max_seq = max(max_seq, header.seq);
            wp += zns_alloc_unit_size();
        }
    }

    global_zns_record_seq = max_seq + 1;
    global_max_vid = max_vid;
    global_vectex_vec.clear();
    if (max_vid >= 0) {
        global_vectex_vec.resize(max_vid + 1);
    }
    for (int32_t vid = 0; vid <= max_vid; vid++) {
        Vertex *v = new Vertex(vid);
        v->pma = pma_creat();
        global_vectex_vec[vid] = v;
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
        v->degree += record.header.valid_count;
        v->pma->ele_num = max(v->pma->ele_num, static_cast<int32_t>(record.header.pma_ele));
        v->pma->arr_cap = max(v->pma->arr_cap, static_cast<int32_t>(record.header.pma_cap));
        mark_zone(vid, blk_idx, record.wp);
    }

    double end_t = get_current_time();
    cout << "ZNS-native recovery records: " << latest_records.size()
         << " max_vid: " << global_max_vid
         << " time: " << (end_t - start_t) << endl;
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
                v->blk_list[global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex].global_wp = wp;
                readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, wp);
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
