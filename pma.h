#ifndef __PMA_H_
#define __PMA_H_

#include "global.h"
// #include "utils.h"

using namespace std;

/* Height-based (as opposed to depth-based) thresholds. */
/* Upper density thresholds. */
static const double t_h = 0.75;  /* root. */
static const double t_0 = 1.00;  /* leaves. */
/* Lower density thresholds. */
static const double p_h = 0.50;  /* root. */
static const double p_0 = 0.25;  /* leaves. */
static const uint8_t max_sparseness =  2; // 1 / p_0;

#define BLK_SCALE 1
// #define PM_BLK_SIZE (32 * BLK_SCALE)  // 32个元素  4字节一个元素 默认256slots

PMA* pma_creat (void);
static void compute_capacity (PMA* pma);
int32_t binary_pma_find(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to);
int32_t pma_insert(PMA* pma, Vertex *src, int32_t des, int32_t index);
void rebalance(PMA* pma, Vertex *src, int32_t index);
void spread (PMA* pma, Vertex *src, int32_t from, int32_t to, int32_t occupancy);
void resize (PMA* pma, Vertex *src);

int load_graph(string path);
void graph_init(int32_t max_vid);
void graph_maintenance_baseline(void);
void print_vertex(int32_t vid, int64_t &total_edges);
void print_graph(void);

// void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp);
void dumpToZone(void);
void dumpToZoneParallel(void);

void graph_init_new(int32_t max_vid);
void graph_maintenance_new(void);
void graph_maintenance_parallel(void);
void graph_maintenance_parallel2(void);

void graph_recovery(void);
int32_t binary_pma_find_recovery(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to);


PMA* pma_creat (void)
{
    PMA* pma = (PMA*)malloc(sizeof(PMA));
    pma->ele_num = 0;
    pma->arr_cap = 4;
    pma->seg_len = 2;
    pma->seg_num = pma->arr_cap / pma->seg_len;
    pma->tree_h = floor_lg (pma->seg_num) + 1;
    pma->delta_t = (t_0 - t_h) / pma->tree_h;  // 稠密度按线性变化
    pma->delta_p = (p_h - p_0) / pma->tree_h;
    // cout << "tree_h: " << pma->tree_h << " seg_num: " << pma->seg_num << endl;
    // pma->array = (element_type *)memkind_malloc(kind, sizeof(element_type) * pma->arr_cap);  // on NVM
    // pma->row_offset = (vertexOff *)malloc(sizeof(vertexOff) * pma->vertex_cap); // on DRAM
    return (pma);
}

static void compute_capacity (PMA* pma) 
{
    /*根据已有元素个数确定seg_len，arr_cap*/
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


int32_t binary_pma_find(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to) 
{
    int32_t index;
    int32_t blk_idx, bit_idx;
    // int32_t cur_blk;
    // int32_t *tmp_blk_ptr = new int32_t[PM_BLK_SIZE];  // 用于临时存储block数据
    vector<int32_t> release_blk_vec;

    while (from <= to) {
        int32_t mid = from + (to - from) / 2;
        int32_t i = mid;
        /* Start scanning left until we find a non-empty slot or we reach past the
        * beginning of the subarray. */
        blk_idx = i / PM_BLK_SIZE;
        bit_idx = i % PM_BLK_SIZE;

        if (src->blk_list[blk_idx].buf == nullptr) {
            readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
            release_blk_vec.push_back(blk_idx);
        }
        
        // cur_blk = blk_idx;
        // while (i >= from && (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)){
        while (i >= from && (src->blk_list[blk_idx].buf[bit_idx] < 0)){
            i--;
            blk_idx = i / PM_BLK_SIZE;
            bit_idx = i % PM_BLK_SIZE;
            // 如果blk更换了，咋需要先将对应blk数据读出来再进行处理
            if (src->blk_list[blk_idx].buf == nullptr) {
                readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
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
                    src->blk_list[release_blk_vec[jj]].buf = nullptr; // metadata用完后都赋为nullptr
                }
                // cout << "error: insert the same neighbor. src: " << src << " des: " << des << endl;
                return index;
            }
            else if (src->blk_list[blk_idx].buf[bit_idx] < des) {
                from = mid + 1;
            } else {  /* pma->array [i].key > key */
                to = i - 1;  // 确保最终找到的是小于des的最大元素位置
            }
        }
    }
    /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
    index = to;
    blk_idx = index / PM_BLK_SIZE;
    bit_idx = index % PM_BLK_SIZE;

    if (src->blk_list[blk_idx].buf == nullptr) {
        readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
        release_blk_vec.push_back(blk_idx);
    }

    // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block

    // cur_blk = blk_idx;
    while (index >= 0 && (src->blk_list[blk_idx].buf[bit_idx] < 0)) {
        index--;  // 确保找到的index必然存了值 也有可能是guard标签位
        blk_idx = index / PM_BLK_SIZE;
        bit_idx = index % PM_BLK_SIZE;
        // 如果blk更换了，咋需要先将对应blk数据读出来再进行处理

        if (src->blk_list[blk_idx].buf == nullptr) {
            readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
            release_blk_vec.push_back(blk_idx);
        }
    }

    for (int jj = 0; jj < release_blk_vec.size(); jj++) {
        delete[] src->blk_list[release_blk_vec[jj]].buf;
        src->blk_list[release_blk_vec[jj]].buf = nullptr; // metadata用完后都赋为nullptr
    }
    return index;
}

int32_t binary_pma_find_recovery(PMA* pma, Vertex *src, int32_t des, int32_t from, int32_t to) 
{
    int32_t index;
    int32_t blk_idx, bit_idx;
    // int32_t cur_blk;
    // int32_t *tmp_blk_ptr = new int32_t[PM_BLK_SIZE];  // 用于临时存储block数据
    vector<int32_t> release_blk_vec;

    while (from <= to) {
        int32_t mid = from + (to - from) / 2;
        int32_t i = mid;
        /* Start scanning left until we find a non-empty slot or we reach past the
        * beginning of the subarray. */
        blk_idx = i / PM_BLK_SIZE;
        bit_idx = i % PM_BLK_SIZE;

        if (src->blk_list[blk_idx].buf == nullptr) {
            readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
            release_blk_vec.push_back(blk_idx);
        }
        
        // cur_blk = blk_idx;
        // while (i >= from && (src->blk_list[blk_idx].blk_ptr[bit_idx] < 0)){
        while (i >= from && (src->blk_list[blk_idx].buf[bit_idx] < 0)){
            i--;
            blk_idx = i / PM_BLK_SIZE;
            bit_idx = i % PM_BLK_SIZE;
            // 如果blk更换了，咋需要先将对应blk数据读出来再进行处理
            if (src->blk_list[blk_idx].buf == nullptr) {
                readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
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
                    src->blk_list[release_blk_vec[jj]].buf = nullptr; // metadata用完后都赋为nullptr
                }
                // cout << "error: insert the same neighbor. src: " << src << " des: " << des << endl;
                return -1;
            }
            else if (src->blk_list[blk_idx].buf[bit_idx] < des) {
                from = mid + 1;
            } else {  /* pma->array [i].key > key */
                to = i - 1;  // 确保最终找到的是小于des的最大元素位置
            }
        }
    }
    /* Didn't find `key'. `to' should hold its predecessor (unless it's empty). */
    index = to;
    blk_idx = index / PM_BLK_SIZE;
    bit_idx = index % PM_BLK_SIZE;

    if (src->blk_list[blk_idx].buf == nullptr) {
        readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
        release_blk_vec.push_back(blk_idx);
    }

    // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block

    // cur_blk = blk_idx;
    while (index >= 0 && (src->blk_list[blk_idx].buf[bit_idx] < 0)) {
        index--;  // 确保找到的index必然存了值 也有可能是guard标签位
        blk_idx = index / PM_BLK_SIZE;
        bit_idx = index % PM_BLK_SIZE;
        // 如果blk更换了，咋需要先将对应blk数据读出来再进行处理

        if (src->blk_list[blk_idx].buf == nullptr) {
            readBlockFromZone(zns.dev_fd, src->blk_list[blk_idx].buf, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
            release_blk_vec.push_back(blk_idx);
        }
    }

    for (int jj = 0; jj < release_blk_vec.size(); jj++) {
        delete[] src->blk_list[release_blk_vec[jj]].buf;
        src->blk_list[release_blk_vec[jj]].buf = nullptr; // metadata用完后都赋为nullptr
    }
    return index;
}


int32_t pma_insert(PMA* pma, Vertex *src, int32_t des, int32_t index)
{
    /*
        insert after index // index 需要进行换算
        缓存相关block到DRAM中
    */
    int32_t j = index + 1;
    int32_t blk_idx, bit_idx;
    // int32_t cur_blk = -1;
    // int32_t *tmp_blk_ptr = new int32_t[PM_BLK_SIZE];

    while (j < pma->arr_cap) { // 向右查找空slot
        blk_idx = j / PM_BLK_SIZE;  // 读取bitmap代替读取NVM
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
                cacheBlock(src, blk_idx);  // 将对应blk缓存至DRAM中
            }

            int32_t pre_blk_idx = (j-1) / PM_BLK_SIZE;  // 将空位移动出来
            int32_t pre_bit_idx = (j-1) % PM_BLK_SIZE;
            if (src->blk_list[pre_blk_idx].buf == nullptr) {
                cacheBlock(src, pre_blk_idx);  // 将对应blk缓存至DRAM中
            }

            src->blk_list[blk_idx].buf[bit_idx] = src->blk_list[pre_blk_idx].buf[pre_bit_idx];
            j--;
        }

        blk_idx = (index + 1) / PM_BLK_SIZE;  // 计算插入位置
        bit_idx = (index + 1) % PM_BLK_SIZE;

        if (src->blk_list[blk_idx].buf == nullptr) {
            cacheBlock(src, blk_idx);  // 将对应blk缓存至DRAM中
        }
        src->blk_list[blk_idx].buf[bit_idx] = des;
        index = index + 1;  // 此处index表示元素插入的实际位置，后续rebalance的时候用
    } else { /* No empty space to the right of i. Try left. */
        // int32_t be_cur_blk = -1;
        // cout << "test --1 des: " << des << " index: " << index << endl;
        j = index - 1;
        while (j > 0) {
            blk_idx = j / PM_BLK_SIZE;  // 读取bitmap代替读取NVM
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
                blk_idx = j / PM_BLK_SIZE;  // 读取bitmap代替读取NVM
                bit_idx = j % PM_BLK_SIZE;

                if (src->blk_list[blk_idx].buf == nullptr) {
                    cacheBlock(src, blk_idx);  // 将对应blk缓存至DRAM中
                }

                int32_t be_blk_idx = (j+1) / PM_BLK_SIZE;  // 读取bitmap代替读取NVM
                int32_t be_bit_idx = (j+1) % PM_BLK_SIZE;

                if (src->blk_list[be_blk_idx].buf == nullptr) {
                    cacheBlock(src, be_blk_idx);  // 将对应blk缓存至DRAM中
                }

                src->blk_list[blk_idx].buf[bit_idx] = src->blk_list[be_blk_idx].buf[be_bit_idx];
                j++;
            }
            blk_idx = (index) / PM_BLK_SIZE; 
            bit_idx = (index) % PM_BLK_SIZE;

            if (src->blk_list[blk_idx].buf == nullptr) {
                cacheBlock(src, blk_idx);  // 将对应blk缓存至DRAM中
            }
            src->blk_list[blk_idx].buf[bit_idx] = des;
        }
    }
    pma->ele_num++;
    src->degree++;
    // print_graph();

    // cout << "test---before rebalance" << endl;
    
    rebalance(pma, src, index);

    return index;
}


void rebalance(PMA* pma, Vertex *src, int32_t index) 
{
    // 查找合适window的过程基于各block的bitmap即可
    int32_t window_start, window_end;
    int32_t height = 0;
    uint32_t occupancy = 1;
    int32_t left_index = index - 1;
    int32_t right_index = index + 1;
    double density, t_height, p_height;
    int32_t blk_idx, bit_idx;

    int32_t cur_blk = -1;
    int32_t *tmp_blk_ptr = new int32_t[PM_BLK_SIZE];
    
    do {
        uint32_t window_size = pma->seg_len * (1 << height);
        uint32_t window = index / window_size;
        window_start = window * window_size;
        window_end = window_start + window_size;

        while (left_index >= window_start) {
            blk_idx = left_index / PM_BLK_SIZE;
            bit_idx = left_index % PM_BLK_SIZE;

            if (src->blk_list[blk_idx].buf != nullptr) {
                if (src->blk_list[blk_idx].buf[bit_idx] >= 0) {
                    occupancy++;
                }
            } else {
                if (blk_idx != cur_blk) {
                    readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
                    cur_blk = blk_idx;
                }     

                if (tmp_blk_ptr[bit_idx] >= 0){
                    occupancy++;
                }
            }
            left_index--;
        }

        cur_blk = -1;
        while (right_index < window_end) {
            blk_idx = right_index / PM_BLK_SIZE;
            bit_idx = right_index % PM_BLK_SIZE;

            if (src->blk_list[blk_idx].buf != nullptr) {
                if (src->blk_list[blk_idx].buf[bit_idx] >= 0) {
                    occupancy++;
                }
            } else {
                if (blk_idx != cur_blk) {
                    readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
                    cur_blk = blk_idx;
                }     
                if (tmp_blk_ptr[bit_idx] >= 0) {
                    occupancy++;
                }
            }
            right_index++;
        }
        density = (double)occupancy / (double)window_size;
        t_height = t_0 - (height * pma->delta_t);
        p_height = p_0 + (height * pma->delta_p);
        height++;
    } while ((density < p_height ||  // 直至找到符合稠密度要求的Window density < p_height可删除
            density >= t_height) &&
           height < pma->tree_h);

    delete[] tmp_blk_ptr;
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


void spread (PMA* pma, Vertex *src, int32_t from, int32_t to, int32_t occupancy)
{
    /*
        将PMA上对应window放置在DRAM中平均分布
    */
    assert (from < to);
    int32_t capacity = to - from; // 不包含to
    int32_t frequency = (capacity << 8) / occupancy;  /* 8-bit fixed point arithmetic. */
    int32_t read_index = to - 1;  // 最后一个元素位置
    int32_t write_index = (to << 8) - frequency;

    int32_t start_blk = from / PM_BLK_SIZE;
    int32_t end_blk = (to - 1) / PM_BLK_SIZE;

    int32_t span_blks = end_blk - start_blk + 1;

    int32_t *tmp_win_arr = (int32_t*)malloc(capacity * sizeof(int32_t));

    if (tmp_win_arr != NULL) {
    // 使用循环逐个元素将数组初始化为 -1
        for (size_t i = 0; i < capacity; ++i) {
            tmp_win_arr[i] = -1;
        }
    }

    int32_t cur_blk = -1;
    int32_t *tmp_blk_ptr = new int32_t[PM_BLK_SIZE];

    if (span_blks == 1) {// 说明当前window在同一blk中
        int32_t base_off = start_blk * PM_BLK_SIZE;
        read_index = read_index - base_off;
        write_index = ((to-base_off) << 8) - frequency;  // 相当于blk内部绝对index
        from = from - base_off;

        if (src->blk_list[start_blk].buf == nullptr) {
            cacheBlock(src, start_blk);  // 将对应blk缓存至DRAM中
        }

        // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[start_blk].global_wp); // 从wp开始读取一个block
        // invalidBlock(src->blk_list[start_blk].global_wp);

        int32_t real_write_idx = (write_index >> 8) - from;
        occupancy -= 1; // 类似read_index初始绝对位置
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
        // copy tmp_win_arr to pma
        memcpy(&src->blk_list[start_blk].buf[from], tmp_win_arr, sizeof(int32_t) * capacity); // 将数组元素复制回PMA
        // 将blk写回zone
        // 必须申请个新blk用来存储该blk数据
        // if (graph_blk_que.size() == 0) {  // 不足一个blk
        //     generate_block_baseline(100);
        // }
        // pm_blk new_blk = graph_blk_que.front();
        // graph_blk_que.pop();
        // writeBlockToZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, new_blk.global_wp);
        // new_blk.src = src->id; // 用于GC
        // new_blk.blk_idx = start_blk;
        // // 修改vertex metadata
        // src->blk_list[start_blk] = new_blk;
        // validBlock(new_blk.global_wp);
        // 无效化被改写的blk

    } else {  // 当前window跨越多个block，挨个block进行处理
        occupancy -= 1; // 类似read_index初始绝对位置
        int32_t real_write_idx = (write_index >> 8) - from;
        for (int i = end_blk; i >= start_blk; i--) {
            int32_t loc_win_from, loc_win_to; // blk中的绝对地址

            if (src->blk_list[i].buf == nullptr) {
                cacheBlock(src, i);  // 将对应blk缓存至DRAM中
            }

            // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[i].global_wp); // 从wp开始读取一个block

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
        int32_t cpy_size = PM_BLK_SIZE - cpy_des;  // 元素数量 不是字节数
        // cout << "start blk: " << start_blk << " end blk: " << end_blk << endl; 
        for (int32_t i = start_blk; i <= end_blk; i++) {
            // readBlockFromZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, src->blk_list[i].global_wp); // 从wp开始读取一个block
            // invalidBlock(src->blk_list[i].global_wp);
            
            memcpy(&src->blk_list[i].buf[cpy_des], &tmp_win_arr[cpy_src], cpy_size * sizeof(int32_t));

            // 写入zone
            // if (graph_blk_que.size() == 0) {  // 不足一个blk
            //     generate_block_baseline(100);
            // }
            // pm_blk new_blk = graph_blk_que.front();
            // graph_blk_que.pop();
            // writeBlockToZone(zns.dev_fd, tmp_blk_ptr, PM_BLK_SIZE, new_blk.global_wp);
            // new_blk.src = src->id; // 用于GC
            // new_blk.blk_idx = i;
            // // 修改vertex metadata
            // src->blk_list[i] = new_blk;
            // validBlock(new_blk.global_wp);
            
            assert((cpy_des+cpy_size-1)<PM_BLK_SIZE);
            
            cpy_des = 0;
            cpy_src += cpy_size;

            if (i == (end_blk-1)) { // 为最后一个blk准备cpy_size
                cpy_size = ((to-1) % PM_BLK_SIZE) + 1;  // +1 是因为index从0开始，copy的数量要+1
            } else {
                cpy_size = PM_BLK_SIZE;
            }
        }
    }
    free(tmp_win_arr);
    // cout << "spread done" << endl;
}

void resize (PMA* pma, Vertex *src)
{
    int32_t read_index = pma->arr_cap - 1;  // 原先pma的末端开始读取
    int32_t readBitmapIdx = src->blk_list.size() - 1;
    // print_vertex(src);
    compute_capacity(pma);
    pma->tree_h = floor_lg (pma->seg_num) + 1;
    pma->delta_t = (t_0 - t_h) / pma->tree_h;
    pma->delta_p = (p_h - p_0) / pma->tree_h;


    // 根据现有blk数量和pma容量计算所需额外blk数量
    int32_t blk_num = ceil((double)pma->arr_cap / PM_BLK_SIZE);  // 所需总的blk数量
    int32_t ex_blk_num = blk_num - src->blk_list.size();
    

    for (int32_t i=0; i < ex_blk_num; i++) {
        pm_blk new_blk;
        new_blk.buf = new int32_t[PM_BLK_SIZE];
        fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
        src->blk_list.push_back(new_blk);
        global_cache_blk++;
    }
    // cout << "ex_blk_num: " << ex_blk_num << " src->blk_list.size: " << src->blk_list.size() << endl;

    // 申请一片内存，用于存储blk_num个block数据
    // int32_t total_slots = blk_num * PM_BLK_SIZE;
    // int32_t *tmp_win_arr = (int32_t*)malloc(total_slots * sizeof(int32_t));
    // if (tmp_win_arr != NULL) {
    //     for (size_t i = 0; i < total_slots; ++i) { // 使用循环逐个元素将数组初始化为 -1
    //         tmp_win_arr[i] = -1;
    //     }
    // }

    int32_t capacity = pma->arr_cap;
    int32_t frequency = (capacity << 8) / pma->ele_num;
    int32_t write_index = (capacity << 8) - frequency;

    int32_t writeBitmapIdx = src->blk_list.size() - 1; // 最后一个blk 索引

    // int32_t cur_read_blk = -1;
    // int32_t cur_write_blk = -1;
    // int32_t *tmp_read_blk_ptr = new int32_t[PM_BLK_SIZE];
    // int32_t *tmp_write_blk_ptr = new int32_t[PM_BLK_SIZE];
    // cout << "test--0" << endl;
    
    while ((write_index >> 8) > read_index && read_index >= 0)
    {
        // cout << "read_idx: " << read_index << " write_index: " << (write_index >> 8) 
        //     << " readBitmapIdx: " << readBitmapIdx << " writeBitmapIdx: " << writeBitmapIdx << endl;
        int32_t write_off = (write_index >> 8) % PM_BLK_SIZE;
        int32_t read_off = read_index % PM_BLK_SIZE;

        if (src->blk_list[readBitmapIdx].buf == nullptr) {
            cacheBlock(src, readBitmapIdx); 
        }

        if (src->blk_list[readBitmapIdx].buf[read_off] >= 0) {
        // if (testBit(&src.blk_list[readBitmapIdx].myBitmap, read_off)) {
            // cout << "read_off: " << read_off << " write_off: " << write_off << endl;
            src->blk_list[writeBitmapIdx].buf[write_off] = src->blk_list[readBitmapIdx].buf[read_off];
            // cout << "test ---0 " << endl;
            // clearBit(&src.blk_list[readBitmapIdx].myBitmap, read_off);
            src->blk_list[readBitmapIdx].buf[read_off] = -1;
            // setBit(&src.blk_list[writeBitmapIdx].myBitmap, write_off);
            read_index--;
            write_index -= frequency;

            
            if ((read_index % PM_BLK_SIZE) > read_off) {
                readBitmapIdx --;
                // cout << "readBitmapIdx: " << readBitmapIdx << endl;
            }
            if (((write_index >> 8) % PM_BLK_SIZE) > write_off) {
                writeBitmapIdx--;
                // cout << "writeBitmapIdx: " << writeBitmapIdx << endl;
            }
        } else {
            // cout << "test ---1 " << endl;
            read_index--;
            if ((read_index % PM_BLK_SIZE) > read_off) {
                readBitmapIdx --;
                // cout << "readBitmapIdx: " << readBitmapIdx << endl;
            }
        }
    }
    // print_graph();
}

int load_graph(string path)
{
    std::ifstream inputFile(path);

    if (!inputFile) {
        std::cerr << "Failed to open file." << std::endl;
        return 1;
    }
    cout << "Load dataset " << path << endl;

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
    return 1;
}

void graph_init(int32_t max_vid)
{
	generate_block_baseline(9000);
	global_vectex_vec.resize(max_vid + 1);

    // cout << "test ---0" << endl;
	for (int32_t vidx = 0; vidx <= max_vid; vidx ++) {
        // cout << "vidx: " << vidx << " global_blk_vec size:" << global_blk_vec.size() << endl;
		Vertex* v = new Vertex(vidx);
        // 每个顶点初始化一个block buf，用于存放邻节点
        pm_blk new_blk;
        new_blk.buf = new int32_t[PM_BLK_SIZE]; // 在DRAM中分配了内存
        fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
        v->blk_list.push_back(new_blk);

		v->pma = pma_creat();
        global_vectex_vec[vidx] = v;
        global_cache_blk++;
	}
    cout << "graph initial complete! total cache blk: " << global_cache_blk << endl;
}

void graph_init_new(int32_t max_vid)
{
	generate_block_baseline((zns.nr_blocks_perzone * global_num_queues));
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

void graph_init_new2(int32_t max_vid)
{
	generate_block_opt();  // 对每个zone中的blk进行单独管理，可实现写zone
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

void graph_maintenance_baseline(void)
{
    graph_init(global_max_vid);
    cout << "Load graph on NVM baseline begin!" << endl;
	double start_t = get_current_time();

	for (int64_t edgeIdx = 0; edgeIdx < global_edge_num; edgeIdx++)
	{
        new_edge e = global_edge_vec[edgeIdx];  // 从nvm上读取
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

void graph_maintenance_new(void)
{
    /*
        相对baseline不同的是，首先分配了缓存行，degree小的顶点不需要再存储于PMA中
    */
    graph_init_new(global_max_vid);
    cout << "Load graph on NVM new begin!" << endl;
	double start_t = get_current_time();

	for (int64_t edgeIdx = 0; edgeIdx < global_edge_num; edgeIdx++)
	{
        new_edge e = global_edge_vec[edgeIdx];  // 从nvm上读取
		int32_t src = e.src;
		int32_t des = e.des;
		Vertex* v = global_vectex_vec[src];
		int32_t index;
        // cout << "insert edge " << edgeIdx << ": (" << src << ", " << des << ")" << endl; 

        if (v->degree < CACHE_BLK_SIZE) { // 直接插入新增边，不需要排序
            // cout << "test-- cache insert" << endl;
            v->cache_buf[v->degree] = des;
            v->degree++;
        } else if (v->degree == CACHE_BLK_SIZE) { // 需要将之前缓存的边和最新边都插入
            // cout << "test-- boardline insert" << endl;
            pm_blk new_blk;
            new_blk.buf = new int32_t[PM_BLK_SIZE]; // 在DRAM中分配了内存
            fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
            v->blk_list.push_back(new_blk);
            global_cache_blk++;
            v->degree = 0;
            for (int i = 0; i < CACHE_BLK_SIZE; i++) {
                des = v->cache_buf[i];
                index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
                pma_insert(v->pma, v, des, index);
            }
            des = e.des;  // 把最新的边插进去
            index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
            pma_insert(v->pma, v, des, index);
        } else {
            // cout << "test-- PMA insert" << endl;
            index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
            pma_insert(v->pma, v, des, index);
            if (global_cache_blk >= cache_blk_threshold) {
                cout << "dump to zone" << endl;
                dumpToZone();
                global_cache_blk = 0;
            }
        }
    }
    // print_vertex(global_vectex_vec[1]);
    // print_pma(global_vectex_vec[1]);
    
	double end_t = get_current_time();
    // print_graph();
    // print_pma(global_vectex_vec[1]);
	cout << "Load graph on NVM new time: " << end_t - start_t << endl;
    // print_graph();
}

void graph_maintenance_parallel(void)
{
    /*
        相对baseline不同的是，首先分配了缓存行，degree小的顶点不需要再存储于PMA中
        dumpToZone没做并发
    */
    graph_init_new(global_max_vid);
    cout << "Load graph on NVM new begin!" << endl;
	double start_t = get_current_time();

    // 每个batch插入完成后检查是否DRAM中cache的block是否达到阈值进行dump
    omp_set_num_threads(THREAD_NUM);
    int64_t start_e = 0;
    int64_t end_e = min(start_e + global_batch_size, global_edge_num);
    // global_batch_size = global_edge_num;
    while (end_e <= global_edge_num) {
        #pragma omp parallel for
        for (int64_t edgeIdx = start_e; edgeIdx < end_e; edgeIdx++)
        {
            new_edge e = global_edge_vec[edgeIdx];  // 从nvm上读取
            int32_t src = e.src;
            int32_t des = e.des;
            Vertex* v = global_vectex_vec[src];
            unique_lock<std::shared_timed_mutex> lock(v->mutex); // (已验证ok)确保同一顶点同时只处理一条边
            int32_t index;
            // cout << "insert edge " << edgeIdx << ": (" << src << ", " << des << ")" << endl; 

            if (v->degree < CACHE_BLK_SIZE) { // 直接插入新增边，不需要排序
                // cout << "test-- cache insert" << endl;
                v->cache_buf[v->degree] = des;
                v->degree++;
            } else if (v->degree == CACHE_BLK_SIZE) { // 需要将之前缓存的边和最新边都插入
                // cout << "test-- boardline insert" << endl;
                pm_blk new_blk;
                new_blk.buf = new int32_t[PM_BLK_SIZE]; // 在DRAM中分配了内存
                fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
                v->blk_list.push_back(new_blk);
                global_cache_blk++;
                v->degree = 0;
                for (int i = 0; i < CACHE_BLK_SIZE; i++) {
                    des = v->cache_buf[i];
                    index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
                    pma_insert(v->pma, v, des, index);
                }
                des = e.des;  // 把最新的边插进去
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

void graph_maintenance_parallel2(void)
{
    /*
        相对baseline不同的是，首先分配了缓存行，degree小的顶点不需要再存储于PMA中
        dumpToZone没做并发
    */
    graph_init_new2(global_max_vid);
    cout << "Load graph on NVM new begin!" << endl;
	double start_t = get_current_time();

    // 每个batch插入完成后检查是否DRAM中cache的block是否达到阈值进行dump
    omp_set_num_threads(THREAD_NUM);
    int64_t start_e = 0;
    global_batch_size = global_edge_num;
    int64_t end_e = min(start_e + global_batch_size, global_edge_num);
    
    while (end_e <= global_edge_num) {
        // cout << "end_e: " << end_e << endl;
        #pragma omp parallel for
        for (int64_t edgeIdx = start_e; edgeIdx < end_e; edgeIdx++)
        {
            new_edge e = global_edge_vec[edgeIdx];  // 从nvm上读取
            int32_t src = e.src;
            int32_t des = e.des;
            Vertex* v = global_vectex_vec[src];
            unique_lock<std::shared_timed_mutex> lock(v->mutex); // (已验证ok)确保同一顶点同时只处理一条边
            int32_t index;
            // cout << "insert edge " << edgeIdx << ": (" << src << ", " << des << ")" << endl; 

            if (v->degree < CACHE_BLK_SIZE) { // 直接插入新增边，不需要排序
                // cout << "test-- cache insert" << endl;
                v->cache_buf[v->degree] = des;
                v->degree++;
            } else if (v->degree == CACHE_BLK_SIZE) { // 需要将之前缓存的边和最新边都插入
                // cout << "test-- boardline insert" << endl;
                pm_blk new_blk;
                new_blk.buf = new int32_t[PM_BLK_SIZE]; // 在DRAM中分配了内存
                fill(new_blk.buf, new_blk.buf + PM_BLK_SIZE, -1);
                v->blk_list.push_back(new_blk);
                global_cache_blk++;
                v->degree = 0;
                for (int i = 0; i < CACHE_BLK_SIZE; i++) {
                    des = v->cache_buf[i];
                    index = binary_pma_find(v->pma, v, des, 0, v->pma->arr_cap-1);
                    pma_insert(v->pma, v, des, index);
                }
                des = e.des;  // 把最新的边插进去
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
        return;
    }
    // cout << "vid: " << vid << " blknum: " << v->blk_list.size() << " {";

    if (v->degree <= CACHE_BLK_SIZE) {
        total_edges += v->degree;
        return;
    }
    
    for (int i = 0; i < v->blk_list.size(); i++) {
        if (v->blk_list[i].buf == nullptr) {
            // cout << "test ==" << endl;
            // cout << "read wp: " << v->blk_list[i].global_wp << " blk: " << (v->blk_list[i].global_wp/ZONE_BLK_SIZE)
            // << " status: " << test_block(0, (v->blk_list[i].global_wp/ZONE_BLK_SIZE)) << endl;
            readBlockFromZone(zns.dev_fd, blk_ptr, PM_BLK_SIZE, v->blk_list[i].global_wp); // 从wp开始读取一个block
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

// // 记录zone中blk在对应顶点blk_list中的位置，用于回溯时精准定位并更改vertex metadata
// void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp)
// {
//     // cout << "zone size: " << zns.info.zone_size << endl;
//     int32_t zone_idx = wp / zns.info.zone_size; // 第几个zone
//     int32_t info_idx = (wp % zns.info.zone_size) / ZONE_BLK_SIZE; // blk_info_vec的index
//     global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status = 1;
//     global_zone_vec[zone_idx]->blk_info_vec[blk_idx].src = vid;
//     global_zone_vec[zone_idx]->blk_info_vec[blk_idx].blkIndex = blk_idx;
// }

// 当DRAM中cache的block达到阈值后将数据写入zone，只写入block数量较小的顶点的block
void dumpToZone(void)
{
    // int total_blks = 0;
    double start_t = get_current_time();
    for (int32_t vidx = 0; vidx < global_vectex_vec.size(); vidx++) {
        Vertex *v = global_vectex_vec[vidx];
        // cout << "dump vid: " << v->id << " degree: " << v->degree << " blks: " << v->blk_list.size() << endl;
        if ((v->degree > 0) && (v->blk_list.size() > 0)) {
            for (int blkidx = 0; blkidx < v->blk_list.size(); blkidx++) {
                if (v->blk_list[blkidx].buf == nullptr) { // 该blk没有cache在DRAM中
                    continue;
                }
                // total_blks++;
                // cout << "dump vid: " << v->id << endl;

                if (graph_blk_que.size() == 0) {
                    garbage_collection();
                }
                
                pm_blk new_blk = graph_blk_que.front();
                graph_blk_que.pop();
                writeBlockToZone(zns.dev_fd, v->blk_list[blkidx].buf, PM_BLK_SIZE, new_blk.global_wp);
                delete[] v->blk_list[blkidx].buf;
                mark_zone(v->id, blkidx, new_blk.global_wp);  // 用于GC
                // 修改vertex metadata
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
                if (v->blk_list[blkidx].buf == nullptr) { // 该blk没有cache在DRAM中
                    continue;
                }

                // if (graph_blk_que.size() == 0) {
                //     garbage_collection();
                // }

                unsigned long long wp = requestMemoryBlocks(threadId); 
                // cout << "dump wp: " << wp << endl;
                
                writeBlockToZone(zns.dev_fd, v->blk_list[blkidx].buf, PM_BLK_SIZE, wp);
                // cout << " write to zone done." << endl;
                delete[] v->blk_list[blkidx].buf;
                mark_zone(v->id, blkidx, wp);  // 用于GC
                // 修改vertex metadata
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
            wp += ZONE_BLK_SIZE;
        }
    }
    cout << "test --0" << endl;
    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for
    for (int64_t edgeIdx = (global_edge_num-1); edgeIdx >= 0; edgeIdx--) {
        new_edge e = global_edge_vec[edgeIdx];  // 从nvm上读取
        int32_t src = e.src;
        int32_t des = e.des;
        Vertex* v = global_vectex_vec[src];
        unique_lock<std::shared_timed_mutex> lock(v->mutex); // (已验证ok)确保同一顶点同时只处理一条边
        int32_t index;
        if (v->degree <= CACHE_BLK_SIZE) { // 直接插入新增边，不需要排序
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
