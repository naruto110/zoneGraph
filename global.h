#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <libzbd/zbd.h>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include <filesystem>
#include <vector>
#include "bitmap.h"
#include <queue>
#include "utils.h"
#include "tools.h"
#include <cmath>
#include <algorithm> 
#include <fstream>
#include <omp.h>
#include <shared_mutex>
#include <mutex>
#include <atomic>
#include <random>
// #include <benchmark.h>

using namespace std;

#define PM_BLK_SIZE 64
#define CACHE_BLK_SIZE 16 // 一个cacheline存放的邻节点数
#define ZONE_BLK_SIZE (PM_BLK_SIZE * sizeof(int32_t))  // 每个blk的字节数
#define NVME_DEVICE_PATH "/dev/nullb0"

int THREAD_NUM = 64;

// struct zns zns;

typedef struct _pma
{
    int32_t ele_num;  // pma_array中元素个数
    int32_t arr_cap;  // pma_array size
    int32_t seg_len;  // segment length
    int32_t seg_num;  // segment 数量
    int32_t tree_h;   // tree height
    double delta_t;  // Delta for upper density threshold.
    double delta_p;  // Delta for lower density threshold
} PMA;

typedef struct pma_blk
{
	// int32_t blk_id;
	// int zone_idx; // blk起始地址，nvm上的地址
	// unsigned long long local_wp;  // write pointer, 局部偏移量
	unsigned long long global_wp; // write pointer, 全局偏移量

	// int32_t src; // 用于存储哪个src顶点的PMA
	// int32_t blk_idx; //在当前src neighbor list中的blk index，方便直接替换
	int32_t *buf; // 为了避免写一个数就要重复制整个乃至多个blk（巨大写当大），先将用到的buf缓存在DRAM中，当DRAM中的缓存blk足够多，再将数据写入zone中
} pm_blk;

typedef struct new_edge
{
	int32_t src;  // blk在graph_blk_vec中的idex下标
	int32_t des;  // 最大degree的顶点
	// bool statue;
	// char type; // edge insertion, graph processing task...
} newEdge;

class Vertex {
public:
	int32_t id;
	int32_t degree;
	// mutable std::shared_timed_mutex mutex;
	vector<pm_blk> blk_list;
	PMA* pma;
	int32_t cache_buf[CACHE_BLK_SIZE];
	mutable std::shared_timed_mutex mutex;

	bool recovery_flag;

	Vertex(int32_t vid){
		id = vid;
		degree = 0;
		recovery_flag = false;
		// cache_buf = new int32_t[CACHE_BLK_SIZE];
		// fill(cache_buf, cache_buf + CACHE_BLK_SIZE, -1);
	}

	~Vertex(){}
};

typedef struct block_information
{
	bool status;
	int32_t src;
	int32_t blkIndex; 
} blk_info;

class ZoneStatus {
public:
	int nr_blocks; // 总blk数量
	int val_blocks; // 有效blk数量
	vector<blk_info> blk_info_vec; // resize成nr_block大小

	ZoneStatus(){
		nr_blocks = 0;
		val_blocks = 0;
	}

	ZoneStatus(int nb, int vb){
		nr_blocks = nb;
		val_blocks = vb;
	}

	~ZoneStatus(){}
};

int nr_blocks_perzone = 0;

class MemoryWPQueue {
public:
    MemoryWPQueue() {
        // refill();
    }

	MemoryWPQueue(unsigned long long wp) {
        refill(wp);
    }

    ~MemoryWPQueue() {
        clear();
    }

    // 获取下一个内存块
    bool getNextBlock(unsigned long long& wp) {
        lock_guard<std::mutex> lock(mtx_);
        if (memoryWPs.empty()) {
            return false; // 没有可用的内存块
        }
        wp = memoryWPs.front();
        memoryWPs.pop();
        return true;
    }

    // 垃圾回收函数
    void reclaimMemory() {
        lock_guard<std::mutex> lock(mtx_);
        // 模拟释放内存并重新填充队列
        // cout << "Reclaiming memory and refilling the queue." << std::endl;
        // clear();
        // refill();
    }

// private:
    queue<unsigned long long> memoryWPs; // 存储内存块的队列,zone中的物理地址
    mutex mtx_; // 保护对队列的访问

    // 每个zone划分成多个block块
    void refill(unsigned long long wp) {
		for (int i=0; i < nr_blocks_perzone; i++) {
			memoryWPs.push(wp);
			wp += ZONE_BLK_SIZE;
		}
    }

    // 释放内存块
    void clear() {

    }
};


typedef struct zone_information
{
	// Bitmap32 myBitmap;  // 用于标记block状态，1:valid，0：invalid
	int nr_blocks; // 总blk数量
	int val_blocks; // 有效blk数量
	vector<blk_info> blk_info_vec; // resize成nr_block大小
} zone_info;

queue<pm_blk> graph_blk_que; // 存储备用blk，串行分配zone blks
vector<MemoryWPQueue*> global_blk_queues; // 用于并发分配zone blocks
/*
	queue中包含多个MemoryWPQueue对象，其中有个成员变量memoryWPs存储每个block的起始地址
	每个zone对应一个MemoryWPQueue对象，用来维护其空间分配
*/

unsigned int global_nr_zones = 0;  // 总的zone数量
unsigned int global_used_zones = 0;  // 被分配完的zone数量
unsigned long long global_threshold_wp;  // 需要初始化为80%总空间大位置
unsigned long long global_blk_wp = 0; // 当前分配到的wp地址
int32_t global_blk_id = 0;
// int32_t global_cache_blk = 0;
atomic<int> global_cache_blk(0);
int32_t cache_blk_threshold = 1;
int64_t global_edge_num = 0;
int32_t global_max_vid = 0;
vector<Vertex *> global_vectex_vec;
vector<new_edge> global_edge_vec;
vector<ZoneStatus *> global_zone_vec;
int32_t global_empty_zone = 0;

int32_t global_num_queues = 0;
int64_t global_batch_size = 10000000;  
double total_dump_time = 0.0;

void initial() { 
	/*
		初始化global参数
	*/
	
	
}

#endif
