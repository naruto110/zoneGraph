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
#define CACHE_BLK_SIZE 16 // 小度顶点直接缓存在 DRAM 中的邻居数量
#define ZONE_BLK_SIZE (PM_BLK_SIZE * sizeof(int32_t))  // 一个逻辑块的字节数
#define NVME_DEVICE_PATH "/dev/nvme0n1"

int THREAD_NUM = 64;

// struct zns zns;

typedef struct _pma
{
    int32_t ele_num;  // PMA 中当前有效元素个数
    int32_t arr_cap;  // PMA 总容量（按 slot 计）
    int32_t seg_len;  // 叶子段长度
    int32_t seg_num;  // 段数量
    int32_t tree_h;   // PMA 树高
    double delta_t;   // 上界密度阈值的层级变化量
    double delta_p;   // 下界密度阈值的层级变化量
} PMA;

typedef struct pma_blk
{
	// int32_t blk_id;
	// int zone_idx;
	// unsigned long long local_wp;
	unsigned long long global_wp; // 该块在设备上的全局偏移

	// int32_t src;
	// int32_t blk_idx;
	int32_t *buf; // DRAM 中的块缓存；非空表示该块当前已换入内存
} pm_blk;

typedef struct new_edge
{
	int32_t src;  // 边的源点
	int32_t des;  // 边的目标点
	// bool statue;
	// char type; // edge insertion, graph processing task...
} newEdge;

class Vertex {
public:
	int32_t id;
	int32_t degree; // 当前顶点邻居数
	vector<pm_blk> blk_list; // PMA 对应的块列表
	PMA* pma; // 顶点的 PMA 元数据
	int32_t cache_buf[CACHE_BLK_SIZE]; // 小度顶点的快速缓存
	mutable std::shared_timed_mutex mutex; // 顶点级锁，用于并发插入/恢复

	bool recovery_flag; // 恢复过程中避免重复处理

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
	int nr_blocks; // zone 内总块数
	int val_blocks; // zone 内有效块数
	vector<blk_info> blk_info_vec; // 每个物理块对应的逻辑映射信息

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

    // 申请下一个可写块地址
    bool getNextBlock(unsigned long long& wp) {
        lock_guard<std::mutex> lock(mtx_);
        if (memoryWPs.empty()) {
            return false; // 当前队列没有可用块
        }
        wp = memoryWPs.front();
        memoryWPs.pop();
        return true;
    }

    // 预留的回收接口，当前未真正实现
    void reclaimMemory() {
        lock_guard<std::mutex> lock(mtx_);
        // cout << "Reclaiming memory and refilling the queue." << std::endl;
        // clear();
        // refill();
    }

// private:
    queue<unsigned long long> memoryWPs; // 存储该 zone 内尚未分配的物理块地址
    mutex mtx_; // 保护队列访问的互斥锁

    // 用一个 zone 的全部块地址初始化队列
    void refill(unsigned long long wp) {
		for (int i=0; i < nr_blocks_perzone; i++) {
			memoryWPs.push(wp);
			wp += ZONE_BLK_SIZE;
		}
    }

    // 这里不负责释放设备空间，保留空实现
    void clear() {

    }
};


typedef struct zone_information
{
	// Bitmap32 myBitmap;  // 旧版本中用于标记块状态的位图
	int nr_blocks; // zone 内总块数
	int val_blocks; // zone 内有效块数
	vector<blk_info> blk_info_vec; // 每个物理块对应的逻辑映射信息
} zone_info;

queue<pm_blk> graph_blk_que; // 基线模式下的空闲块队列
vector<MemoryWPQueue*> global_blk_queues; // 并行模式下按 zone 划分的可写块队列
/*
	每个队列对应一个 zone。
	队列中的 memoryWPs 保存该 zone 内尚未使用的块起始地址。
*/

unsigned int global_nr_zones = 0;  // zone 总数
unsigned int global_used_zones = 0;  // 已使用 zone 数
unsigned long long global_threshold_wp;  // 预留阈值位置，目前未完整使用
unsigned long long global_blk_wp = 0; // 基线模式下当前分配到的全局写指针
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
		预留的全局初始化入口，当前未使用。
	*/
}

#endif
