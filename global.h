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
#define CACHE_BLK_SIZE 16 // һ��cacheline��ŵ��ڽڵ���
#define ZONE_BLK_SIZE (PM_BLK_SIZE * sizeof(int32_t))  // ÿ��blk���ֽ���
#define NVME_DEVICE_PATH "/dev/nullb0"

int THREAD_NUM = 64;

// struct zns zns;

typedef struct _pma
{
    int32_t ele_num;  // pma_array��Ԫ�ظ���
    int32_t arr_cap;  // pma_array size
    int32_t seg_len;  // segment length
    int32_t seg_num;  // segment ����
    int32_t tree_h;   // tree height
    double delta_t;  // Delta for upper density threshold.
    double delta_p;  // Delta for lower density threshold
} PMA;

typedef struct pma_blk
{
	// int32_t blk_id;
	// int zone_idx; // blk��ʼ��ַ��nvm�ϵĵ�ַ
	// unsigned long long local_wp;  // write pointer, �ֲ�ƫ����
	unsigned long long global_wp; // write pointer, ȫ��ƫ����

	// int32_t src; // ���ڴ洢�ĸ�src�����PMA
	// int32_t blk_idx; //�ڵ�ǰsrc neighbor list�е�blk index������ֱ���滻
	int32_t *buf; // Ϊ�˱���дһ������Ҫ�ظ��������������blk���޴�д���󣩣��Ƚ��õ���buf������DRAM�У���DRAM�еĻ���blk�㹻�࣬�ٽ�����д��zone��
} pm_blk;

typedef struct new_edge
{
	int32_t src;  // blk��graph_blk_vec�е�idex�±�
	int32_t des;  // ���degree�Ķ���
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
	int nr_blocks; // ��blk����
	int val_blocks; // ��Чblk����
	vector<blk_info> blk_info_vec; // resize��nr_block��С

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

    // ��ȡ��һ���ڴ��
    bool getNextBlock(unsigned long long& wp) {
        lock_guard<std::mutex> lock(mtx_);
        if (memoryWPs.empty()) {
            return false; // û�п��õ��ڴ��
        }
        wp = memoryWPs.front();
        memoryWPs.pop();
        return true;
    }

    // �������պ���
    void reclaimMemory() {
        lock_guard<std::mutex> lock(mtx_);
        // ģ���ͷ��ڴ沢����������
        // cout << "Reclaiming memory and refilling the queue." << std::endl;
        // clear();
        // refill();
    }

// private:
    queue<unsigned long long> memoryWPs; // �洢�ڴ��Ķ���,zone�е������ַ
    mutex mtx_; // �����Զ��еķ���

    // ÿ��zone���ֳɶ��block��
    void refill(unsigned long long wp) {
		for (int i=0; i < nr_blocks_perzone; i++) {
			memoryWPs.push(wp);
			wp += ZONE_BLK_SIZE;
		}
    }

    // �ͷ��ڴ��
    void clear() {

    }
};


typedef struct zone_information
{
	// Bitmap32 myBitmap;  // ���ڱ��block״̬��1:valid��0��invalid
	int nr_blocks; // ��blk����
	int val_blocks; // ��Чblk����
	vector<blk_info> blk_info_vec; // resize��nr_block��С
} zone_info;

queue<pm_blk> graph_blk_que; // �洢����blk�����з���zone blks
vector<MemoryWPQueue*> global_blk_queues; // ���ڲ�������zone blocks
/*
	queue�а������MemoryWPQueue���������и���Ա����memoryWPs�洢ÿ��block����ʼ��ַ
	ÿ��zone��Ӧһ��MemoryWPQueue��������ά����ռ����
*/

unsigned int global_nr_zones = 0;  // �ܵ�zone����
unsigned int global_used_zones = 0;  // ���������zone����
unsigned long long global_threshold_wp;  // ��Ҫ��ʼ��Ϊ80%�ܿռ��λ��
unsigned long long global_blk_wp = 0; // ��ǰ���䵽��wp��ַ
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
		��ʼ��global����
	*/
	
	
}

#endif
