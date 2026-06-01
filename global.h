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
#include <queue>
#include "utils.h"
#include "tools.h"
#include <cmath>
#include <algorithm> 
#include <fstream>
#include <omp.h>
#include <shared_mutex>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <deque>
#include <atomic>
#include <random>
#include <cstdint>
#include <climits>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
// #include <benchmark.h>

using namespace std;

#define PM_BLK_SIZE 64
#define CACHE_BLK_SIZE 16 // 小度顶点直接缓存在 DRAM 中的邻居数量
#define ZONE_BLK_SIZE (PM_BLK_SIZE * sizeof(int32_t))  // 一个逻辑块的字节数
#define NVME_DEVICE_PATH "/dev/nvme0n1"
#define ZNS_BLOCK_MAGIC 0x5A47504D
#define ZNS_BLOCK_VERSION 3
#define ZNS_RECORD_SIZE 512
#define ZNS_RECORD_TAIL_MAGIC 0x504D475A
#define ZNS_REDO_MAGIC 0x5A475244
#define ZNS_REDO_TAIL_MAGIC 0x4452475A
#define ZNS_REDO_PAGE_MAGIC 0x5A475250
#define ZNS_REDO_PAGE_TAIL_MAGIC 0x5052475A

enum ZNSRecordType : uint32_t {
	ZNS_RECORD_DATA = 1,
	ZNS_RECORD_COMMIT = 2
};

enum EdgeLogOpType : uint32_t {
	EDGE_LOG_INSERT = 1,
	EDGE_LOG_DELETE = 2
};

struct ZNSBlockHeader {
	uint32_t magic;
	uint32_t version;
	uint32_t record_type;
	uint32_t payload_bytes;
	uint64_t epoch;
	int32_t src_vid;
	int32_t pma_blk_idx;
	uint64_t seq;
	uint64_t first_seq;
	uint64_t last_seq;
	uint32_t data_record_count;
	uint32_t aggregate_checksum;
	uint32_t valid_count;
	uint32_t pma_ele;
	uint32_t pma_cap;
	uint64_t committed_op_id;
	uint32_t checksum;
	uint32_t tail_magic;
};

struct EdgeLogRecord {
	uint32_t magic;
	uint32_t version;
	uint32_t op_type;
	uint32_t reserved;
	int32_t src;
	int32_t dst;
	uint32_t checksum;
	uint32_t tail_magic;
};

struct EdgeLogEntry {
	int32_t src;
	int32_t dst;
	uint32_t op_type;
};

struct EdgeLogPageHeader {
	uint32_t magic;
	uint32_t version;
	uint32_t record_count;
	uint32_t payload_bytes;
	uint64_t first_op_id;
	uint64_t last_op_id;
	uint32_t checksum;
	uint32_t tail_magic;
};

static const uint32_t ZNS_REDO_ENTRIES_PER_PAGE =
	(ZNS_RECORD_SIZE - sizeof(EdgeLogPageHeader)) / sizeof(EdgeLogEntry);

static inline size_t zns_alloc_unit_size(void);

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
	unsigned long long global_wp = 0; // 该块在设备上的全局偏移

	// int32_t src;
	// int32_t blk_idx;
	int32_t *buf = nullptr; // DRAM 中的块缓存；非空表示该块当前已换入内存
	int32_t min_val = INT32_MAX; // 当前块内有效邻居最小值，仅内存态
	int32_t max_val = INT32_MIN; // 当前块内有效邻居最大值，仅内存态
	int32_t valid_count = 0; // 当前块内非 -1 slot 数
	bool sorted = true; // 当前块内有效元素按 slot 顺序非降序，可用二分/有序扫描
} pm_blk;

static inline void reset_pm_blk_meta(pm_blk &blk)
{
	blk.min_val = INT32_MAX;
	blk.max_val = INT32_MIN;
	blk.valid_count = 0;
	blk.sorted = true;
}

static inline void refresh_pm_blk_meta(pm_blk &blk, const int32_t *data)
{
	reset_pm_blk_meta(blk);
	if (data == nullptr) {
		return;
	}
	int32_t prev = INT32_MIN;
	for (int32_t i = 0; i < PM_BLK_SIZE; i++) {
		if (data[i] >= 0) {
			blk.valid_count++;
			blk.min_val = min(blk.min_val, data[i]);
			blk.max_val = max(blk.max_val, data[i]);
			if (data[i] < prev) {
				blk.sorted = false;
			}
			prev = data[i];
		}
	}
}

typedef struct new_edge
{
	int32_t src;  // 边的源点
	int32_t des;  // 边的目标点
	// bool statue;
	// char type; // edge insertion, graph processing task...
} newEdge;

enum VertexClass {
	VERTEX_COLD = 0,
	VERTEX_HOT = 1,
	VERTEX_HUB = 2
};

enum ATZVertexType {
	ATZ_VERTEX_COLD = 0,
	ATZ_VERTEX_WARM = 1,
	ATZ_VERTEX_HOT = 2
};

enum ATZZoneType {
	ATZ_ZONE_FREE = 0,
	ATZ_ZONE_MIXED = 1,
	ATZ_ZONE_COLD = 2,
	ATZ_ZONE_WARM = 3,
	ATZ_ZONE_HOT = 4
};

class Vertex {
public:
	int32_t id;
	int32_t degree; // 当前顶点邻居数
	vector<pm_blk> blk_list; // PMA 对应的块列表
		PMA* pma; // 顶点的 PMA 元数据
		int32_t cache_buf[CACHE_BLK_SIZE]; // 小度顶点的快速缓存
		unordered_map<int32_t, vector<int32_t>> pending_insert_by_block; // 按 PMA block 分组的 append-only DRAM 插入缓冲
		vector<int32_t> pending_insert_unplaced; // 暂时无法绑定到稳定 block 的插入缓冲
		mutable std::shared_timed_mutex mutex; // 顶点级锁，用于并发插入/恢复

		bool recovery_flag; // 恢复过程中避免重复处理
		int32_t update_count; // 动态更新次数，用于 hot 顶点分类
		VertexClass vclass; // 可选 vertex-class to zone-class 策略使用
		double temperature; // ATZ-GC 指数衰减更新热度
		int32_t recent_updates; // ATZ-GC 最近一个窗口内的更新次数
		ATZVertexType atz_type; // ATZ-GC hot/warm/cold 分类

		Vertex(int32_t vid){
			id = vid;
			degree = 0;
			recovery_flag = false;
			update_count = 0;
			vclass = VERTEX_COLD;
			temperature = 0.0;
			recent_updates = 0;
			atz_type = ATZ_VERTEX_COLD;
			fill(cache_buf, cache_buf + CACHE_BLK_SIZE, -1);
		}

	~Vertex(){}
};

typedef struct block_information
{
	bool status = false;
	int32_t src = -1;
	int32_t blkIndex = -1; 
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
		blk_info_vec.resize(nb);
	}

	~ZoneStatus(){}
};

int nr_blocks_perzone = 0;
bool global_zns_record_enabled = true;
bool global_redo_log_enabled = false;
bool global_recovery_stream_redo = false;
uint64_t global_zns_record_seq = 1;
uint64_t global_current_epoch = 1;
bool global_zns_epoch_active = false;
uint64_t global_zns_epoch_first_seq = 0;
uint64_t global_zns_epoch_last_seq = 0;
uint32_t global_zns_epoch_data_count = 0;
uint32_t global_zns_epoch_aggregate_checksum = 0;
int32_t global_meta_zone = -1;
unsigned long long global_meta_wp = 0;
int32_t global_redo_zone = -1;
unsigned long long global_redo_wp = 0;
vector<int32_t> global_redo_zones;
int32_t global_current_redo_zone_idx = 0;
uint64_t global_next_op_id = 1;
uint64_t global_last_applied_op_id = 0;
uint64_t global_committed_op_id = 0;
bool global_edge_log_replay_active = false;
size_t global_zns_alloc_unit_size = ZNS_RECORD_SIZE;
std::mutex global_zns_record_mutex;
std::mutex global_redo_log_mutex;
std::mutex global_zns_device_write_mutex;
std::mutex global_gc_mutex;
std::mutex global_pending_merge_mutex;
std::shared_timed_mutex global_query_update_mutex;
vector<EdgeLogEntry> global_redo_pending_entries;
deque<vector<EdgeLogEntry>> global_redo_write_queue;
condition_variable global_redo_writer_cv;
condition_variable global_redo_writer_drained_cv;
condition_variable global_redo_writer_space_cv;
thread global_redo_writer_thread;
bool global_redo_writer_running = false;
bool global_redo_writer_stop = false;
uint32_t global_redo_writer_active_batches = 0;
size_t global_redo_writer_queued_entries = 0;
const size_t global_redo_flush_entry_threshold = 262144;
const size_t global_redo_max_queued_entries = 64 * global_redo_flush_entry_threshold;

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
			wp += zns_alloc_unit_size();
		}
    }

    // 这里不负责释放设备空间，保留空实现
    void clear() {

    }
};


typedef struct zone_information
{
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
int32_t cache_blk_threshold = 200000;
int32_t query_prefetch_block_limit = 0;
atomic<int64_t> global_pending_edge_count(0);
int64_t pending_merge_threshold = 4000000;
bool global_debug_pending_merge = false;
int32_t global_cache_buf_threshold = CACHE_BLK_SIZE;
int64_t global_edge_num = 0;
int32_t global_max_vid = 0;
vector<Vertex *> global_vectex_vec;
vector<new_edge> global_edge_vec;
vector<ZoneStatus *> global_zone_vec;
vector<unsigned long long> global_zone_write_wp;
vector<unsigned long long> global_zone_write_limit;
int32_t global_empty_zone = 0;

int32_t global_num_queues = 0;
int64_t global_batch_size = 1000000;  
double total_dump_time = 0.0;
int64_t total_dump_count = 0;
double total_gc_time = 0.0;
int64_t total_gc_count = 0;
int64_t total_gc_moved_blocks = 0;
double last_graph_maintenance_total_time = 0.0;
double last_graph_maintenance_compute_time = 0.0;
bool global_delete_compaction_enabled = false;
bool global_zone_class_enabled = false;
bool global_atz_gc_enabled = false;
bool global_basic_gc_enabled = false;
bool atz_temperature_initialized = false;
int32_t global_max_data_zones_override = 0;
int64_t global_gc_proactive_free_blocks = 0;
int32_t hot_update_threshold = 8;
int32_t hub_degree_threshold = 1024;
int32_t hub_block_threshold = 16;
vector<int32_t> hot_zone_ids;
vector<int32_t> cold_zone_ids;
vector<int32_t> hub_zone_ids;
vector<ATZZoneType> global_atz_zone_types;
double atz_temperature_alpha = 0.5;
double atz_hot_threshold = 8.0;
double atz_cold_threshold = 2.0;
double atz_cold_live_beta = 8.0;
double atz_cold_gc_threshold = 0.90;
int64_t atz_free_zone_emergency_watermark = 0;
int64_t atz_hot_zone_gc_count = 0;
int64_t atz_warm_zone_gc_count = 0;
int64_t atz_mixed_zone_gc_count = 0;
int64_t atz_cold_zone_gc_count = 0;
int64_t atz_emergency_gc_count = 0;
int64_t atz_migrated_hot_blocks = 0;
int64_t atz_migrated_warm_blocks = 0;
int64_t atz_migrated_cold_blocks = 0;
int64_t atz_cold_zone_skipped_count = 0;
int64_t atz_cold_zone_protected_count = 0;
int64_t atz_invalid_blocks_reclaimed = 0;
int64_t atz_gc_write_bytes = 0;
int64_t atz_user_write_bytes = 0;
int64_t atz_hot_block_update_writes = 0;
int64_t atz_warm_block_update_writes = 0;
int64_t atz_cold_block_update_writes = 0;

static inline size_t zns_alloc_unit_size(void)
{
	return global_zns_alloc_unit_size;
}

void initial() {
	/*
		预留的全局初始化入口，当前未使用。
	*/
}

#endif
