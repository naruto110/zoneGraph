#ifndef __UTILS_H_
#define __UTILS_H_

#include "global.h"	
// #include "pma.h"

// #define NVME_DEVICE_PATH "/dev/nullb0"
#define ZONE_GC_THRESHOLD 0.8

using namespace std;

bool test_block(int zone_idx, int blk_idx);
void readBlockFromZone(int fd, int32_t* &array, size_t size, unsigned long long wp);
void readBlockFromZoneInto(int fd, int32_t* array, size_t size, unsigned long long wp);
void writeBlockToZone(int fd, int32_t* array, size_t size, unsigned long long wp);
void writeZNSRecord(int fd, int32_t src_vid, int32_t pma_blk_idx, PMA *pma,
                    int32_t *array, size_t size, unsigned long long wp);
void beginZNSEpoch(void);
void commitZNSEpoch(void);
uint64_t appendZNSRedoLog(uint32_t op_type, int32_t src, int32_t dst);
bool readZNSRecordHeader(int fd, ZNSBlockHeader &header, unsigned long long wp);
bool readZNSRecordChecked(int fd, ZNSBlockHeader &header, int32_t* array,
                          size_t size, unsigned long long wp);
bool readZNSRedoLogChecked(int fd, EdgeLogRecord &record, unsigned long long wp);
void readZNSRecordPayload(int fd, int32_t* &array, size_t size, unsigned long long wp);
void writeGraphBlockToZone(Vertex *v, int32_t blk_idx, unsigned long long wp);
void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp);
void reset_zns_write_pointers(void);

// ZNS 设备的运行时信息。
// 这里既保存 libzbd 返回的设备元信息，也保存本工程用于块划分的参数。
struct zns
{
    /* 设备信息 */
    const char			*path; 
	int			block_size;
	int			zone_size;
	int			dev_fd;
	struct zbd_info		info;  // libzbd 返回的设备信息
	unsigned int		nr_zones;
	unsigned int		nr_conv_zones;
	struct zbd_zone		*zones;  // 设备返回的每个 zone 的状态信息
	// struct zone_information		*zones_info;
	int nr_blocks_perzone;
};

struct zns zns;

static void zns_fix_zone_values(struct zbd_zone *zbdz, int nrz)
{
	int i;

	if (zns.block_size == 1)
		return;

	for (i = 0; i < nrz; i++) {
		zbdz->start /= zns.block_size;
		zbdz->len /= zns.block_size;
		zbdz->capacity /= zns.block_size;
		if (!zbd_zone_cnv(zbdz))
			zbdz->wp /= zns.block_size;
		zbdz++;
	}
}

// 关闭设备
static void zns_close(void)
{

	if (zns.dev_fd < 0)
		return;

	zbd_close(zns.dev_fd);
	zns.dev_fd = -1;

	free(zns.zones);
}

// 打开设备，并初始化与 zone/block 相关的全局元数据
int zns_open(void) {
	unsigned int i;
	int ret;
	unsigned long long min_zone_capacity = 0;
	int total_blks = 0;

	/* 打开设备文件 */
	zns.dev_fd = zbd_open(zns.path, O_RDWR, &zns.info);
	if (zns.dev_fd < 0)
		return -1;

	// if (zns.block_size > 1 &&
	//     zns.info.zone_size % zns.block_size) {
	// 	fprintf(stderr, "Invalid block size\n");
	// 	return 1;
	// }

	/* 获取设备上的全部 zone 信息 */
	ret = zbd_list_zones(zns.dev_fd, 0, 0, ZBD_RO_ALL,
			     &zns.zones, &zns.nr_zones);
	if (ret != 0)
		goto out;
	if (!zns.nr_zones) {
		fprintf(stderr, "No zones reported\n");
		ret = 1;
		goto out;
	}

	zns.zone_size = zns.info.zone_size;
	{
		size_t payload_unit = global_zns_record_enabled ?
			static_cast<size_t>(ZNS_RECORD_SIZE) :
			static_cast<size_t>(ZONE_BLK_SIZE);
		size_t lblock = zns.info.lblock_size == 0 ? 1 :
			static_cast<size_t>(zns.info.lblock_size);
		global_zns_alloc_unit_size =
			((payload_unit + lblock - 1) / lblock) * lblock;
	}

	min_zone_capacity = static_cast<unsigned long long>(zns.zone_size);
	for (int32_t zidx = 0; zidx < static_cast<int32_t>(zns.nr_zones); zidx++) {
		if (zns.zones[zidx].capacity > 0) {
			min_zone_capacity = min(min_zone_capacity,
			                        static_cast<unsigned long long>(zns.zones[zidx].capacity));
		}
	}
	zns.nr_blocks_perzone = min_zone_capacity / zns_alloc_unit_size();
	if (zns.nr_zones == 0 || zns.nr_blocks_perzone == 0) {
		fprintf(stderr, "Invalid ZNS geometry: nr_zones=%u zone_size=%d block_size=%zu\n",
				zns.nr_zones, zns.zone_size, zns_alloc_unit_size());
		zns_close();
		return -1;
	}
	nr_blocks_perzone = zns.nr_blocks_perzone;
	if (global_zns_record_enabled) {
		if (zns.nr_zones < 4) {
			fprintf(stderr, "ZNS record mode requires at least 4 zones for data, GC, redo log, and metadata\n");
			zns_close();
			return -1;
		}
		global_meta_zone = zns.nr_zones - 1; // 最后一个 zone 专门追加 commit record
		global_redo_zone = zns.nr_zones - 2; // 倒数第二个 zone 专门追加 redo log
		global_empty_zone = zns.nr_zones - 3;  // 倒数第三个 zone 作为 GC 空白缓冲区
		global_num_queues = zns.nr_zones - 3; // 其余 zone 用于 data record 分配
		global_meta_wp = static_cast<unsigned long long>(global_meta_zone) *
		                 static_cast<unsigned long long>(zns.zone_size);
		global_redo_wp = static_cast<unsigned long long>(global_redo_zone) *
		                 static_cast<unsigned long long>(zns.zone_size);
	} else {
		global_meta_zone = -1;
		global_meta_wp = 0;
		global_redo_zone = -1;
		global_redo_wp = 0;
		global_empty_zone = zns.nr_zones - 1;  // 默认把最后一个 zone 作为空白缓冲区
		global_num_queues = zns.nr_zones - 1; // 其余 zone 用于分配可写块
	}

	total_blks = zns.nr_blocks_perzone;

	global_zone_vec.clear();
	for (int i = 0; i < static_cast<int>(zns.nr_zones); i++) {
		ZoneStatus *z = new ZoneStatus(total_blks, 0);
		global_zone_vec.push_back(z);
	}
	reset_zns_write_pointers();

	// zns_fix_zone_values(zns.zones, zns.nr_zones);  

	// for (i = 0; i < zns.nr_zones; i++) {
	// 	if (zbd_zone_cnv(&zns.zones[i]))
	// 		zns.nr_conv_zones++;
	// }

out:
	if (ret)
		zns_close();

	return ret;
}

static unsigned long long zns_zone_start(int32_t zone_idx)
{
	return static_cast<unsigned long long>(zone_idx) *
	       static_cast<unsigned long long>(zns.zone_size);
}

static unsigned long long zns_zone_limit(int32_t zone_idx)
{
	unsigned long long start = zns_zone_start(zone_idx);
	if (zone_idx >= 0 && zone_idx < static_cast<int32_t>(zns.nr_zones) &&
	    zns.zones != nullptr && zns.zones[zone_idx].capacity > 0) {
		return start + static_cast<unsigned long long>(zns.zones[zone_idx].capacity);
	}
	return start + static_cast<unsigned long long>(zns.zone_size);
}

void reset_zns_write_pointers(void)
{
	global_zone_write_wp.assign(zns.nr_zones, 0);
	global_zone_write_limit.assign(zns.nr_zones, 0);
	for (int32_t zone_idx = 0; zone_idx < static_cast<int32_t>(zns.nr_zones); zone_idx++) {
		global_zone_write_wp[zone_idx] = zns_zone_start(zone_idx);
		global_zone_write_limit[zone_idx] = zns_zone_limit(zone_idx);
	}
}

static int32_t zns_zone_idx_from_wp(unsigned long long wp)
{
	if (zns.zone_size <= 0) {
		return -1;
	}
	return static_cast<int32_t>(wp / static_cast<unsigned long long>(zns.zone_size));
}

static void reset_zns_write_pointer_for_zone(int32_t zone_idx)
{
	if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_zone_write_wp.size())) {
		return;
	}
	global_zone_write_wp[zone_idx] = zns_zone_start(zone_idx);
	global_zone_write_limit[zone_idx] = zns_zone_limit(zone_idx);
}

static void checkedZoneWrite(int fd, const void *data, size_t bytes,
                             unsigned long long wp, const char *tag)
{
	if (fd == -1) {
		perror("ZNS device is not open");
		exit(EXIT_FAILURE);
	}
	if (bytes > zns_alloc_unit_size()) {
		std::cerr << tag << " write size " << bytes
		          << " exceeds allocation unit " << zns_alloc_unit_size() << std::endl;
		exit(EXIT_FAILURE);
	}
	if ((wp % zns_alloc_unit_size()) != 0) {
		std::cerr << tag << " unaligned wp=" << wp
		          << " unit=" << zns_alloc_unit_size() << std::endl;
		exit(EXIT_FAILURE);
	}

	int32_t zone_idx = zns_zone_idx_from_wp(wp);
	if (zone_idx < 0 || zone_idx >= static_cast<int32_t>(global_zone_write_wp.size())) {
		std::cerr << tag << " invalid zone for wp=" << wp << std::endl;
		exit(EXIT_FAILURE);
	}

	vector<char> write_buf(zns_alloc_unit_size(), 0);
	memcpy(write_buf.data(), data, bytes);

	lock_guard<std::mutex> lock(global_zns_device_write_mutex);
	if (wp != global_zone_write_wp[zone_idx]) {
		std::cerr << tag << " out-of-order ZNS write: zone=" << zone_idx
		          << " expected_wp=" << global_zone_write_wp[zone_idx]
		          << " actual_wp=" << wp << std::endl;
		exit(EXIT_FAILURE);
	}
	if (wp + zns_alloc_unit_size() > global_zone_write_limit[zone_idx]) {
		std::cerr << tag << " write exceeds zone capacity: zone=" << zone_idx
		          << " wp=" << wp
		          << " limit=" << global_zone_write_limit[zone_idx] << std::endl;
		exit(EXIT_FAILURE);
	}
	if (pwrite(fd, write_buf.data(), zns_alloc_unit_size(), wp) !=
	    static_cast<ssize_t>(zns_alloc_unit_size())) {
		perror(tag);
		exit(EXIT_FAILURE);
	}
	global_zone_write_wp[zone_idx] += zns_alloc_unit_size();
}

void report_device_info(struct zbd_info info) {
	cout << "lblock_size: " << info.lblock_size << endl;
	cout << "max_nr_active_zones: " << info.max_nr_active_zones << endl;
	cout << "max_nr_open_zones: " << info.max_nr_open_zones << endl;
	cout << "model: " << info.model << endl;
	cout << "nr_sectors: " << info.nr_sectors << endl;
	cout << "nr_zones: " << info.nr_zones << endl;
	cout << "zone_size: " << info.zone_size << endl;
	cout << "zone_sectors: " << info.zone_sectors << endl;
	cout << "nr_blocks_perzone: " << zns.nr_blocks_perzone << endl;
}

void report_zones_info(struct zbd_zone	*zones, int nr_zones) {
	for (int i = 0; i < nr_zones; i++) {
		cout << "*** Zone " << i << " info: " << endl;
		cout << "capacity: " << zones->capacity << " ";
		cout << "start: " << zones->start << " ";
		cout << "cond: " << zones->cond << " ";
		cout << "len: " << zones->len << " ";
		cout << "wp: " << zones->wp << " ";
		cout << "type: " << zones->type << endl;
		// cout << "" << zones->start << endl;
		zones++;
	}
}

// 基线模式垃圾回收：
// 把有效块搬迁到当前空白 zone，再 reset 原 zone。
bool garbage_collection() {
	/*
		1. 遍历各个 zone，找到有效块
		2. 将有效块复制到空白 zone
		3. 更新顶点元数据和 zone 元数据
		4. reset 原 zone，并把它设为新的空白 zone
	*/
	lock_guard<std::mutex> gc_lock(global_gc_mutex);
	cout << "garbage collection begin!" << endl;
	unsigned long long wp = 0;
	unsigned long long rp = 0;
	int32_t pre_empty_zone = global_empty_zone;

	int32_t* array = new int32_t[PM_BLK_SIZE];
	for (int zidx = 0; zidx < zns.nr_zones; zidx++) {
		if (zidx == global_empty_zone || zidx == pre_empty_zone ||
		    zidx == global_meta_zone || zidx == global_redo_zone) {
			continue;
		}
		wp = global_empty_zone * zns.zone_size;  // 空白 zone 起始位置
		rp = zns.zone_size * zidx;

		int tmp_used_blk = 0;
		beginZNSEpoch();
		// 将该 zone 的有效块复制到空白 zone 中
		for (int blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
			if (test_block(zidx, blkidx)) { // 当前为有效块
				// cout << "Collect zone: " << zidx << " block: " << blkidx << endl;
				// rp = zns.zone_size * zidx + ZONE_BLK_SIZE * blkidx;
				readBlockFromZone(zns.dev_fd, array, PM_BLK_SIZE, rp);

				// for (int i=0; i < PM_BLK_SIZE; i++) {
				// 	cout << array[i] << ",";
				// }
				// cout << endl;

				// 更新顶点记录的物理地址
				int32_t vid = global_zone_vec[zidx]->blk_info_vec[blkidx].src;
				int32_t ori_blk_idx = global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex;
				Vertex *v = global_vectex_vec[vid];
				unique_lock<std::shared_timed_mutex> vertex_lock(v->mutex, defer_lock);
				if (!vertex_lock.try_lock()) {
					cout << "GC abort: busy vertex " << vid << endl;
					commitZNSEpoch();
					delete[] array;
					return false;
				}
				blk_info cur_info = global_zone_vec[zidx]->blk_info_vec[blkidx];
				if (!cur_info.status || cur_info.src != vid ||
				    cur_info.blkIndex != ori_blk_idx) {
					continue;
				}
				if (global_zns_record_enabled) {
					writeZNSRecord(zns.dev_fd, vid, ori_blk_idx, v->pma, array, PM_BLK_SIZE, wp);
				} else {
					writeBlockToZone(zns.dev_fd, array, PM_BLK_SIZE, wp);
				}
				refresh_pm_blk_meta(v->blk_list[ori_blk_idx], array);
				v->blk_list[ori_blk_idx].global_wp = wp; // 新块已写入设备，因此对应 buf 应保持为空
				v->blk_list[ori_blk_idx].buf = nullptr;
				global_zone_vec[zidx]->blk_info_vec[blkidx].status = false; // 原块失效
				// 标记新位置的映射关系
				mark_zone(vid, ori_blk_idx, wp);
				wp += zns_alloc_unit_size(); 
				tmp_used_blk++;
			}
			rp += zns_alloc_unit_size();
		}
		commitZNSEpoch();
		// 把空白 zone 剩余位置重新加入空闲块队列
		for (int i = tmp_used_blk; i < zns.nr_blocks_perzone; i++) {
			pm_blk blk{};
			blk.global_wp = wp;
			blk.buf = nullptr;
			graph_blk_que.push(blk);
			wp += zns_alloc_unit_size(); 
		}
	// reset 原 zone
		zbd_reset_zones(zns.dev_fd, (zidx * zns.zone_size), zns.zone_size);
		reset_zns_write_pointer_for_zone(zidx);
		global_empty_zone = zidx;
	}
	delete[] array;
	return true;
}

bool garbage_collection_new() {
		/*
			并行模式的垃圾回收版本。
			回收后将空闲块重新放回对应 zone 的分配队列中。
		*/
	std::cerr << "garbage_collection_new is disabled: it swaps the empty zone "
	          << "with data zones and breaks per-zone write queues." << std::endl;
	return false;
	lock_guard<std::mutex> gc_lock(global_gc_mutex);
	double start_t = get_current_time();
	cout << "garbage collection begin! global_empty_zone:" << global_empty_zone << endl;
	unsigned long long wp = 0;
	unsigned long long rp = 0;
	int32_t pre_empty_zone = global_empty_zone;

	// 清空现有 zone 队列
	for (int i = 0; i < global_blk_queues.size(); i++) {
		while (!global_blk_queues[i]->memoryWPs.empty()) {
			global_blk_queues[i]->memoryWPs.pop();
		}
	}

	int32_t* array = new int32_t[PM_BLK_SIZE];
	for (int zidx = 0; zidx < zns.nr_zones; zidx++) {
		if (zidx == global_empty_zone || zidx == pre_empty_zone ||
		    zidx == global_meta_zone || zidx == global_redo_zone) {
			continue;
		}
		// cout << "global empty zone: " << global_empty_zone << " zns.zone_size: " << zns.zone_size << endl;
		wp = (unsigned long long)global_empty_zone * (unsigned long long)zns.zone_size;  // 空白 zone 起始位置
		rp = (unsigned long long)zns.zone_size * (unsigned long long)zidx;

		int tmp_used_blk = 0;
		beginZNSEpoch();
		// 将有效块搬迁到空白 zone 中
		for (int blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
			if (test_block(zidx, blkidx)) { // 当前为有效块
				// cout << "Collect zone: " << zidx << " block: " << blkidx << endl;
				// rp = zns.zone_size * zidx + ZONE_BLK_SIZE * blkidx;
				// cout << "wp: " << wp << endl;
				readBlockFromZone(zns.dev_fd, array, PM_BLK_SIZE, rp);

				// for (int i=0; i < PM_BLK_SIZE; i++) {
				// 	cout << array[i] << ",";
				// }
				// cout << endl;

				// 更新顶点记录的物理地址
				int32_t vid = global_zone_vec[zidx]->blk_info_vec[blkidx].src;
				int32_t ori_blk_idx = global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex;
				Vertex *v = global_vectex_vec[vid];
				unique_lock<std::shared_timed_mutex> vertex_lock(v->mutex, defer_lock);
				if (!vertex_lock.try_lock()) {
					cout << "GC abort: busy vertex " << vid << endl;
					commitZNSEpoch();
					delete[] array;
					return false;
				}
				blk_info cur_info = global_zone_vec[zidx]->blk_info_vec[blkidx];
				if (!cur_info.status || cur_info.src != vid ||
				    cur_info.blkIndex != ori_blk_idx) {
					continue;
				}
				if (global_zns_record_enabled) {
					writeZNSRecord(zns.dev_fd, vid, ori_blk_idx, v->pma, array, PM_BLK_SIZE, wp);
				} else {
					writeBlockToZone(zns.dev_fd, array, PM_BLK_SIZE, wp);
				}
				refresh_pm_blk_meta(v->blk_list[ori_blk_idx], array);
				v->blk_list[ori_blk_idx].global_wp = wp; // 新块已写入设备，因此对应 buf 应保持为空
				v->blk_list[ori_blk_idx].buf = nullptr;
				global_zone_vec[zidx]->blk_info_vec[blkidx].status = false; // 原块失效
				// 标记新位置的映射关系
				mark_zone(vid, ori_blk_idx, wp);
				wp += zns_alloc_unit_size(); 
				tmp_used_blk++;
			}
			rp += zns_alloc_unit_size();
		}
		commitZNSEpoch();
		// 将空白 zone 剩余块重新加入对应分配队列
		for (int i = tmp_used_blk; i < zns.nr_blocks_perzone; i++) {
			global_blk_queues[zidx]->memoryWPs.push(wp);
			// pm_blk blk;
			// blk.global_wp = wp;
			// blk.buf = nullptr;
			// graph_blk_que.push(blk);
			wp += zns_alloc_unit_size(); 
		}
		// reset 原 zone
		zbd_reset_zones(zns.dev_fd, (zidx * zns.zone_size), zns.zone_size);
		reset_zns_write_pointer_for_zone(zidx);
		global_empty_zone = zidx;
	}
	delete[] array;

	double end_t = get_current_time();
    cout << "Garbage collection time: " << end_t - start_t << endl;
	return true;
}


// 基线模式：从连续地址空间中预分配若干逻辑块
int generate_block_baseline(int blk_num) {
	/*
		从当前全局写指针开始切分出若干块。
		该模式不关心顶点之间是否落在同一个 zone。
	*/
	// cout << "generate zone blocks" << endl;
	while (blk_num--) {
		int32_t zone_idx = zns_zone_idx_from_wp(global_blk_wp);
		while (zone_idx >= 0 && zone_idx < global_num_queues &&
		       global_blk_wp + zns_alloc_unit_size() >
		       zns_zone_start(zone_idx) +
		       static_cast<unsigned long long>(zns.nr_blocks_perzone) *
		       zns_alloc_unit_size()) {
			global_blk_wp = zns_zone_start(zone_idx + 1);
			zone_idx = zns_zone_idx_from_wp(global_blk_wp);
		}
		if (zone_idx < 0 || zone_idx >= global_num_queues) {
			std::cerr << "No free baseline blocks before reserved empty zone." << std::endl;
			exit(EXIT_FAILURE);
		}
		pm_blk blk{};
		blk.global_wp = global_blk_wp; // 记录该块在设备上的起始地址
		// blk.blk_id = global_blk_id;
		blk.buf = nullptr;
		global_blk_id++;
		global_blk_wp += zns_alloc_unit_size(); // 以字节为单位推进写指针
		graph_blk_que.push(blk);
		// if (global_blk_wp >= global_threshold_wp) {

		// }
	}
	// cout << "generate zone blocks done, total blocks: " << global_blk_id << endl;

	return 1;
}

void generate_block_opt(void)
{
	// 初始化每个 zone 对应的可写块队列
	global_blk_queues.clear();
	hot_zone_ids.clear();
	cold_zone_ids.clear();
	hub_zone_ids.clear();

	unsigned long long wp = 0;
    for (int zoneidx = 0; zoneidx < global_num_queues; zoneidx++) {
		MemoryWPQueue* b = new MemoryWPQueue(wp);
        global_blk_queues.push_back(b);
		wp += zns.zone_size;
    }

	if (global_num_queues <= 0) {
		return;
	}

	int32_t hub_count = global_num_queues >= 3 ? 1 : 0;
	int32_t hot_count = max(1, global_num_queues / 4);
	if (hot_count + hub_count > global_num_queues) {
		hot_count = global_num_queues - hub_count;
	}
	int32_t hub_start = global_num_queues - hub_count;
	int32_t cold_start = hot_count;
	int32_t cold_end = hub_start;

	for (int32_t zoneidx = 0; zoneidx < hot_count; zoneidx++) {
		hot_zone_ids.push_back(zoneidx);
	}
	for (int32_t zoneidx = cold_start; zoneidx < cold_end; zoneidx++) {
		cold_zone_ids.push_back(zoneidx);
	}
	for (int32_t zoneidx = hub_start; zoneidx < global_num_queues; zoneidx++) {
		hub_zone_ids.push_back(zoneidx);
	}

	if (cold_zone_ids.empty()) {
		cold_zone_ids = hot_zone_ids;
	}
}

// 并行模式申请可写块
unsigned long long requestMemoryBlocks(int thread_id) 
{
    unsigned long long wp = 0;

	if (global_num_queues <= 0 || global_blk_queues.empty()) {
		std::cerr << "No writable zone queues are available." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (global_blk_queues[(thread_id % global_num_queues)]->getNextBlock(wp)) {
		return wp;
	}
    // 如果线程绑定的队列没有空闲块，则遍历所有队列
	int iters = 32;
	bool blockAcquired = false;
	while ((iters>0) && (!blockAcquired)) {
		for (int i = 0; i < global_num_queues; i++) {
			if (global_blk_queues[i]->getNextBlock(wp)) {  // 单个队列内部已加锁，分配是线程安全的
				blockAcquired = true;
				break; // 成功获取块后退出
			}
		}
		iters--;
	}

	// 检查是否成功获取可写块
	if (blockAcquired) {
		return wp;
	} else {
		std::cout << "All memory blocks are in use. Initiating garbage collection." << std::endl;
		exit(1);
	}
}

void writeBlockToZone(int fd, int32_t* array, size_t size, unsigned long long wp) {
    if (fd == -1) {
        perror("Failed to open file for writing");
        exit(EXIT_FAILURE);
    }

    // 按块粒度写入设备
    size_t bytesToWrite = size * sizeof(int32_t);
    checkedZoneWrite(fd, array, bytesToWrite, wp, "Failed to write data block");
}

void readBlockFromZone(int fd, int32_t* &array, size_t size, unsigned long long wp) {
	if (global_zns_record_enabled) {
		readZNSRecordPayload(fd, array, size, wp);
		return;
	}

    // 从设备上读取一个逻辑块
    if (fd == -1) {
        perror("Failed to open file for reading");
        exit(EXIT_FAILURE);
    }

    // 如有需要，先为读取结果分配内存
	if (array == nullptr) {
		array = new int32_t[size];
	}
    // 

    // 将设备数据读入内存
    size_t bytesToRead = size * sizeof(int32_t);
    if (pread(fd, array, bytesToRead, wp) != bytesToRead) {
        perror("Failed to read data");
        delete[] array;
        array = nullptr;
        exit(EXIT_FAILURE);
    }

    // close(fd);
}

void readBlockFromZoneInto(int fd, int32_t* array, size_t size, unsigned long long wp) {
	if (array == nullptr) {
		std::cerr << "readBlockFromZoneInto requires a caller-owned buffer." << std::endl;
		exit(EXIT_FAILURE);
	}
	if (global_zns_record_enabled) {
		int32_t *array_ref = array;
		readZNSRecordPayload(fd, array_ref, size, wp);
		return;
	}
	if (fd == -1) {
		perror("Failed to open file for reading");
		exit(EXIT_FAILURE);
	}
	size_t bytesToRead = size * sizeof(int32_t);
	if (pread(fd, array, bytesToRead, wp) != bytesToRead) {
		perror("Failed to read data");
		exit(EXIT_FAILURE);
	}
}

static uint32_t zns_record_checksum(const ZNSBlockHeader &header,
                                    const int32_t *array,
                                    size_t size)
{
    uint32_t hash = 2166136261u;
    ZNSBlockHeader tmp = header;
    tmp.checksum = 0;
    const unsigned char *header_bytes =
        reinterpret_cast<const unsigned char *>(&tmp);
    for (size_t i = 0; i < sizeof(ZNSBlockHeader); i++) {
        hash ^= header_bytes[i];
        hash *= 16777619u;
    }
    if (array != nullptr && size > 0) {
        const unsigned char *payload_bytes =
            reinterpret_cast<const unsigned char *>(array);
        for (size_t i = 0; i < size * sizeof(int32_t); i++) {
            hash ^= payload_bytes[i];
            hash *= 16777619u;
        }
    }
    return hash;
}

static uint32_t combine_epoch_checksum(uint32_t aggregate,
                                       uint64_t seq,
                                       uint32_t checksum)
{
    uint32_t mixed = checksum;
    mixed ^= static_cast<uint32_t>(seq) * 0x9e3779b9u;
    mixed ^= static_cast<uint32_t>(seq >> 32) * 0x85ebca6bu;
    mixed ^= mixed >> 16;
    mixed *= 0x7feb352du;
    mixed ^= mixed >> 15;
    return aggregate ^ mixed;
}

static uint32_t count_payload_values(const int32_t *array, size_t size)
{
    uint32_t count = 0;
    for (size_t i = 0; i < size; i++) {
        if (array[i] >= 0) {
            count++;
        }
    }
    return count;
}

static uint32_t zns_redo_checksum(const EdgeLogRecord &record)
{
    uint32_t hash = 2166136261u;
    EdgeLogRecord tmp = record;
    tmp.checksum = 0;
    const unsigned char *bytes =
        reinterpret_cast<const unsigned char *>(&tmp);
    for (size_t i = 0; i < sizeof(EdgeLogRecord); i++) {
        hash ^= bytes[i];
        hash *= 16777619u;
    }
    return hash;
}

uint64_t appendZNSRedoLog(uint32_t op_type, int32_t src, int32_t dst)
{
    if (!global_zns_record_enabled || global_edge_log_replay_active) {
        return 0;
    }
    if (zns.dev_fd == -1) {
        perror("Failed to open file for redo log writing");
        exit(EXIT_FAILURE);
    }
    if (global_redo_zone < 0) {
        std::cerr << "ZNS redo log requires a dedicated redo zone." << std::endl;
        exit(EXIT_FAILURE);
    }

    EdgeLogRecord record;
    unsigned long long wp = 0;
    {
        lock_guard<std::mutex> lock(global_redo_log_mutex);
        if (global_redo_wp + zns_alloc_unit_size() >
            global_zone_write_limit[global_redo_zone]) {
            std::cerr << "No free space left in ZNS redo log zone." << std::endl;
            exit(EXIT_FAILURE);
        }

        memset(&record, 0, sizeof(record));
        record.magic = ZNS_REDO_MAGIC;
        record.version = ZNS_BLOCK_VERSION;
        record.op_type = op_type;
        record.op_id = global_next_op_id++;
        record.src = src;
        record.dst = dst;
        record.tail_magic = ZNS_REDO_TAIL_MAGIC;
        record.checksum = zns_redo_checksum(record);

        wp = global_redo_wp;
        global_redo_wp += zns_alloc_unit_size();
        global_last_applied_op_id = record.op_id;

        char raw[ZNS_RECORD_SIZE];
        memset(raw, 0, sizeof(raw));
        memcpy(raw, &record, sizeof(record));
        checkedZoneWrite(zns.dev_fd, raw, ZNS_RECORD_SIZE, wp,
                         "Failed to write ZNS redo log record");
    }
    return record.op_id;
}

void writeZNSRecord(int fd, int32_t src_vid, int32_t pma_blk_idx, PMA *pma,
                    int32_t *array, size_t size, unsigned long long wp)
{
    if (fd == -1) {
        perror("Failed to open file for ZNS record writing");
        exit(EXIT_FAILURE);
    }
    if (sizeof(ZNSBlockHeader) + size * sizeof(int32_t) > ZNS_RECORD_SIZE) {
        std::cerr << "ZNS record payload is larger than ZNS_RECORD_SIZE." << std::endl;
        exit(EXIT_FAILURE);
    }

    char record[ZNS_RECORD_SIZE];
    memset(record, 0, sizeof(record));

    ZNSBlockHeader header;
    memset(&header, 0, sizeof(header));
    header.magic = ZNS_BLOCK_MAGIC;
    header.version = ZNS_BLOCK_VERSION;
    header.record_type = ZNS_RECORD_DATA;
    header.payload_bytes = static_cast<uint32_t>(size * sizeof(int32_t));
    header.src_vid = src_vid;
    header.pma_blk_idx = pma_blk_idx;
    header.valid_count = count_payload_values(array, size);
    header.pma_ele = pma == nullptr ? header.valid_count : static_cast<uint32_t>(pma->ele_num);
    header.pma_cap = pma == nullptr ? static_cast<uint32_t>(size) : static_cast<uint32_t>(pma->arr_cap);
    header.tail_magic = ZNS_RECORD_TAIL_MAGIC;

    {
        lock_guard<std::mutex> lock(global_zns_record_mutex);
        if (!global_zns_epoch_active) {
            global_zns_epoch_active = true;
            global_zns_epoch_first_seq = 0;
            global_zns_epoch_last_seq = 0;
            global_zns_epoch_data_count = 0;
            global_zns_epoch_aggregate_checksum = 0;
        }
        header.epoch = global_current_epoch;
        header.seq = global_zns_record_seq++;
        header.checksum = zns_record_checksum(header, array, size);

        if (global_zns_epoch_data_count == 0) {
            global_zns_epoch_first_seq = header.seq;
        }
        global_zns_epoch_last_seq = header.seq;
        global_zns_epoch_data_count++;
        global_zns_epoch_aggregate_checksum =
            combine_epoch_checksum(global_zns_epoch_aggregate_checksum,
                                   header.seq,
                                   header.checksum);

        memcpy(record, &header, sizeof(header));
        memcpy(record + sizeof(ZNSBlockHeader), array, size * sizeof(int32_t));
        checkedZoneWrite(fd, record, ZNS_RECORD_SIZE, wp,
                         "Failed to write ZNS record");
    }
}

void beginZNSEpoch(void)
{
    if (!global_zns_record_enabled) {
        return;
    }
    lock_guard<std::mutex> lock(global_zns_record_mutex);
    if (global_zns_epoch_active) {
        return;
    }
    global_zns_epoch_active = true;
    global_zns_epoch_first_seq = 0;
    global_zns_epoch_last_seq = 0;
    global_zns_epoch_data_count = 0;
    global_zns_epoch_aggregate_checksum = 0;
}

void commitZNSEpoch(void)
{
    if (!global_zns_record_enabled) {
        return;
    }

    ZNSBlockHeader header;
    memset(&header, 0, sizeof(header));
    unsigned long long commit_wp = 0;
    {
        lock_guard<std::mutex> lock(global_zns_record_mutex);
        if (!global_zns_epoch_active || global_zns_epoch_data_count == 0) {
            global_zns_epoch_active = false;
            return;
        }
        if (global_meta_zone < 0 ||
            global_meta_wp + zns_alloc_unit_size() >
            global_zone_write_limit[global_meta_zone]) {
            std::cerr << "No free space left in ZNS metadata zone for commit record." << std::endl;
            exit(EXIT_FAILURE);
        }

        header.magic = ZNS_BLOCK_MAGIC;
        header.version = ZNS_BLOCK_VERSION;
        header.record_type = ZNS_RECORD_COMMIT;
        header.payload_bytes = 0;
        header.epoch = global_current_epoch;
        header.src_vid = -1;
        header.pma_blk_idx = -1;
        header.seq = global_zns_record_seq++;
        header.first_seq = global_zns_epoch_first_seq;
        header.last_seq = global_zns_epoch_last_seq;
        header.data_record_count = global_zns_epoch_data_count;
        header.aggregate_checksum = global_zns_epoch_aggregate_checksum;
        header.committed_op_id = global_last_applied_op_id;
        header.tail_magic = ZNS_RECORD_TAIL_MAGIC;
        header.checksum = zns_record_checksum(header, nullptr, 0);

        commit_wp = global_meta_wp;
        global_meta_wp += zns_alloc_unit_size();
        global_current_epoch++;
        global_zns_epoch_active = false;
        global_zns_epoch_first_seq = 0;
        global_zns_epoch_last_seq = 0;
        global_zns_epoch_data_count = 0;
        global_zns_epoch_aggregate_checksum = 0;

        char record[ZNS_RECORD_SIZE];
        memset(record, 0, sizeof(record));
        memcpy(record, &header, sizeof(header));
        checkedZoneWrite(zns.dev_fd, record, ZNS_RECORD_SIZE, commit_wp,
                         "Failed to write ZNS commit record");
    }
}

bool readZNSRecordHeader(int fd, ZNSBlockHeader &header, unsigned long long wp)
{
    if (fd == -1) {
        perror("Failed to open file for ZNS record header reading");
        exit(EXIT_FAILURE);
    }
    char record[ZNS_RECORD_SIZE];
    if (pread(fd, record, ZNS_RECORD_SIZE, wp) != ZNS_RECORD_SIZE) {
        perror("Failed to read ZNS record header");
        exit(EXIT_FAILURE);
    }
    memcpy(&header, record, sizeof(ZNSBlockHeader));
    return header.magic == ZNS_BLOCK_MAGIC &&
           header.version == ZNS_BLOCK_VERSION &&
           header.tail_magic == ZNS_RECORD_TAIL_MAGIC;
}

bool readZNSRecordChecked(int fd, ZNSBlockHeader &header, int32_t* array,
                          size_t size, unsigned long long wp)
{
    if (fd == -1) {
        perror("Failed to open file for ZNS record reading");
        exit(EXIT_FAILURE);
    }
    char record[ZNS_RECORD_SIZE];
    if (pread(fd, record, ZNS_RECORD_SIZE, wp) != ZNS_RECORD_SIZE) {
        perror("Failed to read ZNS record");
        exit(EXIT_FAILURE);
    }
    memcpy(&header, record, sizeof(ZNSBlockHeader));
    if (header.magic != ZNS_BLOCK_MAGIC ||
        header.version != ZNS_BLOCK_VERSION ||
        header.tail_magic != ZNS_RECORD_TAIL_MAGIC) {
        return false;
    }
    if (header.record_type == ZNS_RECORD_DATA) {
        if (header.payload_bytes != size * sizeof(int32_t) || array == nullptr) {
            return false;
        }
        memcpy(array, record + sizeof(ZNSBlockHeader), size * sizeof(int32_t));
        return zns_record_checksum(header, array, size) == header.checksum;
    }
    if (header.record_type == ZNS_RECORD_COMMIT) {
        return zns_record_checksum(header, nullptr, 0) == header.checksum;
    }
    return false;
}

bool readZNSRedoLogChecked(int fd, EdgeLogRecord &record, unsigned long long wp)
{
    if (fd == -1) {
        perror("Failed to open file for redo log reading");
        exit(EXIT_FAILURE);
    }
    char raw[ZNS_RECORD_SIZE];
    if (pread(fd, raw, ZNS_RECORD_SIZE, wp) != ZNS_RECORD_SIZE) {
        perror("Failed to read ZNS redo log record");
        exit(EXIT_FAILURE);
    }
    memcpy(&record, raw, sizeof(EdgeLogRecord));
    if (record.magic != ZNS_REDO_MAGIC ||
        record.version != ZNS_BLOCK_VERSION ||
        record.tail_magic != ZNS_REDO_TAIL_MAGIC ||
        (record.op_type != EDGE_LOG_INSERT &&
         record.op_type != EDGE_LOG_DELETE)) {
        return false;
    }
    return zns_redo_checksum(record) == record.checksum;
}

void readZNSRecordPayload(int fd, int32_t* &array, size_t size, unsigned long long wp)
{
    if (fd == -1) {
        perror("Failed to open file for ZNS record payload reading");
        exit(EXIT_FAILURE);
    }
    char record[ZNS_RECORD_SIZE];
    if (pread(fd, record, ZNS_RECORD_SIZE, wp) != ZNS_RECORD_SIZE) {
        perror("Failed to read ZNS record payload");
        exit(EXIT_FAILURE);
    }
    ZNSBlockHeader header;
    memcpy(&header, record, sizeof(ZNSBlockHeader));
    if (header.magic != ZNS_BLOCK_MAGIC ||
        header.version != ZNS_BLOCK_VERSION ||
        header.record_type != ZNS_RECORD_DATA ||
        header.payload_bytes != size * sizeof(int32_t) ||
        header.tail_magic != ZNS_RECORD_TAIL_MAGIC) {
        std::cerr << "Invalid ZNS record header at wp=" << wp << std::endl;
        exit(EXIT_FAILURE);
    }
	if (array == nullptr) {
		array = new int32_t[size];
	}
    memcpy(array, record + sizeof(ZNSBlockHeader), size * sizeof(int32_t));
    uint32_t checksum = zns_record_checksum(header, array, size);
    if (checksum != header.checksum) {
        std::cerr << "ZNS record checksum mismatch at wp=" << wp << std::endl;
        exit(EXIT_FAILURE);
    }
}

void writeGraphBlockToZone(Vertex *v, int32_t blk_idx, unsigned long long wp)
{
    if (global_zns_record_enabled) {
        writeZNSRecord(zns.dev_fd, v->id, blk_idx, v->pma,
                       v->blk_list[blk_idx].buf, PM_BLK_SIZE, wp);
    } else {
        writeBlockToZone(zns.dev_fd, v->blk_list[blk_idx].buf, PM_BLK_SIZE, wp);
    }
}

void validBlock(unsigned long long	wp) {
	int32_t zone_idx = wp / zns.zone_size; 
	int32_t blk_idx = (wp % zns.zone_size) / zns_alloc_unit_size();
	global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status = true;
	// global_zone_vec[zone_idx]
	// zns.zones_info[zone_idx].blk_info_vec[blk_idx].status = 1;
}

void invalidBlock(unsigned long long wp) {
	int32_t zone_idx = wp / zns.zone_size; 
	int32_t blk_idx = (wp % zns.zone_size) / zns_alloc_unit_size();
	// cout << "wp: " << wp << " zone size: " << zns.zone_size << 
	
	global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status = false;
	// cout << "Invalid zone idx: " << zone_idx << " blk idx: " << blk_idx 
	// << " status: " << global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status << endl;
}

bool test_block(int zone_idx, int blk_idx) 
{
	return global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status;
}

// 记录物理块到顶点逻辑块的映射，供 GC 和恢复使用
void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp)
{
    // cout << "zone size: " << zns.info.zone_size << endl;
    int32_t zone_idx = wp / zns.info.zone_size; // 所属 zone 编号
    int32_t info_idx = (wp % zns.info.zone_size) / zns_alloc_unit_size(); // 在 zone 内的块编号
    global_zone_vec[zone_idx]->blk_info_vec[info_idx].status = true;
    global_zone_vec[zone_idx]->blk_info_vec[info_idx].src = vid;
    global_zone_vec[zone_idx]->blk_info_vec[info_idx].blkIndex = blk_idx;
	// cout << "mark zone-- zone: " << zone_idx << " blk: " << 
}

// test_block(int zone_idx, int blk_idx) 


void cacheBlock(Vertex *v, int32_t blk_idx)
{
	v->blk_list[blk_idx].buf = new int32_t[PM_BLK_SIZE]; // 为块缓存分配 DRAM
	readBlockFromZone(zns.dev_fd, v->blk_list[blk_idx].buf, PM_BLK_SIZE, v->blk_list[blk_idx].global_wp); // 从设备读入该块
	refresh_pm_blk_meta(v->blk_list[blk_idx], v->blk_list[blk_idx].buf);
	invalidBlock(v->blk_list[blk_idx].global_wp); // 读回后旧物理位置不再视作当前有效版本
	global_cache_blk++;
}

void releaseBlock(Vertex *v, int32_t blk_idx)
{
	delete[] v->blk_list[blk_idx].buf;
	v->blk_list[blk_idx].buf = nullptr;
	global_cache_blk--;
}

#endif
