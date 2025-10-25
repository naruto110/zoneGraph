#ifndef __UTILS_H_
#define __UTILS_H_

#include "global.h"	
// #include "pma.h"

// #define NVME_DEVICE_PATH "/dev/nullb0"
#define ZONE_GC_THRESHOLD 0.8

using namespace std;

bool test_block(int zone_idx, int blk_idx);
void readBlockFromZone(int fd, int32_t* &array, size_t size, unsigned long long wp);
void writeBlockToZone(int fd, int32_t* array, size_t size, unsigned long long wp);
void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp);

struct zns
{
    /* Device inforamtion */
    const char			*path; 
	int			block_size;
	int			zone_size;
	int			dev_fd;
	struct zbd_info		info;  // zone device information
	unsigned int		nr_zones;
	unsigned int		nr_conv_zones;
	struct zbd_zone		*zones;  // per zone status
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

// close a device
static void zns_close(void)
{

	if (zns.dev_fd < 0)
		return;

	zbd_close(zns.dev_fd);
	zns.dev_fd = -1;

	free(zns.zones);
}

// open a device
int zns_open(void) {
	unsigned int i;
	int ret;

	/* Open device file */
	zns.dev_fd = zbd_open(zns.path, O_RDWR, &zns.info);
	if (zns.dev_fd < 0)
		return -1;

	// if (zns.block_size > 1 &&
	//     zns.info.zone_size % zns.block_size) {
	// 	fprintf(stderr, "Invalid block size\n");
	// 	return 1;
	// }

	zns.nr_zones = zns.info.nr_zones;
	zns.zone_size = zns.info.zone_size;
	zns.nr_blocks_perzone = zns.zone_size / ZONE_BLK_SIZE;
	nr_blocks_perzone = zns.nr_blocks_perzone;
	global_empty_zone = zns.info.nr_zones - 1;  // 最后一个zone位缓冲zone
	global_num_queues = zns.nr_zones - 1; // 留下一个作为buffer

	int total_blks = zns.info.zone_size / ZONE_BLK_SIZE;

	for (int i = 0; i < zns.nr_zones; i++) {
		ZoneStatus *z = new ZoneStatus(total_blks, 0);
		z->blk_info_vec.reserve(total_blks);  // 默认都是无效的
		// cout << "test---0 blk_info_vec[5].status: " << z->blk_info_vec[5].status << endl;
		global_zone_vec.push_back(z);
	}

	/* Get list of all zones */
	ret = zbd_list_zones(zns.dev_fd, 0, 0, ZBD_RO_ALL,
			     &zns.zones, &zns.nr_zones);
	if (ret != 0)
		goto out;
	if (!zns.nr_zones) {
		fprintf(stderr, "No zones reported\n");
		ret = 1;
		goto out;
	}

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

// garbage collection
void garbage_collection() {
	/*
		1. 遍历所有zone，查看有效blk数量
		2. 存在无效blk则进行回收
		3. 先将有效数据迁移至空zone中，再reset原先zone
			a. 需要修改zone bitmap
			b. 需要修改vertex metadata
	*/
	cout << "garbage collection begin!" << endl;
	unsigned long long wp = 0;
	unsigned long long rp = 0;
	int32_t pre_empty_zone = global_empty_zone;

	int32_t* array = new int32_t[PM_BLK_SIZE];
	for (int zidx = 0; zidx < zns.nr_zones; zidx++) {
		if (zidx == global_empty_zone || zidx == pre_empty_zone) {
			continue;
		}
		wp = global_empty_zone * zns.zone_size;  // 缓冲zone起始位置
		rp = zns.zone_size * zidx;

		int tmp_used_blk = 0;
		// 将数据存入empty zone中，并将空余blk存入graph_blk_que
		for (int blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
			if (test_block(zidx, blkidx)) { // 当前为valid blk
				// cout << "Collect zone: " << zidx << " block: " << blkidx << endl;
				// rp = zns.zone_size * zidx + ZONE_BLK_SIZE * blkidx;
				readBlockFromZone(zns.dev_fd, array, PM_BLK_SIZE, rp);
				writeBlockToZone(zns.dev_fd, array, PM_BLK_SIZE, wp);

				// for (int i=0; i < PM_BLK_SIZE; i++) {
				// 	cout << array[i] << ",";
				// }
				// cout << endl;

				// 更行vertex metadata
				int32_t vid = global_zone_vec[zidx]->blk_info_vec[blkidx].src;
				int32_t ori_blk_idx = global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex;
				Vertex *v = global_vectex_vec[vid];
				v->blk_list[ori_blk_idx].global_wp = wp; // 由于该blk已经被写入zone中，其metadata中的buf为nullptr
				global_zone_vec[zidx]->blk_info_vec[blkidx].status = false; // 无效化原先blk
				// 修改更新后zone的标签：vid，blkidx
				mark_zone(vid, ori_blk_idx, wp);
				wp += ZONE_BLK_SIZE; 
				tmp_used_blk++;
			}
			rp += ZONE_BLK_SIZE;
		}
		// 将缓冲zone剩余部分加入graph_blk_que中
		for (int i = tmp_used_blk; i < zns.nr_blocks_perzone; i++) {
			pm_blk blk;
			blk.global_wp = wp;
			blk.buf = nullptr;
			graph_blk_que.push(blk);
			wp += ZONE_BLK_SIZE; 
		}
		// 将该zone reset
		zbd_reset_zones(zns.dev_fd, (zidx * zns.zone_size), zns.zone_size);
		global_empty_zone = zidx;
	}
	delete[] array;
}

void garbage_collection_new() {
	/*
		针对不同zone空间独立分配的特点
	*/
	double start_t = get_current_time();
	cout << "garbage collection begin! global_empty_zone:" << global_empty_zone << endl;
	unsigned long long wp = 0;
	unsigned long long rp = 0;
	int32_t pre_empty_zone = global_empty_zone;

	// 清空现有zone queue
	for (int i = 0; i < global_blk_queues.size(); i++) {
		while (!global_blk_queues[i]->memoryWPs.empty()) {
			global_blk_queues[i]->memoryWPs.pop();
		}
	}

	int32_t* array = new int32_t[PM_BLK_SIZE];
	for (int zidx = 0; zidx < zns.nr_zones; zidx++) {
		if (zidx == global_empty_zone || zidx == pre_empty_zone) {
			continue;
		}
		// cout << "global empty zone: " << global_empty_zone << " zns.zone_size: " << zns.zone_size << endl;
		wp = (unsigned long long)global_empty_zone * (unsigned long long)zns.zone_size;  // 缓冲zone起始位置
		rp = (unsigned long long)zns.zone_size * (unsigned long long)zidx;

		int tmp_used_blk = 0;
		// 将数据存入empty zone中，并将空余blk存入graph_blk_que
		for (int blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
			if (test_block(zidx, blkidx)) { // 当前为valid blk
				// cout << "Collect zone: " << zidx << " block: " << blkidx << endl;
				// rp = zns.zone_size * zidx + ZONE_BLK_SIZE * blkidx;
				// cout << "wp: " << wp << endl;
				readBlockFromZone(zns.dev_fd, array, PM_BLK_SIZE, rp);
				writeBlockToZone(zns.dev_fd, array, PM_BLK_SIZE, wp);

				// for (int i=0; i < PM_BLK_SIZE; i++) {
				// 	cout << array[i] << ",";
				// }
				// cout << endl;

				// 更行vertex metadata
				int32_t vid = global_zone_vec[zidx]->blk_info_vec[blkidx].src;
				int32_t ori_blk_idx = global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex;
				Vertex *v = global_vectex_vec[vid];
				v->blk_list[ori_blk_idx].global_wp = wp; // 由于该blk已经被写入zone中，其metadata中的buf为nullptr
				global_zone_vec[zidx]->blk_info_vec[blkidx].status = false; // 无效化原先blk
				// 修改更新后zone的标签：vid，blkidx
				mark_zone(vid, ori_blk_idx, wp);
				wp += ZONE_BLK_SIZE; 
				tmp_used_blk++;
			}
			rp += ZONE_BLK_SIZE;
		}
		// 将缓冲zone剩余部分加入graph_blk_que中
		for (int i = tmp_used_blk; i < zns.nr_blocks_perzone; i++) {
			global_blk_queues[zidx]->memoryWPs.push(wp);
			// pm_blk blk;
			// blk.global_wp = wp;
			// blk.buf = nullptr;
			// graph_blk_que.push(blk);
			wp += ZONE_BLK_SIZE; 
		}
		// 将该zone reset
		zbd_reset_zones(zns.dev_fd, (zidx * zns.zone_size), zns.zone_size);
		global_empty_zone = zidx;
	}
	delete[] array;

	double end_t = get_current_time();
    cout << "Garbage collection time: " << end_t - start_t << endl;
}


// allocate block
int generate_block_baseline(int blk_num) {
	/*
		从现有zone中申请blocks, 不考虑是否同一zone，邻节点随机分布
	*/
	// cout << "generate zone blocks" << endl;
	while (blk_num--) {
		pm_blk blk;
		blk.global_wp = global_blk_wp; // 根据global_blk_wp可以计算出zone id和局部wp
		// blk.blk_id = global_blk_id;
		blk.buf = nullptr;
		global_blk_id++;
		global_blk_wp += ZONE_BLK_SIZE; // 偏移的是字节数
		graph_blk_que.push(blk);
		// if (global_blk_wp >= global_threshold_wp) { // 需要进行垃圾回收

		// }
	}
	// cout << "generate zone blocks done, total blocks: " << global_blk_id << endl;

	return 1;
}

void generate_block_opt(void)
{
	// 初始化内存分配
	unsigned long long wp = 0;
    for (int zoneidx = 0; zoneidx < global_num_queues; zoneidx++) {
		MemoryWPQueue* b = new MemoryWPQueue(wp);
        global_blk_queues.push_back(b);
		wp += zns.zone_size;
    }
}

// 并发申请内存块
unsigned long long requestMemoryBlocks(int thread_id) 
{
    unsigned long long wp = 0;

	if (global_blk_queues[(thread_id % global_num_queues)]->getNextBlock(wp)) {
		return wp;
	}
    // 尝试从每个队列获取内存块
	int iters = 32;
	bool blockAcquired = false;
	while ((iters>0) && (!blockAcquired)) {
		for (int i = 0; i < global_num_queues; i++) {
			if (global_blk_queues[i]->getNextBlock(wp)) {  // 此时已经申请锁了，该对象（queue）会被lock，在此函数内完成对zone的写操作是串行的
				blockAcquired = true;
				break; // 成功获取内存块，退出循环
			}
		}
		iters--;
	}

	// 检查是否获得了内存块
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
        return;
    }

    // 将数组数据写入文件
    size_t bytesToWrite = size * sizeof(int32_t);
    if (pwrite(fd, array, bytesToWrite, wp) != bytesToWrite) {
        perror("Failed to write data");
        close(fd);
        return;
    }
}

void readBlockFromZone(int fd, int32_t* &array, size_t size, unsigned long long wp) {
    // 打开文件进行读取
    if (fd == -1) {
        perror("Failed to open file for reading");
        return;
    }

    // 分配内存给数组
	if (array == nullptr) {
		array = new int32_t[size];
	}
    // 

    // 读取数组数据到内存中
    size_t bytesToRead = size * sizeof(int32_t);
    if (pread(fd, array, bytesToRead, wp) != bytesToRead) {
        perror("Failed to read data");
        delete[] array;
        close(fd);
        return;
    }

    // 关闭文件
    // close(fd);
}

void validBlock(unsigned long long	wp) {
	int32_t zone_idx = wp / zns.zone_size; 
	int32_t blk_idx = (wp % zns.zone_size) / ZONE_BLK_SIZE;
	global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status = true;
	// global_zone_vec[zone_idx]
	// zns.zones_info[zone_idx].blk_info_vec[blk_idx].status = 1;
}

void invalidBlock(unsigned long long wp) {
	int32_t zone_idx = wp / zns.zone_size; 
	int32_t blk_idx = (wp % zns.zone_size) / ZONE_BLK_SIZE;
	// cout << "wp: " << wp << " zone size: " << zns.zone_size << 
	
	global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status = false;
	// cout << "Invalid zone idx: " << zone_idx << " blk idx: " << blk_idx 
	// << " status: " << global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status << endl;
}

bool test_block(int zone_idx, int blk_idx) 
{
	return global_zone_vec[zone_idx]->blk_info_vec[blk_idx].status;
}

// 记录zone中blk在对应顶点blk_list中的位置，用于回溯时精准定位并更改vertex metadata
void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp)
{
    // cout << "zone size: " << zns.info.zone_size << endl;
    int32_t zone_idx = wp / zns.info.zone_size; // 第几个zone
    int32_t info_idx = (wp % zns.info.zone_size) / ZONE_BLK_SIZE; // blk_info_vec的index
    global_zone_vec[zone_idx]->blk_info_vec[info_idx].status = true;
    global_zone_vec[zone_idx]->blk_info_vec[info_idx].src = vid;
    global_zone_vec[zone_idx]->blk_info_vec[info_idx].blkIndex = blk_idx;
	// cout << "mark zone-- zone: " << zone_idx << " blk: " << 
}

// test_block(int zone_idx, int blk_idx) 


void cacheBlock(Vertex *v, int32_t blk_idx)
{
	v->blk_list[blk_idx].buf = new int32_t[PM_BLK_SIZE]; // 分配内存
	readBlockFromZone(zns.dev_fd, v->blk_list[blk_idx].buf, PM_BLK_SIZE, v->blk_list[blk_idx].global_wp); // 从wp开始读取一个block
	invalidBlock(v->blk_list[blk_idx].global_wp); // 读取出来后无效化zone中的blk
	global_cache_blk++;
}

void releaseBlock(Vertex *v, int32_t blk_idx)
{
	delete[] v->blk_list[blk_idx].buf;
	v->blk_list[blk_idx].buf = nullptr;
	global_cache_blk--;
}

#endif

