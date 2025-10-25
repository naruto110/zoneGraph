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
	global_empty_zone = zns.info.nr_zones - 1;  // ���һ��zoneλ����zone
	global_num_queues = zns.nr_zones - 1; // ����һ����Ϊbuffer

	int total_blks = zns.info.zone_size / ZONE_BLK_SIZE;

	for (int i = 0; i < zns.nr_zones; i++) {
		ZoneStatus *z = new ZoneStatus(total_blks, 0);
		z->blk_info_vec.reserve(total_blks);  // Ĭ�϶�����Ч��
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
		1. ��������zone���鿴��Чblk����
		2. ������Чblk����л���
		3. �Ƚ���Ч����Ǩ������zone�У���resetԭ��zone
			a. ��Ҫ�޸�zone bitmap
			b. ��Ҫ�޸�vertex metadata
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
		wp = global_empty_zone * zns.zone_size;  // ����zone��ʼλ��
		rp = zns.zone_size * zidx;

		int tmp_used_blk = 0;
		// �����ݴ���empty zone�У���������blk����graph_blk_que
		for (int blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
			if (test_block(zidx, blkidx)) { // ��ǰΪvalid blk
				// cout << "Collect zone: " << zidx << " block: " << blkidx << endl;
				// rp = zns.zone_size * zidx + ZONE_BLK_SIZE * blkidx;
				readBlockFromZone(zns.dev_fd, array, PM_BLK_SIZE, rp);
				writeBlockToZone(zns.dev_fd, array, PM_BLK_SIZE, wp);

				// for (int i=0; i < PM_BLK_SIZE; i++) {
				// 	cout << array[i] << ",";
				// }
				// cout << endl;

				// ����vertex metadata
				int32_t vid = global_zone_vec[zidx]->blk_info_vec[blkidx].src;
				int32_t ori_blk_idx = global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex;
				Vertex *v = global_vectex_vec[vid];
				v->blk_list[ori_blk_idx].global_wp = wp; // ���ڸ�blk�Ѿ���д��zone�У���metadata�е�bufΪnullptr
				global_zone_vec[zidx]->blk_info_vec[blkidx].status = false; // ��Ч��ԭ��blk
				// �޸ĸ��º�zone�ı�ǩ��vid��blkidx
				mark_zone(vid, ori_blk_idx, wp);
				wp += ZONE_BLK_SIZE; 
				tmp_used_blk++;
			}
			rp += ZONE_BLK_SIZE;
		}
		// ������zoneʣ�ಿ�ּ���graph_blk_que��
		for (int i = tmp_used_blk; i < zns.nr_blocks_perzone; i++) {
			pm_blk blk;
			blk.global_wp = wp;
			blk.buf = nullptr;
			graph_blk_que.push(blk);
			wp += ZONE_BLK_SIZE; 
		}
		// ����zone reset
		zbd_reset_zones(zns.dev_fd, (zidx * zns.zone_size), zns.zone_size);
		global_empty_zone = zidx;
	}
	delete[] array;
}

void garbage_collection_new() {
	/*
		��Բ�ͬzone�ռ����������ص�
	*/
	double start_t = get_current_time();
	cout << "garbage collection begin! global_empty_zone:" << global_empty_zone << endl;
	unsigned long long wp = 0;
	unsigned long long rp = 0;
	int32_t pre_empty_zone = global_empty_zone;

	// �������zone queue
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
		wp = (unsigned long long)global_empty_zone * (unsigned long long)zns.zone_size;  // ����zone��ʼλ��
		rp = (unsigned long long)zns.zone_size * (unsigned long long)zidx;

		int tmp_used_blk = 0;
		// �����ݴ���empty zone�У���������blk����graph_blk_que
		for (int blkidx = 0; blkidx < zns.nr_blocks_perzone; blkidx++) {
			if (test_block(zidx, blkidx)) { // ��ǰΪvalid blk
				// cout << "Collect zone: " << zidx << " block: " << blkidx << endl;
				// rp = zns.zone_size * zidx + ZONE_BLK_SIZE * blkidx;
				// cout << "wp: " << wp << endl;
				readBlockFromZone(zns.dev_fd, array, PM_BLK_SIZE, rp);
				writeBlockToZone(zns.dev_fd, array, PM_BLK_SIZE, wp);

				// for (int i=0; i < PM_BLK_SIZE; i++) {
				// 	cout << array[i] << ",";
				// }
				// cout << endl;

				// ����vertex metadata
				int32_t vid = global_zone_vec[zidx]->blk_info_vec[blkidx].src;
				int32_t ori_blk_idx = global_zone_vec[zidx]->blk_info_vec[blkidx].blkIndex;
				Vertex *v = global_vectex_vec[vid];
				v->blk_list[ori_blk_idx].global_wp = wp; // ���ڸ�blk�Ѿ���д��zone�У���metadata�е�bufΪnullptr
				global_zone_vec[zidx]->blk_info_vec[blkidx].status = false; // ��Ч��ԭ��blk
				// �޸ĸ��º�zone�ı�ǩ��vid��blkidx
				mark_zone(vid, ori_blk_idx, wp);
				wp += ZONE_BLK_SIZE; 
				tmp_used_blk++;
			}
			rp += ZONE_BLK_SIZE;
		}
		// ������zoneʣ�ಿ�ּ���graph_blk_que��
		for (int i = tmp_used_blk; i < zns.nr_blocks_perzone; i++) {
			global_blk_queues[zidx]->memoryWPs.push(wp);
			// pm_blk blk;
			// blk.global_wp = wp;
			// blk.buf = nullptr;
			// graph_blk_que.push(blk);
			wp += ZONE_BLK_SIZE; 
		}
		// ����zone reset
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
		������zone������blocks, �������Ƿ�ͬһzone���ڽڵ�����ֲ�
	*/
	// cout << "generate zone blocks" << endl;
	while (blk_num--) {
		pm_blk blk;
		blk.global_wp = global_blk_wp; // ����global_blk_wp���Լ����zone id�;ֲ�wp
		// blk.blk_id = global_blk_id;
		blk.buf = nullptr;
		global_blk_id++;
		global_blk_wp += ZONE_BLK_SIZE; // ƫ�Ƶ����ֽ���
		graph_blk_que.push(blk);
		// if (global_blk_wp >= global_threshold_wp) { // ��Ҫ������������

		// }
	}
	// cout << "generate zone blocks done, total blocks: " << global_blk_id << endl;

	return 1;
}

void generate_block_opt(void)
{
	// ��ʼ���ڴ����
	unsigned long long wp = 0;
    for (int zoneidx = 0; zoneidx < global_num_queues; zoneidx++) {
		MemoryWPQueue* b = new MemoryWPQueue(wp);
        global_blk_queues.push_back(b);
		wp += zns.zone_size;
    }
}

// ���������ڴ��
unsigned long long requestMemoryBlocks(int thread_id) 
{
    unsigned long long wp = 0;

	if (global_blk_queues[(thread_id % global_num_queues)]->getNextBlock(wp)) {
		return wp;
	}
    // ���Դ�ÿ�����л�ȡ�ڴ��
	int iters = 32;
	bool blockAcquired = false;
	while ((iters>0) && (!blockAcquired)) {
		for (int i = 0; i < global_num_queues; i++) {
			if (global_blk_queues[i]->getNextBlock(wp)) {  // ��ʱ�Ѿ��������ˣ��ö���queue���ᱻlock���ڴ˺�������ɶ�zone��д�����Ǵ��е�
				blockAcquired = true;
				break; // �ɹ���ȡ�ڴ�飬�˳�ѭ��
			}
		}
		iters--;
	}

	// ����Ƿ������ڴ��
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

    // ����������д���ļ�
    size_t bytesToWrite = size * sizeof(int32_t);
    if (pwrite(fd, array, bytesToWrite, wp) != bytesToWrite) {
        perror("Failed to write data");
        close(fd);
        return;
    }
}

void readBlockFromZone(int fd, int32_t* &array, size_t size, unsigned long long wp) {
    // ���ļ����ж�ȡ
    if (fd == -1) {
        perror("Failed to open file for reading");
        return;
    }

    // �����ڴ������
	if (array == nullptr) {
		array = new int32_t[size];
	}
    // 

    // ��ȡ�������ݵ��ڴ���
    size_t bytesToRead = size * sizeof(int32_t);
    if (pread(fd, array, bytesToRead, wp) != bytesToRead) {
        perror("Failed to read data");
        delete[] array;
        close(fd);
        return;
    }

    // �ر��ļ�
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

// ��¼zone��blk�ڶ�Ӧ����blk_list�е�λ�ã����ڻ���ʱ��׼��λ������vertex metadata
void mark_zone(int32_t vid, int32_t blk_idx, unsigned long long wp)
{
    // cout << "zone size: " << zns.info.zone_size << endl;
    int32_t zone_idx = wp / zns.info.zone_size; // �ڼ���zone
    int32_t info_idx = (wp % zns.info.zone_size) / ZONE_BLK_SIZE; // blk_info_vec��index
    global_zone_vec[zone_idx]->blk_info_vec[info_idx].status = true;
    global_zone_vec[zone_idx]->blk_info_vec[info_idx].src = vid;
    global_zone_vec[zone_idx]->blk_info_vec[info_idx].blkIndex = blk_idx;
	// cout << "mark zone-- zone: " << zone_idx << " blk: " << 
}

// test_block(int zone_idx, int blk_idx) 


void cacheBlock(Vertex *v, int32_t blk_idx)
{
	v->blk_list[blk_idx].buf = new int32_t[PM_BLK_SIZE]; // �����ڴ�
	readBlockFromZone(zns.dev_fd, v->blk_list[blk_idx].buf, PM_BLK_SIZE, v->blk_list[blk_idx].global_wp); // ��wp��ʼ��ȡһ��block
	invalidBlock(v->blk_list[blk_idx].global_wp); // ��ȡ��������Ч��zone�е�blk
	global_cache_blk++;
}

void releaseBlock(Vertex *v, int32_t blk_idx)
{
	delete[] v->blk_list[blk_idx].buf;
	v->blk_list[blk_idx].buf = nullptr;
	global_cache_blk--;
}

#endif

