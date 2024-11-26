#include <iostream>
#include <vector>
#include <queue>
#include <mutex>
#include <cstring>
#include <omp.h>

const int NUM_QUEUES = 3; // ��������
const int NUM_BLOCKS = 5; // ÿ�������е��ڴ������
const int BLOCK_SIZE = 1024; // ÿ���ڴ��Ĵ�С

class MemoryBlockQueue {
public:
    MemoryBlockQueue() {
        refill();
    }

    ~MemoryBlockQueue() {
        clear();
    }

    // ��ȡ��һ���ڴ��
    bool getNextBlock(char*& block) {
        std::lock_guard<std::mutex> lock(mtx_);
        if (memoryBlocks.empty()) {
            return false; // û�п��õ��ڴ��
        }
        block = memoryBlocks.front();
        memoryBlocks.pop();
        return true;
    }

    // �������պ���
    void reclaimMemory() {
        std::lock_guard<std::mutex> lock(mtx_);
        // ģ���ͷ��ڴ沢����������
        std::cout << "Reclaiming memory and refilling the queue." << std::endl;
        clear();
        refill();
    }

private:
    std::queue<char*> memoryBlocks; // �洢�ڴ��Ķ���
    std::mutex mtx_; // �����Զ��еķ���

    // ����ڴ��
    void refill() {
        for (int i = 0; i < NUM_BLOCKS; ++i) {
            char* block = new char[BLOCK_SIZE]; // �����ڴ��
            memoryBlocks.push(block);
        }
    }

    // �ͷ��ڴ��
    void clear() {
        while (!memoryBlocks.empty()) {
            char* block = memoryBlocks.front();
            delete[] block; // �ͷ��ڴ��
            memoryBlocks.pop();
        }
    }
};

void requestMemoryBlocks(std::vector<MemoryBlockQueue>& queues, int threadId) {
    char* block = nullptr;

    // ���Դ�ÿ�����л�ȡ�ڴ��
    while (true) {
        bool blockAcquired = false;
        for (int i = 0; i < NUM_QUEUES; i++) {
            if (queues[i].getNextBlock(block)) {
                blockAcquired = true;
                break; // �ɹ���ȡ�ڴ�飬�˳�ѭ��
            }
        }

        // ����Ƿ������ڴ��
        if (blockAcquired) {
            std::cout << "Thread " << threadId << " is using memory block at " << static_cast<void*>(block) << std::endl;

            // ģ����ڴ��Ĵ���
            #pragma omp critical
            {
                memset(block, threadId, BLOCK_SIZE); // ���߳�ID����ڴ��
            }

            // ������ɺ�����
            std::cout << "Thread " << threadId << " finished using memory block." << std::endl;
        } else {
            // ��������ڴ�鶼��ʹ���꣬������������
            std::cout << "All memory blocks are in use. Thread " << threadId << " initiating garbage collection." << std::endl;
            for (int i = 0; i < NUM_QUEUES; ++i) {
                queues[i].reclaimMemory(); // ������������
            }
        }
    }
}

int main() {
    std::vector<MemoryBlockQueue> queues(NUM_QUEUES);

    #pragma omp parallel
    {
        int threadId = omp_get_thread_num(); // ��ȡ��ǰ�߳�ID
        requestMemoryBlocks(queues, threadId); // �����ڴ��
    }

    return 0;
}
