#include <iostream>
#include <vector>
#include <queue>
#include <mutex>
#include <cstring>
#include <omp.h>

const int NUM_QUEUES = 3; // 队列数量
const int NUM_BLOCKS = 5; // 每个队列中的内存块数量
const int BLOCK_SIZE = 1024; // 每个内存块的大小

class MemoryBlockQueue {
public:
    MemoryBlockQueue() {
        refill();
    }

    ~MemoryBlockQueue() {
        clear();
    }

    // 获取下一个内存块
    bool getNextBlock(char*& block) {
        std::lock_guard<std::mutex> lock(mtx_);
        if (memoryBlocks.empty()) {
            return false; // 没有可用的内存块
        }
        block = memoryBlocks.front();
        memoryBlocks.pop();
        return true;
    }

    // 垃圾回收函数
    void reclaimMemory() {
        std::lock_guard<std::mutex> lock(mtx_);
        // 模拟释放内存并重新填充队列
        std::cout << "Reclaiming memory and refilling the queue." << std::endl;
        clear();
        refill();
    }

private:
    std::queue<char*> memoryBlocks; // 存储内存块的队列
    std::mutex mtx_; // 保护对队列的访问

    // 填充内存块
    void refill() {
        for (int i = 0; i < NUM_BLOCKS; ++i) {
            char* block = new char[BLOCK_SIZE]; // 创建内存块
            memoryBlocks.push(block);
        }
    }

    // 释放内存块
    void clear() {
        while (!memoryBlocks.empty()) {
            char* block = memoryBlocks.front();
            delete[] block; // 释放内存块
            memoryBlocks.pop();
        }
    }
};

void requestMemoryBlocks(std::vector<MemoryBlockQueue>& queues, int threadId) {
    char* block = nullptr;

    // 尝试从每个队列获取内存块
    while (true) {
        bool blockAcquired = false;
        for (int i = 0; i < NUM_QUEUES; i++) {
            if (queues[i].getNextBlock(block)) {
                blockAcquired = true;
                break; // 成功获取内存块，退出循环
            }
        }

        // 检查是否获得了内存块
        if (blockAcquired) {
            std::cout << "Thread " << threadId << " is using memory block at " << static_cast<void*>(block) << std::endl;

            // 模拟对内存块的处理
            #pragma omp critical
            {
                memset(block, threadId, BLOCK_SIZE); // 用线程ID填充内存块
            }

            // 处理完成后的输出
            std::cout << "Thread " << threadId << " finished using memory block." << std::endl;
        } else {
            // 如果所有内存块都被使用完，调用垃圾回收
            std::cout << "All memory blocks are in use. Thread " << threadId << " initiating garbage collection." << std::endl;
            for (int i = 0; i < NUM_QUEUES; ++i) {
                queues[i].reclaimMemory(); // 进行垃圾回收
            }
        }
    }
}

int main() {
    std::vector<MemoryBlockQueue> queues(NUM_QUEUES);

    #pragma omp parallel
    {
        int threadId = omp_get_thread_num(); // 获取当前线程ID
        requestMemoryBlocks(queues, threadId); // 请求内存块
    }

    return 0;
}
