#ifndef __BITMAP_H_
#define __BITMAP_H_

// #include "global.h"

#include <stdint.h>
#include <stdio.h>
#include <iostream>

using namespace std;

#define BITS_PER_WORD 32
#define BITMAP_SIZE 1  // 256 bits / 32 bits per word = 8 words

typedef struct {
    uint32_t word;
} Bitmap32;

// 初始化位图
void initializeBitmap(Bitmap32* bitmap) {
    bitmap->word = 0;
}

// 设置位图中的某一位
void setBit(Bitmap32* bitmap, int64_t bitIndex) {
    int64_t bitOffset = bitIndex % BITS_PER_WORD;
    bitmap->word |= (1U << bitOffset);
}

// 清除位图中的某一位
void clearBit(Bitmap32* bitmap, int64_t bitIndex) {
    int64_t bitOffset = bitIndex % BITS_PER_WORD;
    bitmap->word &= ~(1U << bitOffset);
}

// 检查位图中的某一位是否被设置
bool testBit(Bitmap32* bitmap, int64_t bitIndex) {
    int64_t bitOffset = bitIndex % BITS_PER_WORD;
    return (bitmap->word & (1U << bitOffset)) != 0;
}

// 合并两个位图，按位与
void andBit(Bitmap32* des_bitmap, Bitmap32* src_bitmap) {
    des_bitmap->word = des_bitmap->word | src_bitmap->word;
}

void print_bitmap(Bitmap32* bitmap) {
    for (int i = 0; i < 32; i++) {
        if (testBit(bitmap, i)) {
            cout << 1 << " ";
        } else {
            cout << 0 << " ";
        }
    }
    cout << endl;
}

// int main() {
//     // 示例用法
//     Bitmap256 myBitmap;

//     initializeBitmap(&myBitmap);

//     // 设置第5位和第10位
//     setBit(&myBitmap, 5);
//     setBit(&myBitmap, 10);

//     // 检查第5位和第8位是否被设置
//     printf("Bit 5 is set: %d\n", testBit(&myBitmap, 5));
//     printf("Bit 8 is set: %d\n", testBit(&myBitmap, 8));

//     // 清除第5位
//     clearBit(&myBitmap, 5);
    
//     // 检查第5位是否被设置
//     printf("Bit 5 is set: %d\n", testBit(&myBitmap, 5));

//     return 0;
// }



#endif

