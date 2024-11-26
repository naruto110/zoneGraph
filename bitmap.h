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

// ��ʼ��λͼ
void initializeBitmap(Bitmap32* bitmap) {
    bitmap->word = 0;
}

// ����λͼ�е�ĳһλ
void setBit(Bitmap32* bitmap, int64_t bitIndex) {
    int64_t bitOffset = bitIndex % BITS_PER_WORD;
    bitmap->word |= (1U << bitOffset);
}

// ���λͼ�е�ĳһλ
void clearBit(Bitmap32* bitmap, int64_t bitIndex) {
    int64_t bitOffset = bitIndex % BITS_PER_WORD;
    bitmap->word &= ~(1U << bitOffset);
}

// ���λͼ�е�ĳһλ�Ƿ�����
bool testBit(Bitmap32* bitmap, int64_t bitIndex) {
    int64_t bitOffset = bitIndex % BITS_PER_WORD;
    return (bitmap->word & (1U << bitOffset)) != 0;
}

// �ϲ�����λͼ����λ��
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
//     // ʾ���÷�
//     Bitmap256 myBitmap;

//     initializeBitmap(&myBitmap);

//     // ���õ�5λ�͵�10λ
//     setBit(&myBitmap, 5);
//     setBit(&myBitmap, 10);

//     // ����5λ�͵�8λ�Ƿ�����
//     printf("Bit 5 is set: %d\n", testBit(&myBitmap, 5));
//     printf("Bit 8 is set: %d\n", testBit(&myBitmap, 8));

//     // �����5λ
//     clearBit(&myBitmap, 5);
    
//     // ����5λ�Ƿ�����
//     printf("Bit 5 is set: %d\n", testBit(&myBitmap, 5));

//     return 0;
// }



#endif

