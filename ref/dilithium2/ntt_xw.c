#include <memory.h>
#include <arm_neon.h>
#include "ntt_xw.h"
#include "reduce.h"

#define Q 8380417           // Q
#define Q_INV 58728449      // Q 逆
#define MONT 2365951        // R = 2^32 R^2 MOD Q = 2365951
#define N_H 128
#define N 256
#define ROOT 3073009        // N-root 原根
#define R_INV 6635910       // ROOT 的模拟
#define INV 16382           // 256 的模拟
#define R_ROOT 4450022      // ROOT 转化成蒙哥马利域上元素
#define R_R_INV 4540456     // R_INV 转化成蒙哥马利域上元素



int32_t rev[N] = {
        0, 128, 64, 192, 32, 160, 96, 224,
        16, 144, 80, 208, 48, 176, 112, 240,
        8, 136, 72, 200, 40, 168, 104, 232,
        24, 152, 88, 216, 56, 184, 120, 248,
        4, 132, 68, 196, 36, 164, 100, 228,
        20, 148, 84, 212, 52, 180, 116, 244,
        12, 140, 76, 204, 44, 172, 108, 236,
        28, 156, 92, 220, 60, 188, 124, 252,
        2, 130, 66, 194, 34, 162, 98, 226,
        18, 146, 82, 210, 50, 178, 114, 242,
        10, 138, 74, 202, 42, 170, 106, 234,
        26, 154, 90, 218, 58, 186, 122, 250,
        6, 134, 70, 198, 38, 166, 102, 230,
        22, 150, 86, 214, 54, 182, 118, 246,
        14, 142, 78, 206, 46, 174, 110, 238,
        30, 158, 94, 222, 62, 190, 126, 254,
        1, 129, 65, 193, 33, 161, 97, 225,
        17, 145, 81, 209, 49, 177, 113, 241,
        9, 137, 73, 201, 41, 169, 105, 233,
        25, 153, 89, 217, 57, 185, 121, 249,
        5, 133, 69, 197, 37, 165, 101, 229,
        21, 149, 85, 213, 53, 181, 117, 245,
        13, 141, 77, 205, 45, 173, 109, 237,
        29, 157, 93, 221, 61, 189, 125, 253,
        3, 131, 67, 195, 35, 163, 99, 227,
        19, 147, 83, 211, 51, 179, 115, 243,
        11, 139, 75, 203, 43, 171, 107, 235,
        27, 155, 91, 219, 59, 187, 123, 251,
        7, 135, 71, 199, 39, 167, 103, 231,
        23, 151, 87, 215, 55, 183, 119, 247,
        15, 143, 79, 207, 47, 175, 111, 239,
        31, 159, 95, 223, 63, 191, 127, 255,
};

int32_t g[N];
int32_t gi[N];

void NTT(int32_t *a)
{
    int32_t temp_a[N];
    for (int i = 0; i < N_H; i++) {
        temp_a[i] = (a[rev[2 * i]] + a[rev[2 * i + 1]]);
        temp_a[i + N_H] = (a[rev[2 * i]] - a[rev[2 * i + 1]] + Q);
    }

    memcpy(a, temp_a, sizeof(temp_a));
    for (int i = 1; i < 8; i++){
        int shift = 8 - 1 - i;
        for (int j = 0; j < N_H; j++) {

            int P = (j >> shift) << shift;
            int32_t odd = montgomery_reduce((int64_t)a[2 * j + 1] * g[P]);
            temp_a[j] = (a[2 * j] + odd);
            temp_a[j + N_H] = (a[2 * j] - odd + Q);
        }
        memcpy(a, temp_a, sizeof(temp_a));
    }
}

void iNTT(int32_t *a)
{
    int32_t temp_a[N];
    for (int i = 0; i < N_H; i++) {
        temp_a[i] = (a[rev[2 * i]] + a[rev[2 * i + 1]]);
        temp_a[i + N_H] = (a[rev[2 * i]] - a[rev[2 * i + 1]] + Q);
    }

    memcpy(a, temp_a, sizeof(temp_a));
    for (int i = 1; i < 8; i++){
        int shift = 8 - 1 - i;
        for (int j = 0; j < N_H; j++) {
            int P = (j >> shift) << shift;
            int32_t odd = montgomery_reduce((int64_t)a[2 * j + 1] * gi[P]);
            temp_a[j] = (a[2 * j] + odd);
            temp_a[j + N_H] = (a[2 * j] - odd + Q);
        }
        memcpy(a, temp_a, sizeof(temp_a));
    }

    for (int i = 0; i < N; i++) {
        a[i] = montgomery_reduce((int64_t)a[i] * INV);
    }
}
