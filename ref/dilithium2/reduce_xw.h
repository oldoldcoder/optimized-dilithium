#ifndef REDUCE_XW_H
#define REDUCE_XW_H

#include "stdint.h"

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

int32_t montgomery_reduce_xw(int64_t a);

int32_t reduce32_xw(int32_t a);

int32_t caddq_xw(int32_t a);

int32_t freeze_xw(int32_t a);

#endif