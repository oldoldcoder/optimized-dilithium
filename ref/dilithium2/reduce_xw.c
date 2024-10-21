#include "reduce_xw.h"

// 适配许伟ntt的蒙哥马利约减算法

int32_t montgomery_reduce_xw(int64_t a){
    int32_t t = (int32_t) a * Q_INV;
    t = (a - (int64_t) t * Q) >> 32;
    return t;
}

int32_t reduce32_xw(int32_t a){
    int32_t t;

    t = (a + (1 << 22)) >> 23;
    t = a - t*Q;
    return t;
}

int32_t caddq_xw(int32_t a){
    a += (a >> 31) & Q;
    return a;
}

int32_t freeze_xw(int32_t a){
    a = reduce32_xw(a);
    a = caddq_xw(a);
    return a;
}