#include "barrett_mul.h"

#define Q 456645

int32_t barrett_mul3(int64_t a, int64_t b) {
    int64_t z = a * b;
    int64_t m1 = (z >> 24);
    int64_t m2 = ((m1 << 25) + (m1 << 15) + (m1 << 5) - (m1 << 2));
    int64_t m3 = (m2 >> 24);
    int64_t t = z - ((m3 << 23) - (m3 << 13) + m3);
    if(t>=Q){
        return (int32_t)(t-Q);
    }
    else{
        return (int32_t)t;
    }
}