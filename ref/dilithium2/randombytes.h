// 在window和Linux下不同的实现，Linux下使用内核获取具体的随机数
#ifndef RANDOMBYTES_H
#define RANDOMBYTES_H

#include <stddef.h>
#include <stdint.h>

void randombytes(uint8_t *out, size_t outlen);

#endif
