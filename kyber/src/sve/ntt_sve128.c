#include "ntt_sve128.h"

void ntt_sve128(int16_t r[256]) {
    level0_sve128(r,zetas,0);
    level0_sve128(r,zetas,1);
    levels1t6_sve128(r,zetas,0);
    levels1t6_sve128(r,zetas,1);
}
