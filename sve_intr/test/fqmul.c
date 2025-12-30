#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <arm_sve.h>

#define QINV -3327
#define KYBER_Q 3329

int16_t fqmul(int16_t, int16_t);
int16_t montgomery_reduce(int32_t);

svint16_t fqmul_sv(svint16_t, svint16_t, svbool_t);
svint16_t montgomery_reduce_sv(svint16_t, svint16_t, svbool_t);

uint32_t fqmul_test(uint32_t num_tests) {
    srand(time(NULL));
    uint32_t vio_cnt = 0;
    for(uint32_t i=0; i < num_tests; i++) {
        int16_t a = rand();
        int16_t b = rand();
        int16_t res = fqmul(a,b);
        svbool_t pred = svwhilelt_b16_u32(0,1);
        svint16_t sva = svld1_s16(pred, &a);
        svint16_t svb = svld1_s16(pred, &b);
        svint16_t svresreg = fqmul_sv(sva, svb, pred);
        int16_t svres = 0;
        svst1_s16(pred, &svres, svresreg);

        if (res != svres) {
            printf("int case: %d | sv case %d",res,svres);
            vio_cnt++;
        }
    }
    return vio_cnt;
}


int main(int argc, char** argv) {
    if (argc == 2) {
        unsigned int num_tests = atoi(argv[1]);
        printf("num tests %d",num_tests);
        printf("Final vio count %d",fqmul_test(num_tests));
        printf("\n");
    }
    else if (argc == 3) {
        int16_t a = atoi(argv[1]);
        int16_t b = atoi(argv[2]);
        int16_t res = fqmul(a,b);
        svbool_t pred = svwhilelt_b16_u32(0,1);
        svint16_t sva = svld1_s16(pred, &a);
        svint16_t svb = svld1_s16(pred, &b);
        svint16_t svresreg = fqmul_sv(sva, svb, pred);
        int16_t svres = 0;
        svst1_s16(pred, &svres, svresreg);
        printf("res: %d svres: %d\n",res,svres);
    }
}

int16_t fqmul(int16_t a, int16_t b){
    return montgomery_reduce((int32_t) a*b);
}
int16_t montgomery_reduce(int32_t a)
{
  int16_t t;

  t = (int16_t)a*QINV;
  t = (a - (int32_t)t*KYBER_Q) >> 16;
  return t;
}


// int16_t fqmul_sv(int16_t a, int16_t b){
//     svint16_t au, al;
//     svint16_t sa, sb;
//     svbool_t pred = svwhilelt_b16_u32(0, 1);
//     svint16_t res;
//     int16_t final;
//     sa = svld1(pred, &a);
//     sb = svld1(pred, &b);
//     al = svmul_s16_m(pred, sa, sb);
//     au = svmulh_s16_m(pred, sa, sb);
//     res = montgomery_reduce_sv(al, au);
//     svst1_s16(pred, &final, res);
//     return final;
// }

svint16_t montgomery_reduce_sv(svint16_t ua, svint16_t la, svbool_t pred) {
    svint16_t t;
    t = svmul_n_s16_m(pred, la, QINV);
    t = svmulh_n_s16_m(pred, t, KYBER_Q);
    return svsub_s16_m(pred, ua, t);
}

svint16_t fqmul_sv(svint16_t a, svint16_t b, svbool_t pred) {
  svint16_t uprod, lprod;

  lprod = svmul_s16_m(pred, a, b);
  uprod = svmulh_s16_m(pred, a, b);

  return montgomery_reduce_sv(uprod, lprod, pred);
}
