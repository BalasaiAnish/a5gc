#include "ntt.h"

inline svint16_t montgomery_reduce(svint32_t a, svbool_t pred) {
  svint16_t t_16;
  svint32_t t_32;
  svint16_t ua_16;
  svint32_t u, ua_32;

  ua_32 = svmulh_n_s32_m(pred, a, QINV);
  ua_16 = svreinterpret_s16_s32(ua_32);
  u = svreinterpret_s32_s16(ua_16);
  t_32 = svmulh_n_s32_m(pred, u, KYBER_Q);
  t_32 = svsub_s32_m(pred, a, t_32);
  t_32 = svasr_n_s32_m(pred, t_32, 16);
  t_16 = svreinterpret_s16_s32(t_32);
  return t_16;
}

inline svint16_t fqmul(svint16_t a, svint16_t b, svbool_t pred) {
  svint32_t a_32, b_32;
  svint32_t prod;

  a_32 = svreinterpret_s32_s16(a);
  b_32 = svreinterpret_s32_s16(b);
  prod = svmulh_s32_m(pred, a_32, b_32);

  return montgomery_reduce(prod, pred);
}

inline svint16_t svcat2_s16(svint16_t a, svint16_t b, unsigned int vl) {
    svbool_t half0;
    half0 = svwhilelt_b16_u32(0, vl/2);
    a = svsplice_s16(half0, a, b);
    return a;
}

inline svint16_t svcat4_s16(svint16_t a, svint16_t b, svint16_t c, svint16_t d, unsigned int vl) {
    svbool_t half;
    svbool_t quarter;

    half = svwhilelt_b16_u32(0, vl/2);
    quarter = svwhilelt_b16_u32(0, vl/4);

    a = svsplice_s16(quarter, a, b);
    c = svsplice_s16(quarter, c, d);
    a = svsplice_s16(half, a, c);
    return a;
}

inline svint16_t svcat8_s16(svint16_t reg0, svint16_t reg1, svint16_t reg2, svint16_t reg3,
    svint16_t reg4, svint16_t reg5, svint16_t reg6, svint16_t reg7, unsigned int vl) {
    svbool_t half = svwhilelt_b16_u32(0, vl/2);
    svbool_t quarter = svwhilelt_b16_u32(0, vl/4);
    svbool_t eigth = svwhilelt_b16_u32(0, vl/8);

    reg0 = svsplice_s16(eigth, reg0, reg1);
    reg2 = svsplice_s16(eigth, reg2, reg3);
    reg4 = svsplice_s16(eigth, reg4, reg5);
    reg6 = svsplice_s16(eigth, reg6, reg7);

    reg0 = svsplice_s16(quarter, reg0, reg2);
    reg2 = svsplice_s16(quarter, reg4, reg6);

    reg0 = svsplice_s16(half, reg0, reg2);

    return reg0;
}

inline svint16_t svcat16_s16(svint16_t reg0, svint16_t reg1, svint16_t reg2, svint16_t reg3,
    svint16_t reg4, svint16_t reg5, svint16_t reg6, svint16_t reg7,
    svint16_t reg8, svint16_t reg9, svint16_t reg10, svint16_t reg11,
    svint16_t reg12, svint16_t reg13, svint16_t reg14, svint16_t reg15,
    unsigned int vl) {
    svbool_t half = svwhilelt_b16_u32(0, vl/2);
    svbool_t quarter = svwhilelt_b16_u32(0, vl/4);
    svbool_t eigth = svwhilelt_b16_u32(0, vl/8);
    svbool_t sixteenth = svwhilelt_b16_u32(0, vl/16);

    reg0 = svsplice_s16(sixteenth, reg0, reg1);
    reg2 = svsplice_s16(sixteenth, reg2, reg3);
    reg4 = svsplice_s16(sixteenth, reg4, reg5);
    reg6 = svsplice_s16(sixteenth, reg6, reg7);
    reg8 = svsplice_s16(sixteenth, reg8, reg9);
    reg10 = svsplice_s16(sixteenth, reg10, reg11);
    reg12 = svsplice_s16(sixteenth, reg12, reg13);
    reg14 = svsplice_s16(sixteenth, reg14, reg15);

    reg0 = svsplice_s16(eigth, reg0, reg2);
    reg2 = svsplice_s16(eigth, reg4, reg6);
    reg4 = svsplice_s16(eigth, reg8, reg10);
    reg6 = svsplice_s16(eigth, reg12, reg14);

    reg0 = svsplice_s16(quarter, reg0, reg2);
    reg2 = svsplice_s16(quarter, reg4, reg6);

    reg0 = svsplice_s16(half, reg0, reg2);

    return reg0;
}

void ntt(int16_t r[256]) {
  unsigned int len, start, j, k;
  int16_t t, zeta;

  k = 1;
  for(len = 128; len >= 2; len >>= 1) {
    for(start = 0; start < 256; start = j + len) {
      zeta = zetas[k++];
      for(j = start; j < start + len; j++) {
        t = zeta*r[j+len];
        r[j + len] = r[j] - t;
        r[j] = r[j] + t;
      }
    }
  }
}

inline void butterfly_level0(int16_t r[256], int16_t zeta, unsigned int start, unsigned int len) {
    svbool_t pred;
    svint16_t a, b;
    svint16_t svzeta;
    svint16_t t;
    uint32_t i = 0;
    unsigned int vl = svcnth();

    pred = svwhilelt_b16_u32(i, len);
    svzeta = svdup_s16(zeta);
    while (svptest_first(svptrue_b16(), pred)) {
        a = svld1(pred,&r[i+start]);
        b = svld1(pred,&r[i+start+len]);
        t = fqmul(svzeta, b, pred);
        b = svsub_s16_m(pred, a, t);
        a = svadd_s16_m(pred, a, t);
        svst1(pred, &r[i+start], a);
        svst1(pred, &r[i+start+len], b);
        i += vl;
        pred = svwhilelt_b16_u32(i, len);
    }
}

inline void butterfly_level1(int16_t r[256], int16_t zeta0, int16_t zeta1, unsigned int start0, unsigned int start1, unsigned int len) {
    svbool_t pred;
    svbool_t pred_half0, pred_half1;
    svint16_t a, a0, a1, b, b0, b1;
    svint16_t svzeta, svzeta0, svzeta1;
    svint16_t t;
    uint32_t i = 0;
    unsigned int vl = svcnth();
    // looping logic might be suspect, check pred stuff
    pred = svwhilelt_b16_u32(i, len);
    pred_half0 = svwhilelt_b16_u32(0, len/2);
    pred_half1 = svwhilege_b16_u32(0, len/2);
    svzeta = svcat2_s16(svzeta0, svzeta1, vl);

    while (svptest_first(svptrue_b16(), pred)) {
        a0 = svld1(pred_half0,&r[i+start0]);
        b0 = svld1(pred_half0,&r[i+start0+len]);
        a1 = svld1(pred_half0,&r[i+start1]);
        b1 = svld1(pred_half0,&r[i+start1+len]);
        a = svcat2_s16(a0, a1, vl);
        b = svcat2_s16(b0, b1, vl);

        t = fqmul(svzeta, b, pred);
        b = svsub_s16_m(pred, a, t);
        a = svadd_s16_m(pred, a, t);

        svst1(pred_half0, &r[i+start0], a);
        svst1(pred_half0, &r[i+start0+len], b);
        svst1(pred_half1, &r[i+start1], a);
        svst1(pred_half1, &r[i+start1+len], b);
        i += vl;
        pred = svwhilelt_b16_u32(i, len);
    }
}

inline void butterfly_level2(int16_t r[256], int16_t zeta0, int16_t zeta1, int16_t zeta2, int16_t zeta3,
    unsigned int start0, unsigned int start1, unsigned int start2, unsigned int start3, unsigned int len) {
    svbool_t pred;
    svbool_t pred_half0, pred_quarter0, pred_quarter1, pred_quarter2, pred_quarter3;
    svint16_t a, a0, a1, a2, a3, b, b0, b1, b2, b3;
    svint16_t svzeta, svzeta0, svzeta1, svzeta2, svzeta3;
    svint16_t t;
    uint32_t i = 0;
    unsigned int vl = svcnth();
    // looping logic might be suspect, check pred stuff
    pred = svwhilelt_b16_u32(i, len);
    pred_half0 = svwhilelt_b16_u32(0, len/2);
    pred_quarter0 = svwhilelt_b16_u32(0, len/4);
    pred_quarter1 = svwhilelt_b16_u32(len/4, len/2);
    pred_quarter2 = svwhilelt_b16_u32(len/2, 3*len/4);
    pred_quarter3 = svwhilelt_b16_u32(3*len/4, len);
    svzeta0 = svdup_s16(zeta0);
    svzeta1 = svdup_s16(zeta1);
    svzeta2 = svdup_s16(zeta2);
    svzeta3 = svdup_s16(zeta3);
    svzeta = svcat4_s16(svzeta0, svzeta1, svzeta2, svzeta3, len);

    while (svptest_first(svptrue_b16(), pred)) {
        a0 = svld1(pred_quarter0,&r[i+start0]);
        b0 = svld1(pred_quarter0,&r[i+start0+len]);
        a1 = svld1(pred_quarter0,&r[i+start1]);
        b1 = svld1(pred_quarter0,&r[i+start1+len]);
        a2 = svld1(pred_quarter0,&r[i+start2]);
        b2 = svld1(pred_quarter0,&r[i+start2+len]);
        a3 = svld1(pred_quarter0,&r[i+start3]);
        b3 = svld1(pred_quarter0,&r[i+start3+len]);
        a = svcat4_s16(a0, a1, a2, a3, vl);
        b = svcat4_s16(b0, b1, b2, b3, vl);

        t = fqmul(svzeta, b, pred);
        b = svsub_s16_m(pred, a, t);
        a = svadd_s16_m(pred, a, t);

        svst1(pred_quarter0,&r[i+start0],a);
        svst1(pred_quarter0,&r[i+start0+len],b);
        svst1(pred_quarter1,&r[i+start1],a);
        svst1(pred_quarter1,&r[i+start1+len],b);
        svst1(pred_quarter2,&r[i+start2],a);
        svst1(pred_quarter2,&r[i+start2+len],b);
        svst1(pred_quarter3,&r[i+start3],a);
        svst1(pred_quarter3,&r[i+start3+len],b);

        i += vl;
        pred = svwhilelt_b16_u32(i, len);
    }
}

inline void butterfly_level3(int16_t r[256], int16_t zeta0, int16_t zeta1, int16_t zeta2, int16_t zeta3,
    int16_t zeta4, int16_t zeta5, int16_t zeta6, int16_t zeta7,
    unsigned int start0, unsigned int start1, unsigned int start2, unsigned int start3,
    unsigned int start4, unsigned int start5, unsigned int start6, unsigned int start7,
    unsigned int len) {
    svbool_t pred;
    svbool_t pred_half0, pred_eigth0, pred_eigth1, pred_eigth2, pred_eigth3;
    svbool_t pred_eigth4, pred_eigth5, pred_eigth6, pred_eigth7;
    svint16_t a, a0, a1, a2, a3, a4, a5, a6, a7;
    svint16_t b, b0, b1, b2, b3, b4, b5, b6, b7;
    svint16_t svzeta, svzeta0, svzeta1, svzeta2, svzeta3, svzeta4, svzeta5, svzeta6, svzeta7;
    svint16_t t;
    uint32_t i = 0;
    unsigned int vl = svcnth();
    // looping logic might be suspect, check pred stuff
    pred = svwhilelt_b16_u32(i, len);
    pred_half0 = svwhilelt_b16_u32(0, len/2);
    pred_eigth0 = svwhilelt_b16_u32(0, len/8);
    pred_eigth1 = svwhilelt_b16_u32(len/8, len/4);
    pred_eigth2 = svwhilelt_b16_u32(len/4, 3*len/8);
    pred_eigth3 = svwhilelt_b16_u32(3*len/8, len/2);
    pred_eigth4 = svwhilelt_b16_u32(len/2, 5*len/8);
    pred_eigth5 = svwhilelt_b16_u32(5*len/8, 3*len/4);
    pred_eigth6 = svwhilelt_b16_u32(3*len/4, 7*len/8);
    pred_eigth7 = svwhilelt_b16_u32(7*len/8, len);
    svzeta0 = svdup_s16(zeta0);
    svzeta1 = svdup_s16(zeta1);
    svzeta2 = svdup_s16(zeta2);
    svzeta3 = svdup_s16(zeta3);
    svzeta4 = svdup_s16(zeta4);
    svzeta5 = svdup_s16(zeta5);
    svzeta6 = svdup_s16(zeta6);
    svzeta7 = svdup_s16(zeta7);
    svzeta = svcat8_s16(svzeta0, svzeta1, svzeta2, svzeta3, svzeta4, svzeta5, svzeta6, svzeta7, vl);

    while (svptest_first(svptrue_b16(), pred)) {
        a0 = svld1(pred_eigth0,&r[i+start0]);
        b0 = svld1(pred_eigth0,&r[i+start0+len]);
        a1 = svld1(pred_eigth0,&r[i+start1]);
        b1 = svld1(pred_eigth0,&r[i+start1+len]);
        a2 = svld1(pred_eigth0,&r[i+start2]);
        b2 = svld1(pred_eigth0,&r[i+start2+len]);
        a3 = svld1(pred_eigth0,&r[i+start3]);
        b3 = svld1(pred_eigth0,&r[i+start3+len]);
        a4 = svld1(pred_eigth0,&r[i+start4]);
        b4 = svld1(pred_eigth0,&r[i+start4+len]);
        a5 = svld1(pred_eigth0,&r[i+start5]);
        b5 = svld1(pred_eigth0,&r[i+start5+len]);
        a6 = svld1(pred_eigth0,&r[i+start6]);
        b6 = svld1(pred_eigth0,&r[i+start6+len]);
        a7 = svld1(pred_eigth0,&r[i+start7]);
        b7 = svld1(pred_eigth0,&r[i+start7+len]);
        a = svcat8_s16(a0, a1, a2, a3, a4, a5, a6, a7, vl);
        b = svcat8_s16(b0, b1, b2, b3, b4, b5, b6, b7, vl);

        t = fqmul(svzeta, b, pred);
        b = svsub_s16_m(pred, a, t);
        a = svadd_s16_m(pred, a, t);

        svst1(pred_eigth0,&r[i+start0],a);
        svst1(pred_eigth0,&r[i+start0+len],b);
        svst1(pred_eigth1,&r[i+start1],a);
        svst1(pred_eigth1,&r[i+start1+len],b);
        svst1(pred_eigth2,&r[i+start2],a);
        svst1(pred_eigth2,&r[i+start2+len],b);
        svst1(pred_eigth3,&r[i+start3],a);
        svst1(pred_eigth3,&r[i+start3+len],b);
        svst1(pred_eigth4,&r[i+start4],a);
        svst1(pred_eigth4,&r[i+start4+len],b);
        svst1(pred_eigth5,&r[i+start5],a);
        svst1(pred_eigth5,&r[i+start5+len],b);
        svst1(pred_eigth6,&r[i+start6],a);
        svst1(pred_eigth6,&r[i+start6+len],b);
        svst1(pred_eigth7,&r[i+start7],a);
        svst1(pred_eigth7,&r[i+start7+len],b);

        i += vl;
        pred = svwhilelt_b16_u32(i, len);
    }
}
