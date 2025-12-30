#include <arm_sve.h>
#include <stdint.h>
#include <stdio.h>

#define QINV -3327
#define KYBER_Q 3329

static inline svint16_t montgomery_reduce(svint16_t ua, svint16_t la, svbool_t pred) {
    svint16_t t, qial, u;
    qial = svmul_n_s16_z(pred, la, QINV);
    u = svmulh_n_s16_z(pred, qial, KYBER_Q);
    t = svsub_s16_m(pred, ua, u);
    return t;
}

static inline svint16_t fqmul(svint16_t a, svint16_t b, svbool_t pred) {
  svint16_t uprod, lprod, res;
  lprod = svmul_s16_z(pred, a, b);
  uprod = svmulh_s16_z(pred, a, b);
  res = montgomery_reduce(uprod, lprod, pred);
  return res;
}

static inline void level0_sve128(int16_t r[256], const int16_t zetas[128], const unsigned int half) {
    svint16_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
    svint16_t r16, r17, r18, r19, r20, r21, r22, r23;
    svint16_t zeta;
    // svint16_t t;
    svbool_t pred = svptrue_b16();

    zeta = svdup_s16(zetas[1]);

    r8 = svld1_s16(pred, &r[half*64+128+0*8]);
    r9 = svld1_s16(pred, &r[half*64+128+1*8]);
    r10 = svld1_s16(pred, &r[half*64+128+2*8]);
    r11 = svld1_s16(pred, &r[half*64+128+3*8]);
    r12 = svld1_s16(pred, &r[half*64+128+4*8]);
    r13 = svld1_s16(pred, &r[half*64+128+5*8]);
    r14 = svld1_s16(pred, &r[half*64+128+6*8]);
    r15 = svld1_s16(pred, &r[half*64+128+7*8]);

    r16 = fqmul(zeta, r8, pred);
    r17 = fqmul(zeta, r9, pred);
    r18 = fqmul(zeta, r10, pred);
    r19 = fqmul(zeta, r11, pred);
    r20 = fqmul(zeta, r12, pred);
    r21 = fqmul(zeta, r13, pred);
    r22 = fqmul(zeta, r14, pred);
    r23 = fqmul(zeta, r15, pred);

    r0 = svld1_s16(pred, &r[half*64+0*8]);
    r1 = svld1_s16(pred, &r[half*64+1*8]);
    r2 = svld1_s16(pred, &r[half*64+2*8]);
    r3 = svld1_s16(pred, &r[half*64+3*8]);
    r4 = svld1_s16(pred, &r[half*64+4*8]);
    r5 = svld1_s16(pred, &r[half*64+5*8]);
    r6 = svld1_s16(pred, &r[half*64+6*8]);
    r7 = svld1_s16(pred, &r[half*64+7*8]);

    r8 = svsub_s16_x(pred, r0, r16);
    r9 = svsub_s16_x(pred, r1, r17);
    r10 = svsub_s16_x(pred, r2, r18);
    r11 = svsub_s16_x(pred, r3, r19);
    r12 = svsub_s16_x(pred, r4, r20);
    r13 = svsub_s16_x(pred, r5, r21);
    r14 = svsub_s16_x(pred, r6, r22);
    r15 = svsub_s16_x(pred, r7, r23);

    r0 = svadd_s16_x(pred, r0, r16);
    r1 = svadd_s16_x(pred, r1, r17);
    r2 = svadd_s16_x(pred, r2, r18);
    r3 = svadd_s16_x(pred, r3, r19);
    r4 = svadd_s16_x(pred, r4, r20);
    r5 = svadd_s16_x(pred, r5, r21);
    r6 = svadd_s16_x(pred, r6, r22);
    r7 = svadd_s16_x(pred, r7, r23);

    svst1_s16(pred, &r[half*64+0*8], r0);
    svst1_s16(pred, &r[half*64+1*8], r1);
    svst1_s16(pred, &r[half*64+2*8], r2);
    svst1_s16(pred, &r[half*64+3*8], r3);
    svst1_s16(pred, &r[half*64+4*8], r4);
    svst1_s16(pred, &r[half*64+5*8], r5);
    svst1_s16(pred, &r[half*64+6*8], r6);
    svst1_s16(pred, &r[half*64+7*8], r7);

    svst1_s16(pred, &r[half*64+128+0*8], r8);
    svst1_s16(pred, &r[half*64+128+1*8], r9);
    svst1_s16(pred, &r[half*64+128+2*8], r10);
    svst1_s16(pred, &r[half*64+128+3*8], r11);
    svst1_s16(pred, &r[half*64+128+4*8], r12);
    svst1_s16(pred, &r[half*64+128+5*8], r13);
    svst1_s16(pred, &r[half*64+128+6*8], r14);
    svst1_s16(pred, &r[half*64+128+7*8], r15);
}

static inline void levels1t6_sve128(int16_t r[256], const int16_t zetas[128], const unsigned int half) {
    svint16_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
    svint16_t r16, r17, r18, r19, r20, r21, r22, r23, r24, r25, r26, r27, r28, r29, r30, r31;
    svint16_t r32, r33, r34, r35, r36, r37, r38, r39, r40, r41, r42, r43, r44, r45, r46, r47;
    svint16_t r48, r49, r50, r51, r52, r53, r54, r55, r56, r57, r58, r59, r60, r61, r62, r63;
    svint16_t r64, r65, r66, r67, r68, r69, r70, r71, r72, r73, r74, r75, r76, r77, r78, r79;
    svint16_t zeta, zeta0, zeta1, zeta2, zeta3, zeta4, zeta5, zeta6, zeta7;
    svint16_t zeta8, zeta9, zeta10, zeta11, zeta12, zeta13, zeta14, zeta15;
    svint16_t zeta16, zeta17, zeta18, zeta19, zeta20, zeta21, zeta22, zeta23;
    svint16_t zeta24, zeta25, zeta26, zeta27, zeta28, zeta29, zeta30, zeta31;
    // svint16_t t;
    svbool_t pred = svptrue_b16();

    // Level 1

    zeta = svdup_s16(zetas[2+half*1]);

    r8 = svld1_s16(pred, &r[half*128+64+0*8]);
    r9 = svld1_s16(pred, &r[half*128+64+1*8]);
    r10 = svld1_s16(pred, &r[half*128+64+2*8]);
    r11 = svld1_s16(pred, &r[half*128+64+3*8]);
    r12 = svld1_s16(pred, &r[half*128+64+4*8]);
    r13 = svld1_s16(pred, &r[half*128+64+5*8]);
    r14 = svld1_s16(pred, &r[half*128+64+6*8]);
    r15 = svld1_s16(pred, &r[half*128+64+7*8]);

    r16 = fqmul(zeta, r8, pred);
    r17 = fqmul(zeta, r9, pred);
    r18 = fqmul(zeta, r10, pred);
    r19 = fqmul(zeta, r11, pred);
    r20 = fqmul(zeta, r12, pred);
    r21 = fqmul(zeta, r13, pred);
    r22 = fqmul(zeta, r14, pred);
    r23 = fqmul(zeta, r15, pred);

    r0 = svld1_s16(pred, &r[half*128+0*8]);
    r1 = svld1_s16(pred, &r[half*128+1*8]);
    r2 = svld1_s16(pred, &r[half*128+2*8]);
    r3 = svld1_s16(pred, &r[half*128+3*8]);
    r4 = svld1_s16(pred, &r[half*128+4*8]);
    r5 = svld1_s16(pred, &r[half*128+5*8]);
    r6 = svld1_s16(pred, &r[half*128+6*8]);
    r7 = svld1_s16(pred, &r[half*128+7*8]);

    r8 = svsub_s16_x(pred, r0, r16);
    r9 = svsub_s16_x(pred, r1, r17);
    r10 = svsub_s16_x(pred, r2, r18);
    r11 = svsub_s16_x(pred, r3, r19);
    r12 = svsub_s16_x(pred, r4, r20);
    r13 = svsub_s16_x(pred, r5, r21);
    r14 = svsub_s16_x(pred, r6, r22);
    r15 = svsub_s16_x(pred, r7, r23);

    r0 = svadd_s16_x(pred, r0, r16);
    r1 = svadd_s16_x(pred, r1, r17);
    r2 = svadd_s16_x(pred, r2, r18);
    r3 = svadd_s16_x(pred, r3, r19);
    r4 = svadd_s16_x(pred, r4, r20);
    r5 = svadd_s16_x(pred, r5, r21);
    r6 = svadd_s16_x(pred, r6, r22);
    r7 = svadd_s16_x(pred, r7, r23);

    // Level 2

    svbool_t pred_1st_half = svwhilelt_b16_u32(0, 4);
    svbool_t pred_2nd_half = svnot_b_z(pred, pred_1st_half);
    zeta0 = svdup_s16(zetas[4+half*2]);
    zeta1 = svdup_s16(zetas[5+half*2]);

    r16 = fqmul(zeta0, r4, pred);
    r17 = fqmul(zeta0, r5, pred);
    r18 = fqmul(zeta0, r6, pred);
    r19 = fqmul(zeta0, r7, pred);
    r20 = fqmul(zeta1, r12, pred);
    r21 = fqmul(zeta1, r13, pred);
    r22 = fqmul(zeta1, r14, pred);
    r23 = fqmul(zeta1, r15, pred);

    r4 = svsub_s16_x(pred, r0, r16);
    r5 = svsub_s16_x(pred, r1, r17);
    r6 = svsub_s16_x(pred, r2, r18);
    r7 = svsub_s16_x(pred, r3, r19);
    r12 = svsub_s16_x(pred, r8, r20);
    r13 = svsub_s16_x(pred, r9, r21);
    r14 = svsub_s16_x(pred, r10, r22);
    r15 = svsub_s16_x(pred, r11, r23);

    r0 = svadd_s16_x(pred, r0, r16);
    r1 = svadd_s16_x(pred, r1, r17);
    r2 = svadd_s16_x(pred, r2, r18);
    r3 = svadd_s16_x(pred, r3, r19);
    r8 = svadd_s16_x(pred, r8, r20);
    r9 = svadd_s16_x(pred, r9, r21);
    r10 = svadd_s16_x(pred, r10, r22);
    r11 = svadd_s16_x(pred, r11, r23);

    // Level 3

    svbool_t pred_1st_quarter = svwhilelt_b16_u32(0, 2);
    svbool_t pred_2nd_quarter = sveor_b_z(pred_1st_half, pred_1st_quarter, pred_1st_half);
    svbool_t pred_3_quarters = svwhilelt_b16_u32(0, 6);
    svbool_t pred_3rd_quarter = svand_b_z(pred, pred_3_quarters, pred_2nd_half);
    svbool_t pred_4th_quarter = sveor_b_z(pred_2nd_half, pred_3rd_quarter, pred_2nd_half);
    svbool_t pred_last_3_quarters = svnot_b_z(pred, pred_1st_quarter);

    zeta0 = svdup_s16(zetas[8+half*4]);
    zeta1 = svdup_s16(zetas[9+half*4]);
    zeta2 = svdup_s16(zetas[10+half*4]);
    zeta3 = svdup_s16(zetas[11+half*4]);

    r16 = fqmul(zeta0, r2, pred);
    r17 = fqmul(zeta0, r3, pred);
    r18 = fqmul(zeta1, r6, pred);
    r19 = fqmul(zeta1, r7, pred);
    r20 = fqmul(zeta2, r10, pred);
    r21 = fqmul(zeta2, r11, pred);
    r22 = fqmul(zeta3, r14, pred);
    r23 = fqmul(zeta3, r15, pred);

    r2 = svsub_s16_x(pred, r0, r16);
    r3 = svsub_s16_x(pred, r1, r17);
    r6 = svsub_s16_x(pred, r4, r18);
    r7 = svsub_s16_x(pred, r5, r19);
    r10 = svsub_s16_x(pred, r8, r20);
    r11 = svsub_s16_x(pred, r9, r21);
    r14 = svsub_s16_x(pred, r12, r22);
    r15 = svsub_s16_x(pred, r13, r23);

    r0 = svadd_s16_x(pred, r0, r16);
    r1 = svadd_s16_x(pred, r1, r17);
    r4 = svadd_s16_x(pred, r4, r18);
    r5 = svadd_s16_x(pred, r5, r19);
    r8 = svadd_s16_x(pred, r8, r20);
    r9 = svadd_s16_x(pred, r9, r21);
    r12 = svadd_s16_x(pred, r12, r22);
    r13 = svadd_s16_x(pred, r13, r23);

    // Level 4

    zeta0 = svdup_s16(zetas[16+half*8]);
    zeta1 = svdup_s16(zetas[17+half*8]);
    zeta2 = svdup_s16(zetas[18+half*8]);
    zeta3 = svdup_s16(zetas[19+half*8]);
    zeta4 = svdup_s16(zetas[20+half*8]);
    zeta5 = svdup_s16(zetas[21+half*8]);
    zeta6 = svdup_s16(zetas[22+half*8]);
    zeta7 = svdup_s16(zetas[23+half*8]);

    r16 = fqmul(zeta0, r1, pred);
    r17 = fqmul(zeta1, r3, pred);
    r18 = fqmul(zeta2, r5, pred);
    r19 = fqmul(zeta3, r7, pred);
    r20 = fqmul(zeta4, r9, pred);
    r21 = fqmul(zeta5, r11, pred);
    r22 = fqmul(zeta6, r13, pred);
    r23 = fqmul(zeta7, r15, pred);

    r1 = svsub_s16_x(pred, r0, r16);
    r3 = svsub_s16_x(pred, r2, r17);
    r5 = svsub_s16_x(pred, r4, r18);
    r7 = svsub_s16_x(pred, r6, r19);
    r9 = svsub_s16_x(pred, r8, r20);
    r11 = svsub_s16_x(pred, r10, r21);
    r13 = svsub_s16_x(pred, r12, r22);
    r15 = svsub_s16_x(pred, r14, r23);

    r0 = svadd_s16_x(pred, r0, r16);
    r2 = svadd_s16_x(pred, r2, r17);
    r4 = svadd_s16_x(pred, r4, r18);
    r6 = svadd_s16_x(pred, r6, r19);
    r8 = svadd_s16_x(pred, r8, r20);
    r10 = svadd_s16_x(pred, r10, r21);
    r12 = svadd_s16_x(pred, r12, r22);
    r14 = svadd_s16_x(pred, r14, r23);

    // Level 5

    zeta0 = svdup_n_s16_x(pred_2nd_half, zetas[32+half*16]);
    zeta1 = svdup_n_s16_x(pred_2nd_half, zetas[33+half*16]);
    zeta2 = svdup_n_s16_x(pred_2nd_half, zetas[34+half*16]);
    zeta3 = svdup_n_s16_x(pred_2nd_half, zetas[35+half*16]);
    zeta4 = svdup_n_s16_x(pred_2nd_half, zetas[36+half*16]);
    zeta5 = svdup_n_s16_x(pred_2nd_half, zetas[37+half*16]);
    zeta6 = svdup_n_s16_x(pred_2nd_half, zetas[38+half*16]);
    zeta7 = svdup_n_s16_x(pred_2nd_half, zetas[39+half*16]);
    zeta8 = svdup_n_s16_x(pred_2nd_half, zetas[40+half*16]);
    zeta9 = svdup_n_s16_x(pred_2nd_half, zetas[41+half*16]);
    zeta10 = svdup_n_s16_x(pred_2nd_half, zetas[42+half*16]);
    zeta11 = svdup_n_s16_x(pred_2nd_half, zetas[43+half*16]);
    zeta12 = svdup_n_s16_x(pred_2nd_half, zetas[44+half*16]);
    zeta13 = svdup_n_s16_x(pred_2nd_half, zetas[45+half*16]);
    zeta14 = svdup_n_s16_x(pred_2nd_half, zetas[46+half*16]);
    zeta15 = svdup_n_s16_x(pred_2nd_half, zetas[47+half*16]);

    r16 = fqmul(zeta0, r0, pred_2nd_half);
    r17 = fqmul(zeta1, r1, pred_2nd_half);
    r18 = fqmul(zeta2, r2, pred_2nd_half);
    r19 = fqmul(zeta3, r3, pred_2nd_half);
    r20 = fqmul(zeta4, r4, pred_2nd_half);
    r21 = fqmul(zeta5, r5, pred_2nd_half);
    r22 = fqmul(zeta6, r6, pred_2nd_half);
    r23 = fqmul(zeta7, r7, pred_2nd_half);
    r24 = fqmul(zeta8, r8, pred_2nd_half);
    r25 = fqmul(zeta9, r9, pred_2nd_half);
    r26 = fqmul(zeta10, r10, pred_2nd_half);
    r27 = fqmul(zeta11, r11, pred_2nd_half);
    r28 = fqmul(zeta12, r12, pred_2nd_half);
    r29 = fqmul(zeta13, r13, pred_2nd_half);
    r30 = fqmul(zeta14, r14, pred_2nd_half);
    r31 = fqmul(zeta15, r15, pred_2nd_half);

    r16 = svsplice_s16(pred_2nd_half, r16, r16);
    r17 = svsplice_s16(pred_2nd_half, r17, r17);
    r18 = svsplice_s16(pred_2nd_half, r18, r18);
    r19 = svsplice_s16(pred_2nd_half, r19, r19);
    r20 = svsplice_s16(pred_2nd_half, r20, r20);
    r21 = svsplice_s16(pred_2nd_half, r21, r21);
    r22 = svsplice_s16(pred_2nd_half, r22, r22);
    r23 = svsplice_s16(pred_2nd_half, r23, r23);
    r24 = svsplice_s16(pred_2nd_half, r24, r24);
    r25 = svsplice_s16(pred_2nd_half, r25, r25);
    r26 = svsplice_s16(pred_2nd_half, r26, r26);
    r27 = svsplice_s16(pred_2nd_half, r27, r27);
    r28 = svsplice_s16(pred_2nd_half, r28, r28);
    r29 = svsplice_s16(pred_2nd_half, r29, r29);
    r30 = svsplice_s16(pred_2nd_half, r30, r30);
    r31 = svsplice_s16(pred_2nd_half, r31, r31);

    // update b
    r32 = svsub_s16_m(pred_1st_half, r0, r16);
    r33 = svsub_s16_m(pred_1st_half, r1, r17);
    r34 = svsub_s16_m(pred_1st_half, r2, r18);
    r35 = svsub_s16_m(pred_1st_half, r3, r19);
    r36 = svsub_s16_m(pred_1st_half, r4, r20);
    r37 = svsub_s16_m(pred_1st_half, r5, r21);
    r38 = svsub_s16_m(pred_1st_half, r6, r22);
    r39 = svsub_s16_m(pred_1st_half, r7, r23);
    r40 = svsub_s16_m(pred_1st_half, r8, r24);
    r41 = svsub_s16_m(pred_1st_half, r9, r25);
    r42 = svsub_s16_m(pred_1st_half, r10, r26);
    r43 = svsub_s16_m(pred_1st_half, r11, r27);
    r44 = svsub_s16_m(pred_1st_half, r12, r28);
    r45 = svsub_s16_m(pred_1st_half, r13, r29);
    r46 = svsub_s16_m(pred_1st_half, r14, r30);
    r47 = svsub_s16_m(pred_1st_half, r15, r31);

    r0 = svsplice_s16(pred_1st_half, r0, r32);
    r1 = svsplice_s16(pred_1st_half, r1, r33);
    r2 = svsplice_s16(pred_1st_half, r2, r34);
    r3 = svsplice_s16(pred_1st_half, r3, r35);
    r4 = svsplice_s16(pred_1st_half, r4, r36);
    r5 = svsplice_s16(pred_1st_half, r5, r37);
    r6 = svsplice_s16(pred_1st_half, r6, r38);
    r7 = svsplice_s16(pred_1st_half, r7, r39);
    r8 = svsplice_s16(pred_1st_half, r8, r40);
    r9 = svsplice_s16(pred_1st_half, r9, r41);
    r11 = svsplice_s16(pred_1st_half, r11, r43);
    r10 = svsplice_s16(pred_1st_half, r10, r42);
    r12 = svsplice_s16(pred_1st_half, r12, r44);
    r13 = svsplice_s16(pred_1st_half, r13, r45);
    r14 = svsplice_s16(pred_1st_half, r14, r46);
    r15 = svsplice_s16(pred_1st_half, r15, r47);

    // update a
    r0 = svadd_s16_m(pred_1st_half, r0, r16);
    r1 = svadd_s16_m(pred_1st_half, r1, r17);
    r2 = svadd_s16_m(pred_1st_half, r2, r18);
    r3 = svadd_s16_m(pred_1st_half, r3, r19);
    r4 = svadd_s16_m(pred_1st_half, r4, r20);
    r5 = svadd_s16_m(pred_1st_half, r5, r21);
    r6 = svadd_s16_m(pred_1st_half, r6, r22);
    r7 = svadd_s16_m(pred_1st_half, r7, r23);
    r8 = svadd_s16_m(pred_1st_half, r8, r24);
    r9 = svadd_s16_m(pred_1st_half, r9, r25);
    r10 = svadd_s16_m(pred_1st_half, r10, r26);
    r11 = svadd_s16_m(pred_1st_half, r11, r27);
    r12 = svadd_s16_m(pred_1st_half, r12, r28);
    r13 = svadd_s16_m(pred_1st_half, r13, r29);
    r14 = svadd_s16_m(pred_1st_half, r14, r30);
    r15 = svadd_s16_m(pred_1st_half, r15, r31);

    // Level 6

    zeta0 = svdup_n_s16_x(pred_2nd_quarter, zetas[64+half*32]);
    zeta1 = svdup_n_s16_x(pred_4th_quarter, zetas[65+half*32]);
    zeta2 = svdup_n_s16_x(pred_2nd_quarter, zetas[66+half*32]);
    zeta3 = svdup_n_s16_x(pred_4th_quarter, zetas[67+half*32]);
    zeta4 = svdup_n_s16_x(pred_2nd_quarter, zetas[68+half*32]);
    zeta5 = svdup_n_s16_x(pred_4th_quarter, zetas[69+half*32]);
    zeta6 = svdup_n_s16_x(pred_2nd_quarter, zetas[70+half*32]);
    zeta7 = svdup_n_s16_x(pred_4th_quarter, zetas[71+half*32]);
    zeta8 = svdup_n_s16_x(pred_2nd_quarter, zetas[72+half*32]);
    zeta9 = svdup_n_s16_x(pred_4th_quarter, zetas[73+half*32]);
    zeta10 = svdup_n_s16_x(pred_2nd_quarter, zetas[74+half*32]);
    zeta11 = svdup_n_s16_x(pred_4th_quarter, zetas[75+half*32]);
    zeta12 = svdup_n_s16_x(pred_2nd_quarter, zetas[76+half*32]);
    zeta13 = svdup_n_s16_x(pred_4th_quarter, zetas[77+half*32]);
    zeta14 = svdup_n_s16_x(pred_2nd_quarter, zetas[78+half*32]);
    zeta15 = svdup_n_s16_x(pred_4th_quarter, zetas[79+half*32]);
    zeta16 = svdup_n_s16_x(pred_2nd_quarter, zetas[80+half*32]);
    zeta17 = svdup_n_s16_x(pred_4th_quarter, zetas[81+half*32]);
    zeta18 = svdup_n_s16_x(pred_2nd_quarter, zetas[82+half*32]);
    zeta19 = svdup_n_s16_x(pred_4th_quarter, zetas[83+half*32]);
    zeta20 = svdup_n_s16_x(pred_2nd_quarter, zetas[84+half*32]);
    zeta21 = svdup_n_s16_x(pred_4th_quarter, zetas[85+half*32]);
    zeta22 = svdup_n_s16_x(pred_2nd_quarter, zetas[86+half*32]);
    zeta23 = svdup_n_s16_x(pred_4th_quarter, zetas[87+half*32]);
    zeta24 = svdup_n_s16_x(pred_2nd_quarter, zetas[88+half*32]);
    zeta25 = svdup_n_s16_x(pred_4th_quarter, zetas[89+half*32]);
    zeta26 = svdup_n_s16_x(pred_2nd_quarter, zetas[90+half*32]);
    zeta27 = svdup_n_s16_x(pred_4th_quarter, zetas[91+half*32]);
    zeta28 = svdup_n_s16_x(pred_2nd_quarter, zetas[92+half*32]);
    zeta29 = svdup_n_s16_x(pred_4th_quarter, zetas[93+half*32]);
    zeta30 = svdup_n_s16_x(pred_2nd_quarter, zetas[94+half*32]);
    zeta31 = svdup_n_s16_x(pred_4th_quarter, zetas[95+half*32]);

    r16 = fqmul(zeta0, r0, pred_2nd_quarter);
    r17 = fqmul(zeta1, r0, pred_4th_quarter);
    r18 = fqmul(zeta2, r1, pred_2nd_quarter);
    r19 = fqmul(zeta3, r1, pred_4th_quarter);
    r20 = fqmul(zeta4, r2, pred_2nd_quarter);
    r21 = fqmul(zeta5, r2, pred_4th_quarter);
    r22 = fqmul(zeta6, r3, pred_2nd_quarter);
    r23 = fqmul(zeta7, r3, pred_4th_quarter);
    r24 = fqmul(zeta8, r4, pred_2nd_quarter);
    r25 = fqmul(zeta9, r4, pred_4th_quarter);
    r26 = fqmul(zeta10, r5, pred_2nd_quarter);
    r27 = fqmul(zeta11, r5, pred_4th_quarter);
    r28 = fqmul(zeta12, r6, pred_2nd_quarter);
    r29 = fqmul(zeta13, r6, pred_4th_quarter);
    r30 = fqmul(zeta14, r7, pred_2nd_quarter);
    r31 = fqmul(zeta15, r7, pred_4th_quarter);
    r32 = fqmul(zeta16, r8, pred_2nd_quarter);
    r33 = fqmul(zeta17, r8, pred_4th_quarter);
    r34 = fqmul(zeta18, r9, pred_2nd_quarter);
    r35 = fqmul(zeta19, r9, pred_4th_quarter);
    r36 = fqmul(zeta20, r10, pred_2nd_quarter);
    r37 = fqmul(zeta21, r10, pred_4th_quarter);
    r38 = fqmul(zeta22, r11, pred_2nd_quarter);
    r39 = fqmul(zeta23, r11, pred_4th_quarter);
    r40 = fqmul(zeta24, r12, pred_2nd_quarter);
    r41 = fqmul(zeta25, r12, pred_4th_quarter);
    r42 = fqmul(zeta26, r13, pred_2nd_quarter);
    r43 = fqmul(zeta27, r13, pred_4th_quarter);
    r44 = fqmul(zeta28, r14, pred_2nd_quarter);
    r45 = fqmul(zeta29, r14, pred_4th_quarter);
    r46 = fqmul(zeta30, r15, pred_2nd_quarter);
    r47 = fqmul(zeta31, r15, pred_4th_quarter);

    // shift from 2nd -> 1st and 4th -> 3rd
    r16 = svsplice_s16(pred_2nd_quarter, r16, r16);
    r17 = svsplice_s16(pred_last_3_quarters, r17, r17);
    r18 = svsplice_s16(pred_2nd_quarter, r18, r18);
    r19 = svsplice_s16(pred_last_3_quarters, r19, r19);
    r20 = svsplice_s16(pred_2nd_quarter, r20, r20);
    r21 = svsplice_s16(pred_last_3_quarters, r21, r21);
    r22 = svsplice_s16(pred_2nd_quarter, r22, r22);
    r23 = svsplice_s16(pred_last_3_quarters, r23, r23);
    r24 = svsplice_s16(pred_2nd_quarter, r24, r24);
    r25 = svsplice_s16(pred_last_3_quarters, r25, r25);
    r26 = svsplice_s16(pred_2nd_quarter, r26, r26);
    r27 = svsplice_s16(pred_last_3_quarters, r27, r27);
    r28 = svsplice_s16(pred_2nd_quarter, r28, r28);
    r29 = svsplice_s16(pred_last_3_quarters, r29, r29);
    r30 = svsplice_s16(pred_2nd_quarter, r30, r30);
    r31 = svsplice_s16(pred_last_3_quarters, r31, r31);
    r32 = svsplice_s16(pred_2nd_quarter, r32, r32);
    r33 = svsplice_s16(pred_last_3_quarters, r33, r33);
    r34 = svsplice_s16(pred_2nd_quarter, r34, r34);
    r35 = svsplice_s16(pred_last_3_quarters, r35, r35);
    r36 = svsplice_s16(pred_2nd_quarter, r36, r36);
    r37 = svsplice_s16(pred_last_3_quarters, r37, r37);
    r38 = svsplice_s16(pred_2nd_quarter, r38, r38);
    r39 = svsplice_s16(pred_last_3_quarters, r39, r39);
    r40 = svsplice_s16(pred_2nd_quarter, r40, r40);
    r41 = svsplice_s16(pred_last_3_quarters, r41, r41);
    r42 = svsplice_s16(pred_2nd_quarter, r42, r42);
    r43 = svsplice_s16(pred_last_3_quarters, r43, r43);
    r44 = svsplice_s16(pred_2nd_quarter, r44, r44);
    r45 = svsplice_s16(pred_last_3_quarters, r45, r45);
    r46 = svsplice_s16(pred_2nd_quarter, r46, r46);
    r47 = svsplice_s16(pred_last_3_quarters, r47, r47);

    // update b
    r48 = svsub_s16_m(pred_1st_quarter, r0, r16);
    r49 = svsub_s16_m(pred_3rd_quarter, r0, r17);
    r50 = svsub_s16_m(pred_1st_quarter, r1, r18);
    r51 = svsub_s16_m(pred_3rd_quarter, r1, r19);
    r52 = svsub_s16_m(pred_1st_quarter, r2, r20);
    r53 = svsub_s16_m(pred_3rd_quarter, r2, r21);
    r54 = svsub_s16_m(pred_1st_quarter, r3, r22);
    r55 = svsub_s16_m(pred_3rd_quarter, r3, r23);
    r56 = svsub_s16_m(pred_1st_quarter, r4, r24);
    r57 = svsub_s16_m(pred_3rd_quarter, r4, r25);
    r58 = svsub_s16_m(pred_1st_quarter, r5, r26);
    r59 = svsub_s16_m(pred_3rd_quarter, r5, r27);
    r60 = svsub_s16_m(pred_1st_quarter, r6, r28);
    r61 = svsub_s16_m(pred_3rd_quarter, r6, r29);
    r62 = svsub_s16_m(pred_1st_quarter, r7, r30);
    r63 = svsub_s16_m(pred_3rd_quarter, r7, r31);
    r64 = svsub_s16_m(pred_1st_quarter, r8, r32);
    r65 = svsub_s16_m(pred_3rd_quarter, r8, r33);
    r66 = svsub_s16_m(pred_1st_quarter, r9, r34);
    r67 = svsub_s16_m(pred_3rd_quarter, r9, r35);
    r68 = svsub_s16_m(pred_1st_quarter, r10, r36);
    r69 = svsub_s16_m(pred_3rd_quarter, r10, r37);
    r70 = svsub_s16_m(pred_1st_quarter, r11, r38);
    r71 = svsub_s16_m(pred_3rd_quarter, r11, r39);
    r72 = svsub_s16_m(pred_1st_quarter, r12, r40);
    r73 = svsub_s16_m(pred_3rd_quarter, r12, r41);
    r74 = svsub_s16_m(pred_1st_quarter, r13, r42);
    r75 = svsub_s16_m(pred_3rd_quarter, r13, r43);
    r76 = svsub_s16_m(pred_1st_quarter, r14, r44);
    r77 = svsub_s16_m(pred_3rd_quarter, r14, r45);
    r78 = svsub_s16_m(pred_1st_quarter, r15, r46);
    r79 = svsub_s16_m(pred_3rd_quarter, r15, r47);

    // shift from 1st -> 2nd and 3rd -> 4th
    r48 = svsplice_s16(pred_1st_quarter, r48, r48);
    r49 = svsplice_s16(pred_1st_quarter, r49, r49);
    r50 = svsplice_s16(pred_1st_quarter, r50, r50);
    r51 = svsplice_s16(pred_1st_quarter, r51, r51);
    r52 = svsplice_s16(pred_1st_quarter, r52, r52);
    r53 = svsplice_s16(pred_1st_quarter, r53, r53);
    r54 = svsplice_s16(pred_1st_quarter, r54, r54);
    r55 = svsplice_s16(pred_1st_quarter, r55, r55);
    r56 = svsplice_s16(pred_1st_quarter, r56, r56);
    r57 = svsplice_s16(pred_1st_quarter, r57, r57);
    r58 = svsplice_s16(pred_1st_quarter, r58, r58);
    r59 = svsplice_s16(pred_1st_quarter, r59, r59);
    r60 = svsplice_s16(pred_1st_quarter, r60, r60);
    r61 = svsplice_s16(pred_1st_quarter, r61, r61);
    r62 = svsplice_s16(pred_1st_quarter, r62, r62);
    r63 = svsplice_s16(pred_1st_quarter, r63, r63);
    r64 = svsplice_s16(pred_1st_quarter, r64, r64);
    r65 = svsplice_s16(pred_1st_quarter, r65, r65);
    r66 = svsplice_s16(pred_1st_quarter, r66, r66);
    r67 = svsplice_s16(pred_1st_quarter, r67, r67);
    r68 = svsplice_s16(pred_1st_quarter, r68, r68);
    r69 = svsplice_s16(pred_1st_quarter, r69, r69);
    r70 = svsplice_s16(pred_1st_quarter, r70, r70);
    r71 = svsplice_s16(pred_1st_quarter, r71, r71);
    r72 = svsplice_s16(pred_1st_quarter, r72, r72);
    r73 = svsplice_s16(pred_1st_quarter, r73, r73);
    r74 = svsplice_s16(pred_1st_quarter, r74, r74);
    r75 = svsplice_s16(pred_1st_quarter, r75, r75);
    r76 = svsplice_s16(pred_1st_quarter, r76, r76);
    r77 = svsplice_s16(pred_1st_quarter, r77, r77);
    r78 = svsplice_s16(pred_1st_quarter, r78, r78);
    r79 = svsplice_s16(pred_1st_quarter, r79, r79);

    r0 = svsel_s16(pred_2nd_quarter, r48, r0);
    r0 = svsel_s16(pred_4th_quarter, r49, r0);
    r1 = svsel_s16(pred_2nd_quarter, r50, r1);
    r1 = svsel_s16(pred_4th_quarter, r51, r1);
    r2 = svsel_s16(pred_2nd_quarter, r52, r2);
    r2 = svsel_s16(pred_4th_quarter, r53, r2);
    r3 = svsel_s16(pred_2nd_quarter, r54, r3);
    r3 = svsel_s16(pred_4th_quarter, r55, r3);
    r4 = svsel_s16(pred_2nd_quarter, r56, r4);
    r4 = svsel_s16(pred_4th_quarter, r57, r4);
    r5 = svsel_s16(pred_2nd_quarter, r58, r5);
    r5 = svsel_s16(pred_4th_quarter, r59, r5);
    r6 = svsel_s16(pred_2nd_quarter, r60, r6);
    r6 = svsel_s16(pred_4th_quarter, r61, r6);
    r7 = svsel_s16(pred_2nd_quarter, r62, r7);
    r7 = svsel_s16(pred_4th_quarter, r63, r7);
    r8 = svsel_s16(pred_2nd_quarter, r64, r8);
    r8 = svsel_s16(pred_4th_quarter, r65, r8);
    r9 = svsel_s16(pred_2nd_quarter, r66, r9);
    r9 = svsel_s16(pred_4th_quarter, r67, r9);
    r10 = svsel_s16(pred_2nd_quarter, r68, r10);
    r10 = svsel_s16(pred_4th_quarter, r69, r10);
    r11 = svsel_s16(pred_2nd_quarter, r70, r11);
    r11 = svsel_s16(pred_4th_quarter, r71, r11);
    r12 = svsel_s16(pred_2nd_quarter, r72, r12);
    r12 = svsel_s16(pred_4th_quarter, r73, r12);
    r13 = svsel_s16(pred_2nd_quarter, r74, r13);
    r13 = svsel_s16(pred_4th_quarter, r75, r13);
    r14 = svsel_s16(pred_2nd_quarter, r76, r14);
    r14 = svsel_s16(pred_4th_quarter, r77, r14);
    r15 = svsel_s16(pred_2nd_quarter, r78, r15);
    r15 = svsel_s16(pred_4th_quarter, r79, r15);

    // update a
    r0 = svadd_s16_m(pred_1st_quarter, r0, r16);
    r0 = svadd_s16_m(pred_3rd_quarter, r0, r17);
    r1 = svadd_s16_m(pred_1st_quarter, r1, r18);
    r1 = svadd_s16_m(pred_3rd_quarter, r1, r19);
    r2 = svadd_s16_m(pred_1st_quarter, r2, r20);
    r2 = svadd_s16_m(pred_3rd_quarter, r2, r21);
    r3 = svadd_s16_m(pred_1st_quarter, r3, r22);
    r3 = svadd_s16_m(pred_3rd_quarter, r3, r23);
    r4 = svadd_s16_m(pred_1st_quarter, r4, r24);
    r4 = svadd_s16_m(pred_3rd_quarter, r4, r25);
    r5 = svadd_s16_m(pred_1st_quarter, r5, r26);
    r5 = svadd_s16_m(pred_3rd_quarter, r5, r27);
    r6 = svadd_s16_m(pred_1st_quarter, r6, r28);
    r6 = svadd_s16_m(pred_3rd_quarter, r6, r29);
    r7 = svadd_s16_m(pred_1st_quarter, r7, r30);
    r7 = svadd_s16_m(pred_3rd_quarter, r7, r31);
    r8 = svadd_s16_m(pred_1st_quarter, r8, r32);
    r8 = svadd_s16_m(pred_3rd_quarter, r8, r33);
    r9 = svadd_s16_m(pred_1st_quarter, r9, r34);
    r9 = svadd_s16_m(pred_3rd_quarter, r9, r35);
    r10 = svadd_s16_m(pred_1st_quarter, r10, r36);
    r10 = svadd_s16_m(pred_3rd_quarter, r10, r37);
    r11 = svadd_s16_m(pred_1st_quarter, r11, r38);
    r11 = svadd_s16_m(pred_3rd_quarter, r11, r39);
    r12 = svadd_s16_m(pred_1st_quarter, r12, r40);
    r12 = svadd_s16_m(pred_3rd_quarter, r12, r41);
    r13 = svadd_s16_m(pred_1st_quarter, r13, r42);
    r13 = svadd_s16_m(pred_3rd_quarter, r13, r43);
    r14 = svadd_s16_m(pred_1st_quarter, r14, r44);
    r14 = svadd_s16_m(pred_3rd_quarter, r14, r45);
    r15 = svadd_s16_m(pred_1st_quarter, r15, r46);
    r15 = svadd_s16_m(pred_3rd_quarter, r15, r47);

    // Store result
    svst1_s16(pred, &r[half*128+0*8],r0);
    svst1_s16(pred, &r[half*128+1*8],r1);
    svst1_s16(pred, &r[half*128+2*8],r2);
    svst1_s16(pred, &r[half*128+3*8],r3);
    svst1_s16(pred, &r[half*128+4*8],r4);
    svst1_s16(pred, &r[half*128+5*8],r5);
    svst1_s16(pred, &r[half*128+6*8],r6);
    svst1_s16(pred, &r[half*128+7*8],r7);
    svst1_s16(pred, &r[half*128+64+0*8],r8);
    svst1_s16(pred, &r[half*128+64+1*8],r9);
    svst1_s16(pred, &r[half*128+64+2*8],r10);
    svst1_s16(pred, &r[half*128+64+3*8],r11);
    svst1_s16(pred, &r[half*128+64+4*8],r12);
    svst1_s16(pred, &r[half*128+64+5*8],r13);
    svst1_s16(pred, &r[half*128+64+6*8],r14);
    svst1_s16(pred, &r[half*128+64+7*8],r15);
}

void ntt_sve128(int16_t [256], const int16_t [128]);
