#include <arm_sve.h>
#include <stdio.h>

int main() {
    svint16_t a, b;
    int16_t arr[4]  = {1,2,3,4};
    svbool_t full = svptrue_b16();

    svbool_t quarter = svwhilelt_b16(0,1);
    svbool_t half = svwhilelt_b16(0,2);
    svbool_t half2 = svnot_b_z(full, half);
    svbool_t quarter1 = sveor_b_z(full, quarter, half);
    a = svld1_s16(quarter1,&arr[0]);
    b = svld1_s16(full,&arr[0]);
    a = svadd_s16_m(half2,b,a);
    // a = svsplice_s16(quarter1,a,a);
    svst1(full,&arr[0],a);

    for(unsigned int i=0; i < 4; i++) {
        printf("arr[%d] = %d\n",i,arr[i]);
    }
}
