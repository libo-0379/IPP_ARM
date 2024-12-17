/*
串行版本
*/
#ifndef IPP_S_H
#define IPP_S_H

#include <arm_neon.h>


#define EPSILON 1e-6
#define IPP_PI    ( 3.14159265358979323846 )  /* ANSI C does not support M_PI */
#define IPP_2PI   ( 6.28318530717958647692 )  /* 2*pi                         */


#define ippStsToneMagnErr                  -46 /* Tone magnitude is less than, or equal to zero. */
#define ippStsToneFreqErr                  -45 /* Tone frequency is negative, or greater than, or equal to 0.5. */
#define ippStsTonePhaseErr                 -44 /* Tone phase is negative, or greater than, or equal to 2*PI. */
#define ippStsTrnglMagnErr                 -43 /* Triangle magnitude is less than, or equal to zero. */
#define ippStsTrnglFreqErr                 -42 /* Triangle frequency is negative, or greater than, or equal to 0.5. */
#define ippStsTrnglPhaseErr                -41 /* Triangle phase is negative, or greater than, or equal to 2*PI. */
#define ippStsTrnglAsymErr                 -40 /* Triangle asymmetry is less than -PI, or greater than, or equal to PI. */
#define ippStsNullPtrErr                    -8 /* Null pointer error. */
#define ippStsSizeErr                       -6 /* Incorrect value for data size. */
#define ippStsNoErr                          0 /* No errors. */
#define ippStsOverflow                      12 /* Overflow in the operation. */
#define ippStsUnderflow                     17 /* Underflow in the operation. */


typedef enum {
 ippAlgHintNone,
 ippAlgHintFast,
 ippAlgHintAccurate
} IppHintAlgorithm;


#define CHECK_NULL(ptr) if(!ptr)\return ippStsNullPtrErr;

typedef float32_t Ipp32f;
typedef float64_t Ipp64f;
typedef int32_t Ipp32s;
typedef int16_t Ipp16s;
typedef int32_t IppStatus;
typedef struct {
    Ipp32f  re;
    Ipp32f  im;
} Ipp32fc;



//page 392;A24保证小数点后6位十进制精度
extern IppStatus ippsExp_32f_A24_s(const Ipp32f* pSrc,Ipp32f* dst,Ipp32s len);
//page 400
extern IppStatus ippsCos_32f_A24_s(const Ipp32f* pSrc,Ipp32f* dst,Ipp32s len);
extern IppStatus ippsSin_32f_A24_s(const Ipp32f* pSrc,Ipp32f* dst,Ipp32s len);
//page 54
extern IppStatus ippsTone_16s_s(Ipp16s* pDst, int len, Ipp16s magn, Ipp32f rFreq, Ipp32f*
pPhase, IppHintAlgorithm hint);
extern IppStatus ippsTone_32f_s(Ipp32f* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f*
pPhase, IppHintAlgorithm hint);
extern IppStatus ippsTriangle_32fc_s(Ipp32fc* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f
asym, Ipp32f* pPhase);


//float r[i]=op1[i]/op2[i]
extern float32x4_t vdivq_f32_neon(float32x4_t op1,float32x4_t op2);
//int r[i]=op1[i]/op2[i]
// extern int32x4_t vdivq_s32_neon(int32x4_t op1,int32x4_t op2);

/*
1. 三角函数自变量 x 未转换到 [0,2π] 区间，如果 x 过大计算结果可能超出数值范围
*/




#endif