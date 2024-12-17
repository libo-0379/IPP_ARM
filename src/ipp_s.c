#include <stdio.h>
#include <math.h>
#include <arm_neon.h>
#include "../include/ipp_s.h"

IppStatus ippsExp_32f_A24_s(const Ipp32f *pSrc, Ipp32f *dst, Ipp32s len)
{
    // 初始化第一个值
    int n = 0;
    float32_t prior = 1.0; //前一项数值初始化为第一项的数值
    float32_t sum = prior; //求和保存结果
    Ipp32f x = *pSrc;
    while (1)
    {
        float32_t cur = prior * x / ++n;   //当前加项的数值 = 前项*(src[i]/n)
        // printf("cur %f\n",cur);
        if (fabs(cur)-EPSILON < 0) //如果小于精度值停止计算
            break;
        sum += cur;
        prior = cur;
    }
    // printf("n %d\n",n);
    *dst = sum;
    return ippStsNoErr;
}

IppStatus ippsCos_32f_A24_s(const Ipp32f *pSrc,Ipp32f *dst,Ipp32s len)
{
    // 初始化第一个值
    int n = 0;
    float32_t prior = 1.0; //前一项数值初始化为第一项的数值
    float32_t sum = prior; //求和保存结果
    Ipp64f x = *pSrc;
    //将x转换到 [-2π,2π]
    Ipp32s k = x/IPP_2PI;
    x -= k*IPP_2PI;
    while (1)
    {
        int32_t temp1 = 2* ++n;
        float32_t cur = prior * (x * x) /(temp1-- * temp1) ;   //当前加项的数值 = 前项*(x^2/(2n*(2n-1)))
        prior = cur;
        if (fabs(cur) - EPSILON <= 0) //如果小于精度值停止计算
            break;
        sum += (n%2==0?1:-1)*cur;   //符号位
    }
    // printf("n %d\n",n);
    *dst = sum;
    return ippStsNoErr;
}

IppStatus ippsSin_32f_A24_s(const Ipp32f *pSrc, Ipp32f *dst, Ipp32s len)
{
    // 初始化第一个值
    int n = 0;
    Ipp64f x = *pSrc;
    //将x转换到 [-2π,2π]
    Ipp32s k = x/IPP_2PI;
    x -= k*IPP_2PI;
    float32_t prior = x; //前一项数值初始化为第一项的数值
    float32_t sum = prior; //求和保存结果
    printf("x %f\n",x);
    while (1)
    {
        int32_t temp1 = 2* ++n;
        float32_t cur = prior * (x * x) /(temp1++ * temp1) ;   //当前加项的数值 = 前项*(x^2/(2n*(2n+1)))
        prior = cur;
        // printf("cur %f\n",cur);
        if (fabs(cur) - EPSILON <= 0) //如果小于精度值停止计算
            break;
        sum += (n%2==0?1:-1)*cur;   //符号位
    }
    printf("n %d\n",n);
    *dst = sum;
    return ippStsNoErr;
}

// x[n] = magn * cth(2π* rFreq*n + phase), n = 0, 1, 2,...
IppStatus ippsTone_32f_s(Ipp32f *pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f *pPhase, 
IppHintAlgorithm hint)
{
    if(rFreq<0.0 || (rFreq-0.5)>EPSILON)
        return ippStsToneFreqErr;
    if(*pPhase<0.0 || (*pPhase-IPP_2PI)>=EPSILON)
        return ippStsTonePhaseErr;
    if(magn<=0.0)
        return ippStsToneMagnErr;
    Ipp32f phase = *pPhase;
    for(int i=0;i<len;i++)
    {
        float32_t result = 0.0;
        Ipp32f x = IPP_2PI*rFreq*i + phase;
        ippsCos_32f_A24_s(&x,&result,1);
        *(pDst+i) = result * magn;  //类型强转
    }
    return ippStsNoErr;
}

// x[n] = magn * [cth(2π* rFreq*n + phase) + j * sth(2π* rFreq*n + phase)], n = 0, 1, 2,...
IppStatus ippsTriangle_32fc_s(Ipp32fc* pDst, int len, Ipp32f magn, Ipp32f rFreq, 
Ipp32f asym, Ipp32f* pPhase)
{
    if(rFreq<0.0 || (rFreq-0.5)>EPSILON)//[0.0, 0.5)
        return ippStsToneFreqErr;
    if(*pPhase<0.0 || (*pPhase-IPP_2PI)>=EPSILON)//[0.0, 2π)
        return ippStsTonePhaseErr;
    if(magn<=0.0)
        return ippStsToneMagnErr;
    if((asym+IPP_PI)<EPSILON || (asym-IPP_PI) >=EPSILON )   //[-π，π)
        return ippStsTrnglAsymErr;

    Ipp32f H = IPP_PI+asym;
    for(int i=0;i<len;i++)
    {
        Ipp32f x = IPP_2PI*rFreq*i + *pPhase;
        printf("x %f\n",x);
        //因为x>0 可以转换到 [0,2π]
        Ipp32s k = x/IPP_2PI;
        x -= k*IPP_2PI;
        Ipp32fc result;
        //计算实部
        if(x>=0 && (H-x)>=EPSILON)
             result.re = -2*x/H + 1.0;  // -2*x/H + 1 
        else    
            result.re = (2*(x-IPP_PI)-H)/(IPP_2PI-H);// (2(x-PI)-H)/(2PI-H)
        //计算虚部
        if(x>=0.0 && ((IPP_2PI-H)/2)-x >= EPSILON)  
            result.im = 2*x/(IPP_2PI-H);    // 2*x/(2PI-H)
        else if(x-(IPP_2PI-H)/2 >=EPSILON && (IPP_2PI+H)/2-x>=EPSILON)
            result.im = -2*(x-IPP_PI)/H;    //-2*(x-PI)/H
        else
            result.im = 2*(x-IPP_2PI)/(IPP_2PI-H);  //2(x-2PI)/(2PI-H)
        result.re *= magn;
        result.im *= magn;
        *(pDst+i) = result;
    }
    return ippStsNoErr;
}

//r[i]=op1[i]/op2[i]
float32x4_t vdivq_f32_neon(float32x4_t op1,float32x4_t op2)
{
    float32x4_t res1 = {0};
    float32x4_t res2 = {0};
    float32x4_t res3 = {0};
    res1 = vrecpeq_f32(op2);
    res2 = vmulq_f32(vrecpsq_f32(op2,res1),res1);
    res3 = vmulq_f32(vrecpsq_f32(op2,res2),res2);
    return vmulq_f32(res3,op1);
}
