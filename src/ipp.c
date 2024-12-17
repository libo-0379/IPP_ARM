#include <stdio.h>
#include <math.h>
#include "../include/ipp.h"

IppStatus ippsExp_32f_A24(const Ipp32f *pSrc, Ipp32f *dst, Ipp32s len)
{
    for(int i=0;i<len/4*4;i+=4)
    {   
        //每个元素泰勒展开项计数，初始化为0
        float32x4_t vn = vdupq_n_f32(0);
        float32x4_t vx = vld1q_f32(pSrc+i);
        //记录当前项的展开值，初始化为第一项的值
        float32x4_t vprior = vdupq_n_f32(1.0);
        //保存每个元素展开的求和
        float32x4_t vsum = vprior;
        while(1)
        {
            //计算当前展开项
            // 1. 计算 prior * x / ++n
            vn = vaddq_f32(vn,vdupq_n_f32(1));
            // 自己实现的 vdivq_f32_neon 对 -20.055 计算误差太大，改用系统提供接口 vdivq_f32
            // float32x4_t vcur = vdivq_f32_neon(vmulq_f32(vprior,vx),vn);
            float32x4_t vcur = vdivq_f32(vmulq_f32(vprior,vx),vn);
            // 2. 将当前项赋值为前项            
            vprior = vcur;
            // printf("0 %f\n",vgetq_lane_f32(vcur,0));
            // 3. 如果所有元素的最大展开项小于精度控制，不再计算
            if(vmaxvq_f32(vabsq_f32(vcur)) - EPSILON < 0)
                break;
            // 4. 将当前展开项加和
            vsum = vaddq_f32(vsum,vcur);
        }
        // printf("n %f\n",vgetq_lane_f32(vn,0));
        vst1q_f32(dst+i,vsum);
    }

    //处理尾部数据
    for(int i=len/4*4;i<len;i++)
    {
        // 初始化第一个值
        int n = 0;
        double prior = 1.0; //前一项数值初始化为第一项的数值
        double sum = prior; //求和保存结果
        Ipp32f x = *(pSrc+i);
        while (1)
        {
            double cur = prior * x / ++n;   //当前加项的数值 = 前项*(src[i]/n)
            prior = cur;
            // printf("%f\n",cur);
            if (fabs(cur)-EPSILON <= 0) //如果小于精度值停止计算
                break;
            sum += cur;
        }
        // printf("n %d\n",n);
        *(dst+i) = sum;
    }
    return ippStsNoErr;
}

IppStatus ippsCos_32f_A24(const Ipp32f *pSrc, Ipp32f *dst, Ipp32s len)
{
    float32x4_t v2PI = vdupq_n_f32(IPP_2PI);
    for(int i=0;i<len/4*4;i+=4)
    {
        //每个元素泰勒展开项计数，初始化为0
        float32x4_t vn = vdupq_n_f32(0);
        float32x4_t vx = vld1q_f32(pSrc+i);
        //记录当前项的展开值，初始化为第一项的值
        float32x4_t vprior = vdupq_n_f32(1.0);
        //保存每个元素展开的求和
        float32x4_t vsum = vprior;
        //将x转换到 [-2π,2π]
        //类型转换向下取整,注意这里是浮点数除法
        float32x4_t vk = vcvtq_f32_s32(vcvtq_s32_f32(vdivq_f32_neon(vx,v2PI)));   
        vx = vfmsq_f32(vx,v2PI,vk);
        while(1)
        {
            //计算当前展开项
            // 1. 计算 prior * (x^2) / (2(++n)*(2n-1))
            vn = vaddq_f32(vn,vdupq_n_f32(1));
            float32x4_t vcur = vmulq_n_f32(vn,2);
            vcur = vdivq_f32_neon(vmulq_f32(vprior,vmulq_f32(vx,vx)),vmulq_f32(vcur,vsubq_f32(vcur,vdupq_n_f32(1.0))));
            //将当前项赋值为前项,不带符号          
            vprior = vcur;
            // 2. 添加正负号
            float32_t n = vgetq_lane_f32(vn,0);
            vcur = vmulq_n_f32(vcur,(int)n%2==0?1:-1);
            // 3. 如果所有元素的展开项绝对值最大项小于精度控制，不再计算
            if(vmaxvq_f32(vabsq_f32(vcur)) - EPSILON <= 0)
                break;
            // 4. 将当前展开项加和
            vsum = vaddq_f32(vsum,vcur);
        }
        printf("n %f\n",vgetq_lane_f32(vn,0));
        vst1q_f32(dst+i,vsum);
    }
    //处理尾部数据
    for(int i=len/4*4;i<len;i++)
    {
        // 初始化第一个值
        int n = 0;
        double prior = 1.0; //前一项数值初始化为第一项的数值
        double sum = prior; //求和保存结果
        Ipp64f x = *(pSrc+i);
        //对x转换到 [-2π,2π]
        Ipp32s k = x/IPP_2PI;
        x -= k*IPP_2PI;
        while (1)
        {
            int32_t temp1 = 2* ++n;
            double cur = prior * (x * x) /(temp1-- * temp1) ;   //当前加项的数值 = 前项*(x^2/(2n*(2n-1)))
            prior = cur;
            if (fabs(cur) - EPSILON <= 0) //如果小于精度值停止计算
                break;
            sum += (n%2==0?1:-1)*cur;   //符号位
        }
        // printf("n %d\n",n);
        *(dst+i) = sum;
    }
    return ippStsNoErr;
}

IppStatus ippsSin_32f_A24(const Ipp32f *pSrc, Ipp32f *dst, Ipp32s len)
{    
    float32x4_t v2PI = vdupq_n_f32(IPP_2PI);
    for(int i=0;i<len/4*4;i+=4)
    {
        //每个元素泰勒展开项计数，初始化为0
        float32x4_t vn = vdupq_n_f32(0);
        float32x4_t vx = vld1q_f32(pSrc+i);
        //将x转换到 [-2π,2π]
        //类型转换向下取整,注意这里是浮点数除法
        float32x4_t vk = vcvtq_f32_s32(vcvtq_s32_f32(vdivq_f32_neon(vx,v2PI)));   
        vx = vfmsq_f32(vx,v2PI,vk);
        // !! 注意这里易错点，应在区间转换后再赋值给 vprior 
        //记录当前项的展开值，初始化为第一项的值
        float32x4_t vprior = vx;
        //保存每个元素展开的求和
        float32x4_t vsum = vprior;
        printf("x0 %f\n",vgetq_lane_f32(vx,0));
        while(1)
        {
            // 1. 计算 prior * (x^2) / (2(++n)*(2n+1))
            vn = vaddq_f32(vn,vdupq_n_f32(1));
            float32x4_t vcur = vmulq_n_f32(vn,2);
            vcur = vdivq_f32_neon(vmulq_f32(vprior,vmulq_f32(vx,vx)),vmulq_f32(vcur,vaddq_f32(vcur,vdupq_n_f32(1))));
            // printf("cur0 %f\n",vgetq_lane_f32(vcur,0));
            //将当前项赋值为前项,不带符号          
            vprior = vcur;
            // 2. 添加正负号
            float32_t n = vgetq_lane_f32(vn,0);
            vcur = vmulq_n_f32(vcur,(int)n%2==0?1:-1);
            // 3.如果所有元素的展开项绝对值最大项小于精度控制，不再计算
            if(vmaxvq_f32(vabsq_f32(vcur)) - EPSILON <= 0)
                break;
            // 4.将当前展开项加和
            vsum = vaddq_f32(vsum,vcur);
        }
        printf("n %f\n",vgetq_lane_f32(vn,0));
        vst1q_f32(dst+i,vsum);
    }
    //处理尾部数据
    for(int i=len/4*4;i<len;i++)
    {
        // 初始化第一个值
        int n = 0;
        Ipp64f x = *(pSrc+i);
        //将x转换到 [-2π,2π]
        Ipp32s k = x/IPP_2PI;
        x -= k*IPP_2PI;
        double prior = x; //前一项数值初始化为第一项的数值
        double sum = prior; //求和保存结果
        while (1)
        {
            int32_t temp1 = 2* ++n;
            double cur = prior * (x * x) /(temp1++ * temp1) ;   //当前加项的数值 = 前项*(x^2/(2n*(2n+1)))
            prior = cur;
            if (fabs(cur) - EPSILON <= 0) //如果小于精度值停止计算
                break;
            sum += (n%2==0?1:-1)*cur;   //符号位
        }
        // printf("n %d\n",n);
        *(dst+i) = sum;
    }
    return ippStsNoErr;
}

IppStatus ippsTone_32f(Ipp32f* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f*
pPhase, IppHintAlgorithm hint)
{
    if(rFreq<0.0 || (rFreq-0.5)>EPSILON)
        return ippStsToneFreqErr;
    if(*pPhase<0.0 || (*pPhase-IPP_2PI)>=EPSILON)
        return ippStsTonePhaseErr;
    if(magn<=0.0)
        return ippStsToneMagnErr;
    
    float32x4_t vphase = vdupq_n_f32(*pPhase);
    float32x4_t vrFreq= vdupq_n_f32(rFreq);
    for(int i=0;i<len/4*4;i+=4)
    {
        //保存结果
        float32_t result[4] = {0.0};
        float32x4_t vresult = {0.0};
        // IPP_2PI*rFreq*i + phase
        float32_t x[4] = {i,i+1,i+2,i+3};
        float32x4_t vx = vld1q_f32(x);
        vx = vaddq_f32(vphase,vmulq_n_f32(vmulq_f32(vrFreq,vx),IPP_2PI));
        vst1q_f32(x,vx);
        ippsCos_32f_A24(x,result,4);
        vresult = vld1q_f32(result);
        vst1q_f32(pDst+i,vmulq_n_f32(vresult,magn));
    }
    for(int i=len/4*4;i<len;i++)
    {
        float32_t result = 0.0;
        Ipp32f x = IPP_2PI*rFreq*i + (*pPhase);
        ippsCos_32f_A24_s(&x,&result,1);
        *(pDst+i) = result * magn;  //类型强转
    }
    return ippStsNoErr;
}

IppStatus ippsTriangle_32fc(Ipp32fc* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f
asym, Ipp32f* pPhase)
{
    if(rFreq<0.0 || (rFreq-0.5)>EPSILON)
        return ippStsToneFreqErr;
    if(*pPhase<0.0 || (*pPhase-IPP_2PI)>=EPSILON)
        return ippStsTonePhaseErr;
    if(magn<=0.0)
        return ippStsToneMagnErr;
    
    float32x4_t vH = vdupq_n_f32(IPP_PI+asym);
    printf("h %f\n",IPP_PI+asym);

    float32x4_t vPI = vdupq_n_f32(IPP_PI);
    float32x4_t v2PI = vdupq_n_f32(IPP_2PI);
    float32x4_t v2PI_H = vsubq_f32(v2PI,vH);
    float32x4_t vphase = vdupq_n_f32(*pPhase);
    float32x4_t vrFreq= vdupq_n_f32(rFreq);
    for(int i=0;i<len/4*4;i+=4)
    {
        float32_t x[4] = {i,i+1,i+2,i+3};
        float32x4_t vx = vld1q_f32(x);
        // 计算位置 x = 2*PI*rFreq*n + phase
        vx = vaddq_f32(vphase,vmulq_n_f32(vmulq_f32(vrFreq,vx),IPP_2PI));
        //因为x>0 可以转换到 [0,2π]
        //类型转换向下取整,注意这里是浮点数除法
        float32x4_t vk = vcvtq_f32_s32(vcvtq_s32_f32(vdivq_f32_neon(vx,v2PI)));   
        vx = vfmsq_f32(vx,v2PI,vk);

        printf("x %f\n",vgetq_lane_f32(vx,0));
        
        vst1q_f32(x,vx);
        // for(int i=0;i<4;i++)
        //     printf("x %f\n",x[i]);
        //保存结果包含实部虚部
        float32x4x2_t vresult;
        //实部
        float32x4_t vresultRE;
        float32x4_t vresultRE1;
        float32x4_t vresultRE2;
        //虚部
        float32x4_t vresultIM;
        float32x4_t vresultIM_1;
        float32x4_t vresultIM1;
        float32x4_t vresultIM2;
        float32x4_t vresultIM3;
        //对实部和虚部所有区间都计算最后根据 x 的范围取舍
        //计算实部
        //[0,H] -2*x/H + 1 
        vresultRE1 = vaddq_f32(vdupq_n_f32(1.0),vdivq_f32_neon(vmulq_n_f32(vx,-2),vH));
        //[H,2π] (2(x-PI)-H)/(2PI-H)
        vresultRE2 = vdivq_f32_neon(vsubq_f32(vmulq_n_f32(vsubq_f32(vx,vPI),2),vH),v2PI_H);

        //计算虚部
        //[0,π-H/2] 2*x/(2PI-H)
        vresultIM1 = vdivq_f32_neon(vmulq_n_f32(vx,2),v2PI_H);
        //[π-H/2,π+H/2] -2*(x-PI)/H
        vresultIM2 = vdivq_f32_neon(vmulq_n_f32(vsubq_f32(vx,vPI),-2),vH);
        //[π+H/2,2π] 2(x-2PI)/(2PI-H)
        vresultIM3 = vdivq_f32_neon(vmulq_n_f32(vsubq_f32(vx,v2PI),2),v2PI_H);

        //实部取舍
        uint32x4_t vcompare = vcltq_f32(vx,vH);
    
        printf("c %u\n",vgetq_lane_u32(vcompare,0));

        vresultRE = vbslq_f32(vcompare,vresultRE1,vresultRE2);
        vresultRE = vmulq_n_f32(vresultRE,magn);
        //虚部取舍
        float32x4_t v2_H = vdivq_f32_neon(vH,vdupq_n_f32(2));
        // x< π-H/2 取 vresultIM1,否则取 vresultIM2
        vcompare = vcltq_f32(vx,vsubq_f32(vPI,v2_H));
        vresultIM_1 = vbslq_f32(vcompare,vresultIM1,vresultIM2);
        // x> π+H/2 取 vresultIM3,否则取 vresultIM
        vcompare = vcgtq_f32(vx,vaddq_f32(vPI,v2_H));
        vresultIM = vbslq_f32(vcompare,vresultIM3,vresultIM_1);
        vresultIM = vmulq_n_f32(vresultIM,magn);

        vresult.val[0] = vresultRE;
        vresult.val[1] = vresultIM;
        printf("re %f\n",vgetq_lane_f32(vresultRE,0));
        printf("im %f\n",vgetq_lane_f32(vresultIM,0));
        vst2q_f32((float32_t*)(pDst+i),vresult);
    }
    Ipp32f H = IPP_PI+asym;
    for(int i=len/4*4;i<len;i++)
    {
        Ipp32f x = IPP_2PI*rFreq*i + *pPhase;
        //x 需要转换到 [0,2π]
        Ipp32s k = x/IPP_2PI;
        x -= k*IPP_2PI;
        Ipp32fc result;
        // printf("x %f\n",x);
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
