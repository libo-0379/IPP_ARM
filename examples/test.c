#include <stdio.h>
// #include "../include/ipp.h"
#include "../include/ipp.h"
#include "../include/ipp_s.h"

#define N 5

int main()
{
    Ipp32f src[N] = {-20.055, -5.456, 7.835, 11.89, 18.76};
    Ipp32f dst[N] ={0.0};
    for(int i=0;i<N;i++)
    {
        ippsExp_32f_A24_s(&src[i],&dst[i],1);
        printf("result %f\n",dst[i]);
    }
    printf("======\n");
    ippsExp_32f_A24(src,dst,N);
    for(int i=0;i<N;i++)
        printf("result %f\n",dst[i]);
    printf("======\n");

    for(int i=0;i<N;i++)
    {
        ippsCos_32f_A24_s(&src[i],&dst[i],1);
        printf("result %f\n",dst[i]);
    }
    printf("======\n");
    ippsCos_32f_A24(src,dst,N);
    for(int i=0;i<N;i++)
        printf("result %f\n",dst[i]);
    printf("======\n");

    for(int i=0;i<N;i++)
    {
        ippsSin_32f_A24_s(&src[i],&dst[i],1);
        printf("result %f\n",dst[i]);
    }
    printf("======\n");
    ippsSin_32f_A24(src,dst,N);
    for(int i=0;i<N;i++)
        printf("result %f\n",dst[i]);
    printf("======\n");

    Ipp32f phase = 5.3865;
    // ippsTone_32f_s(dst,N,3,0.2456,&phase,ippAlgHintFast);
    // for(int i=0;i<N;i++)
    //     printf("result %f\n",dst[i]);
    // printf("======\n");
    // ippsTone_32f(dst,N,3,0.2456,&phase,ippAlgHintFast);
    // for(int i=0;i<N;i++)
    //     printf("result %f\n",dst[i]);
    // printf("======\n");



    // Ipp32fc dst1[N];
    // ippsTriangle_32fc_s(dst1,N,3,0.2456,1.47,&phase);
    // for(int i=0;i<N;i++)
    //     printf("result %f %f\n",dst1[i].re,dst1[i].im);
    // printf("======\n");

    // ippsTriangle_32fc(dst1,N,3,0.2456,1.47,&phase);
    // for(int i=0;i<N;i++)
    //     printf("result %f %f\n",dst1[i].re,dst1[i].im);
    // printf("======\n");

    return 0;
}
