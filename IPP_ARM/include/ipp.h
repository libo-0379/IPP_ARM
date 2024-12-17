#ifndef IPP_H
#define IPP_H

#include <arm_neon.h>
#include "ipp_s.h"

//page 392
extern IppStatus ippsExp_32f_A24(const Ipp32f* pSrc,Ipp32f* dst,Ipp32s len);
//page 400
extern IppStatus ippsCos_32f_A24(const Ipp32f* pSrc,Ipp32f* dst,Ipp32s len);
extern IppStatus ippsSin_32f_A24(const Ipp32f* pSrc,Ipp32f* dst,Ipp32s len);
//page 54
extern IppStatus ippsTone_32f(Ipp32f* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f*
pPhase, IppHintAlgorithm hint);
extern IppStatus ippsTriangle_32fc(Ipp32fc* pDst, int len, Ipp32f magn, Ipp32f rFreq, Ipp32f
asym, Ipp32f* pPhase);





#endif