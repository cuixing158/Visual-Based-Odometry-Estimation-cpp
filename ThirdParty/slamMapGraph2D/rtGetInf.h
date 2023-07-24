///
/// @file           : rtGetInf.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef RTGETINF_H
#define RTGETINF_H

/// @include file    : Include Files
#include "rtwtypes.h"

#ifdef __cplusplus
extern "C" {
#endif

extern real_T rtGetInf(void);
extern real32_T rtGetInfF(void);
extern real_T rtGetMinusInf(void);
extern real32_T rtGetMinusInfF(void);

#ifdef __cplusplus
}
#endif
#endif
///
/// File trailer for rtGetInf.h
///
/// [EOF]
///
