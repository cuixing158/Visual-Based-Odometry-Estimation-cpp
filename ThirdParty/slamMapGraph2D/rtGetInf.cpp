///
/// @file           : rtGetInf.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

// Abstract:
//       MATLAB for code generation function to initialize non-finite, Inf and
//       MinusInf
/// @include file    : Include Files
#include "rtGetInf.h"
#include "rt_nonfinite.h"

// Function: rtGetInf
// ==================================================================
//  Abstract:
//  Initialize rtInf needed by the generated code.
real_T rtGetInf(void) {
    return rtInf;
}

// Function: rtGetInfF
// =================================================================
//  Abstract:
//  Initialize rtInfF needed by the generated code.
real32_T rtGetInfF(void) {
    return rtInfF;
}

// Function: rtGetMinusInf
// =============================================================
//  Abstract:
//  Initialize rtMinusInf needed by the generated code.
real_T rtGetMinusInf(void) {
    return rtMinusInf;
}

// Function: rtGetMinusInfF
// ============================================================
//  Abstract:
//  Initialize rtMinusInfF needed by the generated code.
real32_T rtGetMinusInfF(void) {
    return rtMinusInfF;
}

///
/// File trailer for rtGetInf.cpp
///
/// [EOF]
///
