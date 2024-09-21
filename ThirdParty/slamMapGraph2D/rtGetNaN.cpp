///
/// @file           : rtGetNaN.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

// Abstract:
//       MATLAB for code generation function to initialize non-finite, NaN
/// @include file    : Include Files
#include "rtGetNaN.h"
#include "rt_nonfinite.h"

// Function: rtGetNaN
// ======================================================================
//  Abstract:
// Initialize rtNaN needed by the generated code.
//  NaN is initialized as non-signaling. Assumes IEEE.
real_T rtGetNaN(void) {
    return rtNaN;
}

// Function: rtGetNaNF
// =====================================================================
//  Abstract:
//  Initialize rtNaNF needed by the generated code.
//  NaN is initialized as non-signaling. Assumes IEEE
real32_T rtGetNaNF(void) {
    return rtNaNF;
}

///
/// File trailer for rtGetNaN.cpp
///
/// [EOF]
///
