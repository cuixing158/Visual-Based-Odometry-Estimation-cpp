///
/// @file           : xswap.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "xswap.h"
#include "rt_nonfinite.h"

/// Function Definitions
///
/// @fn             : xswap
/// @brief          :
/// @param          : double x[9]
///                   int ix0
///                   int iy0
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
void xswap(double x[9], int ix0, int iy0) {
    double temp;
    temp = x[ix0 - 1];
    x[ix0 - 1] = x[iy0 - 1];
    x[iy0 - 1] = temp;
    temp = x[ix0];
    x[ix0] = x[iy0];
    x[iy0] = temp;
    temp = x[ix0 + 1];
    x[ix0 + 1] = x[iy0 + 1];
    x[iy0 + 1] = temp;
}

}  // namespace blas
}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for xswap.cpp
///
/// [EOF]
///
