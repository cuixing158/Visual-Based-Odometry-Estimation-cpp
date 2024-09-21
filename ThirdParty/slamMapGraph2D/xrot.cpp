///
/// @file           : xrot.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "xrot.h"
#include "rt_nonfinite.h"

/// Function Definitions
///
/// @fn             : xrot
/// @brief          :
/// @param          : double x[9]
///                   int ix0
///                   int iy0
///                   double c
///                   double s
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
void xrot(double x[9], int ix0, int iy0, double c, double s) {
    double temp;
    double temp_tmp;
    temp = x[iy0 - 1];
    temp_tmp = x[ix0 - 1];
    x[iy0 - 1] = c * temp - s * temp_tmp;
    x[ix0 - 1] = c * temp_tmp + s * temp;
    temp = c * x[ix0] + s * x[iy0];
    x[iy0] = c * x[iy0] - s * x[ix0];
    x[ix0] = temp;
    temp = x[iy0 + 1];
    temp_tmp = x[ix0 + 1];
    x[iy0 + 1] = c * temp - s * temp_tmp;
    x[ix0 + 1] = c * temp_tmp + s * temp;
}

}  // namespace blas
}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for xrot.cpp
///
/// [EOF]
///
