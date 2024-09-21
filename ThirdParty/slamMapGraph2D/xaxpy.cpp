///
/// @file           : xaxpy.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "xaxpy.h"
#include "rt_nonfinite.h"

/// Function Definitions
///
/// @fn             : xaxpy
/// @brief          :
/// @param          : double a
///                   const double x[9]
///                   int ix0
///                   double y[3]
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
void xaxpy(double a, const double x[9], int ix0, double y[3]) {
    if (!(a == 0.0)) {
        for (int k{0}; k < 2; k++) {
            y[k + 1] += a * x[(ix0 + k) - 1];
        }
    }
}

///
/// @fn             : xaxpy
/// @brief          :
/// @param          : double a
///                   const double x[3]
///                   double y[9]
///                   int iy0
/// @return         : void
///
void xaxpy(double a, const double x[3], double y[9], int iy0) {
    if (!(a == 0.0)) {
        for (int k{0}; k < 2; k++) {
            int i;
            i = (iy0 + k) - 1;
            y[i] += a * x[k + 1];
        }
    }
}

///
/// @fn             : xaxpy
/// @brief          :
/// @param          : int n
///                   double a
///                   int ix0
///                   double y[9]
///                   int iy0
/// @return         : void
///
void xaxpy(int n, double a, int ix0, double y[9], int iy0) {
    if (!(a == 0.0)) {
        int i;
        i = n - 1;
        for (int k{0}; k <= i; k++) {
            int i1;
            i1 = (iy0 + k) - 1;
            y[i1] += a * y[(ix0 + k) - 1];
        }
    }
}

}  // namespace blas
}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for xaxpy.cpp
///
/// [EOF]
///
