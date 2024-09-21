///
/// @file           : xnrm2.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "xnrm2.h"
#include "rt_nonfinite.h"
#include <cmath>

/// Function Definitions
///
/// @fn             : xnrm2
/// @brief          :
/// @param          : int n
///                   const double x[9]
///                   int ix0
/// @return         : double
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
double xnrm2(int n, const double x[9], int ix0) {
    double scale;
    double y;
    int kend;
    y = 0.0;
    scale = 3.3121686421112381E-170;
    kend = (ix0 + n) - 1;
    for (int k{ix0}; k <= kend; k++) {
        double absxk;
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
            double t;
            t = scale / absxk;
            y = y * t * t + 1.0;
            scale = absxk;
        } else {
            double t;
            t = absxk / scale;
            y += t * t;
        }
    }
    return scale * std::sqrt(y);
}

///
/// @fn             : xnrm2
/// @brief          :
/// @param          : const double x[3]
/// @return         : double
///
double xnrm2(const double x[3]) {
    double scale;
    double y;
    y = 0.0;
    scale = 3.3121686421112381E-170;
    for (int k{2}; k < 4; k++) {
        double absxk;
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
            double t;
            t = scale / absxk;
            y = y * t * t + 1.0;
            scale = absxk;
        } else {
            double t;
            t = absxk / scale;
            y += t * t;
        }
    }
    return scale * std::sqrt(y);
}

}  // namespace blas
}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for xnrm2.cpp
///
/// [EOF]
///
