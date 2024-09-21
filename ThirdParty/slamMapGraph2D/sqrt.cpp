///
/// @file           : sqrt.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "sqrt.h"
#include "rt_nonfinite.h"
#include <cmath>

/// Function Declarations
namespace SlamGraph2D {
static double rt_hypotd_snf(double u0, double u1);

}

/// Function Definitions
///
/// @fn             : rt_hypotd_snf
/// @brief          :
/// @param          : double u0
///                   double u1
/// @return         : double
///
namespace SlamGraph2D {
static double rt_hypotd_snf(double u0, double u1) {
    double a;
    double b;
    double y;
    a = std::abs(u0);
    b = std::abs(u1);
    if (a < b) {
        a /= b;
        y = b * std::sqrt(a * a + 1.0);
    } else if (a > b) {
        b /= a;
        y = a * std::sqrt(b * b + 1.0);
    } else if (std::isnan(b)) {
        y = rtNaN;
    } else {
        y = a * 1.4142135623730951;
    }
    return y;
}

///
/// @fn             : b_sqrt
/// @brief          :
/// @param          : creal_T &x
/// @return         : void
///
namespace coder {
namespace internal {
namespace scalar {
void b_sqrt(creal_T &x) {
    double absxi;
    double absxr;
    double xi;
    double xr;
    xr = x.re;
    xi = x.im;
    if (xi == 0.0) {
        if (xr < 0.0) {
            absxr = 0.0;
            absxi = std::sqrt(-xr);
        } else {
            absxr = std::sqrt(xr);
            absxi = 0.0;
        }
    } else if (xr == 0.0) {
        if (xi < 0.0) {
            absxr = std::sqrt(-xi / 2.0);
            absxi = -absxr;
        } else {
            absxr = std::sqrt(xi / 2.0);
            absxi = absxr;
        }
    } else if (std::isnan(xr)) {
        absxr = rtNaN;
        absxi = rtNaN;
    } else if (std::isnan(xi)) {
        absxr = rtNaN;
        absxi = rtNaN;
    } else if (std::isinf(xi)) {
        absxr = std::abs(xi);
        absxi = xi;
    } else if (std::isinf(xr)) {
        if (xr < 0.0) {
            absxr = 0.0;
            absxi = xi * -xr;
        } else {
            absxr = xr;
            absxi = 0.0;
        }
    } else {
        absxr = std::abs(xr);
        absxi = std::abs(xi);
        if ((absxr > 4.4942328371557893E+307) ||
            (absxi > 4.4942328371557893E+307)) {
            absxr *= 0.5;
            absxi = rt_hypotd_snf(absxr, absxi * 0.5);
            if (absxi > absxr) {
                absxr = std::sqrt(absxi) * std::sqrt(absxr / absxi + 1.0);
            } else {
                absxr = std::sqrt(absxi) * 1.4142135623730951;
            }
        } else {
            absxr = std::sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
        }
        if (xr > 0.0) {
            absxi = 0.5 * (xi / absxr);
        } else {
            if (xi < 0.0) {
                absxi = -absxr;
            } else {
                absxi = absxr;
            }
            absxr = 0.5 * (xi / absxi);
        }
    }
    x.re = absxr;
    x.im = absxi;
}

}  // namespace scalar
}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for sqrt.cpp
///
/// [EOF]
///
