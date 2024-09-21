///
/// @file           : slamMapGraph2D_rtwutil.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "slamMapGraph2D_rtwutil.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include <cmath>

/// Function Definitions
///
/// @fn             : rt_atan2d_snf
/// @brief          :
/// @param          : double u0
///                   double u1
/// @return         : double
///
namespace SlamGraph2D {
double rt_atan2d_snf(double u0, double u1) {
    double y;
    if (std::isnan(u0) || std::isnan(u1)) {
        y = rtNaN;
    } else if (std::isinf(u0) && std::isinf(u1)) {
        int i;
        int i1;
        if (u0 > 0.0) {
            i = 1;
        } else {
            i = -1;
        }
        if (u1 > 0.0) {
            i1 = 1;
        } else {
            i1 = -1;
        }
        y = std::atan2(static_cast<double>(i), static_cast<double>(i1));
    } else if (u1 == 0.0) {
        if (u0 > 0.0) {
            y = RT_PI / 2.0;
        } else if (u0 < 0.0) {
            y = -(RT_PI / 2.0);
        } else {
            y = 0.0;
        }
    } else {
        y = std::atan2(u0, u1);
    }
    return y;
}

}  // namespace SlamGraph2D

///
/// File trailer for slamMapGraph2D_rtwutil.cpp
///
/// [EOF]
///
