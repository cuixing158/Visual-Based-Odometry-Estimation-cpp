///
/// @file           : norm.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef NORM_H
#define NORM_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
double b_norm(const ::coder::array<double, 1U> &x);

}
}  // namespace SlamGraph2D

#endif
///
/// File trailer for norm.h
///
/// [EOF]
///
