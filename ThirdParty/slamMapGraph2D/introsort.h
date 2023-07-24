///
/// @file           : introsort.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef INTROSORT_H
#define INTROSORT_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Type Declarations
namespace SlamGraph2D {
namespace coder {
class anonymous_function;

}
} // namespace SlamGraph2D

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
namespace internal {
void introsort(::coder::array<int, 1U> &x, int xend,
               const anonymous_function &cmp);

}
} // namespace coder
} // namespace SlamGraph2D

#endif
///
/// File trailer for introsort.h
///
/// [EOF]
///
