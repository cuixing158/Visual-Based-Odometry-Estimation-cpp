///
/// @file           : sqrt.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef SQRT_H
#define SQRT_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace scalar {
void b_sqrt(creal_T &x);

}
} // namespace internal
} // namespace coder
} // namespace SlamGraph2D

#endif
///
/// File trailer for sqrt.h
///
/// [EOF]
///
