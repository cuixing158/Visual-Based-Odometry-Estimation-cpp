///
/// @file           : xswap.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef XSWAP_H
#define XSWAP_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
void xswap(double x[9], int ix0, int iy0);

}
}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for xswap.h
///
/// [EOF]
///
