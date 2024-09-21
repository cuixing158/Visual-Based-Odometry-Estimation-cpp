///
/// @file           : xrot.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef XROT_H
#define XROT_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
void xrot(double x[9], int ix0, int iy0, double c, double s);

}
}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for xrot.h
///
/// [EOF]
///
