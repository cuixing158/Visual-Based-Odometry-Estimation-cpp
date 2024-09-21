///
/// @file           : svd.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef SVD_H
#define SVD_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
namespace internal {
void svd(const double A[9], double U[9], double s[3], double V[9]);

}
}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for svd.h
///
/// [EOF]
///
