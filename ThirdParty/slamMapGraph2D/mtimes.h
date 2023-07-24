///
/// @file           : mtimes.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef MTIMES_H
#define MTIMES_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
void mtimes(const double A_data[], const int A_size[2],
            const ::coder::array<double, 2U> &B, double C_data[],
            int C_size[2]);

void mtimes(const double A_data[], const int A_size[2], const double B_data[],
            const int B_size[2], double C_data[], int C_size[2]);

} // namespace blas
} // namespace internal
} // namespace coder
} // namespace SlamGraph2D

#endif
///
/// File trailer for mtimes.h
///
/// [EOF]
///
