///
/// @file           : sparse.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef SPARSE_H
#define SPARSE_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Type Declarations
namespace SlamGraph2D {
namespace coder {
class sparse;

}
} // namespace SlamGraph2D

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
void b_sparse(const ::coder::array<double, 1U> &varargin_1,
              const ::coder::array<double, 1U> &varargin_2,
              const ::coder::array<double, 1U> &varargin_3, sparse &y);

}
} // namespace SlamGraph2D

#endif
///
/// File trailer for sparse.h
///
/// [EOF]
///
