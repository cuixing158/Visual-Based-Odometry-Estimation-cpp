///
/// @file           : sparse1.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef SPARSE1_H
#define SPARSE1_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
class sparse {
public:
  void ctranspose(sparse &y) const;
  void fillIn();
  ::coder::array<double, 1U> d;
  ::coder::array<int, 1U> colidx;
  ::coder::array<int, 1U> rowidx;
  int m;
  int n;
};

} // namespace coder
} // namespace SlamGraph2D

#endif
///
/// File trailer for sparse1.h
///
/// [EOF]
///
