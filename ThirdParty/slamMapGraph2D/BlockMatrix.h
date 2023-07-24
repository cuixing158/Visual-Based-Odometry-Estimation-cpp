///
/// @file           : BlockMatrix.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
class BlockMatrix {
public:
  void replaceBlock(double i, const double blockij[9]);
  void extractBlock(double i, ::coder::array<double, 2U> &B) const;
  ::coder::array<double, 2U> Matrix;
  double NumRowBlocks;
  double NumColBlocks;
  double BlockSize[2];
};

class b_BlockMatrix {
public:
  ::coder::array<double, 2U> Matrix;
  double BlockSize[2];
};

} // namespace internal
} // namespace core
} // namespace robotics
} // namespace coder
} // namespace SlamGraph2D

#endif
///
/// File trailer for BlockMatrix.h
///
/// [EOF]
///
