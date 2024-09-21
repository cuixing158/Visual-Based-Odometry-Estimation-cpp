///
/// @file           : BlockInserter2.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef BLOCKINSERTER2_H
#define BLOCKINSERTER2_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace nav {
namespace algs {
namespace internal {
class BlockInserter2 {
   public:
    void insertGradientBlock(double i, const double blocki_data[],
                             int blocki_size);
    ::coder::array<double, 1U> Gradient;
    ::coder::array<double, 1U> NodeDims;
    ::coder::array<double, 1U> NodeMap;
    ::coder::array<double, 2U> HessianCSC;
    double HessianCSCCount;
};

}  // namespace internal
}  // namespace algs
}  // namespace nav
}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for BlockInserter2.h
///
/// [EOF]
///
