///
/// @file           : BlockMatrix.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "BlockMatrix.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

/// Function Definitions
///
/// @fn             : extractBlock
/// @brief          :
/// @param          : double i
///                   ::coder::array<double, 2U> &B
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
void BlockMatrix::extractBlock(double i, ::coder::array<double, 2U> &B) const {
    double colStart;
    double d;
    double rowStart;
    int b_i;
    int b_loop_ub;
    int i1;
    int loop_ub;
    rowStart = BlockSize[0] * (i - 1.0) + 1.0;
    colStart = BlockSize[1] * 0.0 + 1.0;
    d = (rowStart + BlockSize[0]) - 1.0;
    if (rowStart > d) {
        b_i = 0;
        i1 = 0;
    } else {
        b_i = static_cast<int>(rowStart) - 1;
        i1 = static_cast<int>(d);
    }
    d = (colStart + BlockSize[1]) - 1.0;
    if (colStart > d) {
        loop_ub = 0;
    } else {
        loop_ub = static_cast<int>(d);
    }
    b_loop_ub = i1 - b_i;
    B.set_size(b_loop_ub, loop_ub);
    for (i1 = 0; i1 < loop_ub; i1++) {
        for (int i2{0}; i2 < b_loop_ub; i2++) {
            B[i2 + B.size(0) * i1] = Matrix[(b_i + i2) + Matrix.size(0) * i1];
        }
    }
}

///
/// @fn             : replaceBlock
/// @brief          :
/// @param          : double i
///                   const double blockij[9]
/// @return         : void
///
void BlockMatrix::replaceBlock(double i, const double blockij[9]) {
    double colStart;
    double d;
    double rowStart;
    int b_i;
    int i1;
    int loop_ub;
    int unnamed_idx_0;
    rowStart = BlockSize[0] * (i - 1.0) + 1.0;
    colStart = BlockSize[1] * 0.0 + 1.0;
    d = (rowStart + BlockSize[0]) - 1.0;
    if (rowStart > d) {
        b_i = 0;
        i1 = 0;
    } else {
        b_i = static_cast<int>(rowStart) - 1;
        i1 = static_cast<int>(d);
    }
    d = (colStart + BlockSize[1]) - 1.0;
    unnamed_idx_0 = i1 - b_i;
    if (colStart > d) {
        loop_ub = 0;
    } else {
        loop_ub = static_cast<int>(d);
    }
    for (i1 = 0; i1 < loop_ub; i1++) {
        for (int i2{0}; i2 < unnamed_idx_0; i2++) {
            Matrix[(b_i + i2) + Matrix.size(0) * i1] =
                blockij[i2 + unnamed_idx_0 * i1];
        }
    }
}

}  // namespace internal
}  // namespace core
}  // namespace robotics
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for BlockMatrix.cpp
///
/// [EOF]
///
