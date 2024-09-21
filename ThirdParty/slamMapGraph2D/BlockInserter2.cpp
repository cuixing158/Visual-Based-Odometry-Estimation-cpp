///
/// @file           : BlockInserter2.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "BlockInserter2.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

/// Function Declarations
namespace SlamGraph2D {
static void binary_expand_op(coder::nav::algs::internal::BlockInserter2 *in1,
                             int in2, int in4, int in5, const double in6_data[],
                             const int &in6_size);

}

/// Function Definitions
///
/// @fn             : binary_expand_op
/// @brief          :
/// @param          : coder::nav::algs::internal::BlockInserter2 *in1
///                   int in2
///                   int in4
///                   int in5
///                   const double in6_data[]
///                   const int &in6_size
/// @return         : void
///
namespace SlamGraph2D {
static void binary_expand_op(coder::nav::algs::internal::BlockInserter2 *in1,
                             int in2, int in4, int in5, const double in6_data[],
                             const int &in6_size) {
    ::coder::array<double, 1U> b_in1;
    int stride_0_0;
    b_in1.set_size(in6_size);
    stride_0_0 = ((in5 - in4) + 1 != 1);
    for (int i{0}; i < in6_size; i++) {
        b_in1[i] = in1->Gradient[in4 + i * stride_0_0] + in6_data[i];
    }
    stride_0_0 = b_in1.size(0);
    for (int i{0}; i < stride_0_0; i++) {
        in1->Gradient[in2 + i] = b_in1[i];
    }
}

///
/// @fn             : insertGradientBlock
/// @brief          :
/// @param          : double i
///                   const double blocki_data[]
///                   int blocki_size
/// @return         : void
///
namespace coder {
namespace nav {
namespace algs {
namespace internal {
void BlockInserter2::insertGradientBlock(double i, const double blocki_data[],
                                         int blocki_size) {
    ::coder::array<double, 1U> obj;
    double d;
    double rowStart;
    int b_i;
    int i1;
    int i2;
    int loop_ub;
    rowStart = NodeMap[static_cast<int>(i) - 1];
    d = (rowStart + NodeDims[static_cast<int>(i) - 1]) - 1.0;
    if (rowStart > d) {
        b_i = 0;
        i1 = 0;
        i2 = 0;
    } else {
        b_i = static_cast<int>(rowStart) - 1;
        i1 = static_cast<int>(d);
        i2 = static_cast<int>(rowStart) - 1;
    }
    loop_ub = i1 - b_i;
    if (loop_ub == blocki_size) {
        obj.set_size(loop_ub);
        for (i1 = 0; i1 < loop_ub; i1++) {
            obj[i1] = Gradient[b_i + i1] + blocki_data[i1];
        }
        loop_ub = obj.size(0);
        for (b_i = 0; b_i < loop_ub; b_i++) {
            Gradient[i2 + b_i] = obj[b_i];
        }
    } else {
        binary_expand_op(this, i2, b_i, i1 - 1, blocki_data, blocki_size);
    }
}

}  // namespace internal
}  // namespace algs
}  // namespace nav
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for BlockInserter2.cpp
///
/// [EOF]
///
