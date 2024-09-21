///
/// @file           : TrustRegionIndefiniteDogLegInterface.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "TrustRegionIndefiniteDogLegInterface.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

/// Function Definitions
///
/// @fn             : binary_expand_op
/// @brief          :
/// @param          : ::coder::array<double, 1U> &in1
///                   const ::coder::array<double, 1U> &in2
///                   double in3
/// @return         : void
///
namespace SlamGraph2D {
void binary_expand_op(::coder::array<double, 1U> &in1,
                      const ::coder::array<double, 1U> &in2, double in3) {
    ::coder::array<double, 1U> b_in2;
    int loop_ub;
    int stride_0_0;
    int stride_1_0;
    if (in1.size(0) == 1) {
        loop_ub = in2.size(0);
    } else {
        loop_ub = in1.size(0);
    }
    b_in2.set_size(loop_ub);
    stride_0_0 = (in2.size(0) != 1);
    stride_1_0 = (in1.size(0) != 1);
    for (int i{0}; i < loop_ub; i++) {
        b_in2[i] = in2[i * stride_0_0] + in3 * in1[i * stride_1_0];
    }
    in1.set_size(b_in2.size(0));
    loop_ub = b_in2.size(0);
    for (int i{0}; i < loop_ub; i++) {
        in1[i] = b_in2[i];
    }
}

///
/// @fn             : binary_expand_op
/// @brief          :
/// @param          : double in1
///                   double in2
///                   const ::coder::array<double, 1U> &in3
///                   const ::coder::array<double, 1U> &in4
///                   const ::coder::array<double, 1U> &in5
/// @return         : double
///
double binary_expand_op(double in1, double in2,
                        const ::coder::array<double, 1U> &in3,
                        const ::coder::array<double, 1U> &in4,
                        const ::coder::array<double, 1U> &in5) {
    ::coder::array<double, 1U> b_in4;
    double b_in3;
    int loop_ub;
    int stride_0_0;
    int stride_1_0;
    if (in5.size(0) == 1) {
        loop_ub = in4.size(0);
    } else {
        loop_ub = in5.size(0);
    }
    b_in4.set_size(loop_ub);
    stride_0_0 = (in4.size(0) != 1);
    stride_1_0 = (in5.size(0) != 1);
    for (int i{0}; i < loop_ub; i++) {
        b_in4[i] = in4[i * stride_0_0] + in5[i * stride_1_0];
    }
    b_in3 = 0.0;
    loop_ub = in3.size(0);
    for (int i{0}; i < loop_ub; i++) {
        b_in3 += -in3[i] * b_in4[i];
    }
    return (in1 - in2) / b_in3;
}

///
/// @fn             : minus
/// @brief          :
/// @param          : ::coder::array<double, 1U> &in1
///                   const ::coder::array<double, 1U> &in2
///                   const ::coder::array<double, 1U> &in3
/// @return         : void
///
void minus(::coder::array<double, 1U> &in1,
           const ::coder::array<double, 1U> &in2,
           const ::coder::array<double, 1U> &in3) {
    int loop_ub;
    int stride_0_0;
    int stride_1_0;
    if (in3.size(0) == 1) {
        loop_ub = in2.size(0);
    } else {
        loop_ub = in3.size(0);
    }
    in1.set_size(loop_ub);
    stride_0_0 = (in2.size(0) != 1);
    stride_1_0 = (in3.size(0) != 1);
    for (int i{0}; i < loop_ub; i++) {
        in1[i] = in2[i * stride_0_0] - in3[i * stride_1_0];
    }
}

}  // namespace SlamGraph2D

///
/// File trailer for TrustRegionIndefiniteDogLegInterface.cpp
///
/// [EOF]
///
