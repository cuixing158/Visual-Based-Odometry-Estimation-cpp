///
/// @file           : TrustRegionIndefiniteDogLegInterface.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef TRUSTREGIONINDEFINITEDOGLEGINTERFACE_H
#define TRUSTREGIONINDEFINITEDOGLEGINTERFACE_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Function Declarations
namespace SlamGraph2D {
void binary_expand_op(::coder::array<double, 1U> &in1,
                      const ::coder::array<double, 1U> &in2, double in3);

double binary_expand_op(double in1, double in2,
                        const ::coder::array<double, 1U> &in3,
                        const ::coder::array<double, 1U> &in4,
                        const ::coder::array<double, 1U> &in5);

void minus(::coder::array<double, 1U> &in1,
           const ::coder::array<double, 1U> &in2,
           const ::coder::array<double, 1U> &in3);

} // namespace SlamGraph2D

#endif
///
/// File trailer for TrustRegionIndefiniteDogLegInterface.h
///
/// [EOF]
///
