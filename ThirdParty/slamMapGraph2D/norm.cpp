///
/// @file           : norm.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

/// @include file    : Include Files
#include "norm.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

/// Function Definitions
///
/// @fn             : b_norm
/// @brief          :
/// @param          : const ::coder::array<double, 1U> &x
/// @return         : double
///
namespace SlamGraph2D {
namespace coder {
double b_norm(const ::coder::array<double, 1U> &x)
{
  double y;
  if (x.size(0) == 0) {
    y = 0.0;
  } else {
    y = 0.0;
    if (x.size(0) == 1) {
      y = std::abs(x[0]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = x.size(0);
      for (int k{0}; k < kend; k++) {
        double absxk;
        absxk = std::abs(x[k]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * std::sqrt(y);
    }
  }
  return y;
}

} // namespace coder
} // namespace SlamGraph2D

///
/// File trailer for norm.cpp
///
/// [EOF]
///
