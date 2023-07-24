///
/// @file           : xdotc.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

/// @include file    : Include Files
#include "xdotc.h"
#include "rt_nonfinite.h"

/// Function Definitions
///
/// @fn             : xdotc
/// @brief          :
/// @param          : int n
///                   const double x[9]
///                   int ix0
///                   const double y[9]
///                   int iy0
/// @return         : double
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
double xdotc(int n, const double x[9], int ix0, const double y[9], int iy0)
{
  double d;
  int i;
  d = 0.0;
  i = static_cast<unsigned char>(n);
  for (int k{0}; k < i; k++) {
    d += x[(ix0 + k) - 1] * y[(iy0 + k) - 1];
  }
  return d;
}

} // namespace blas
} // namespace internal
} // namespace coder
} // namespace SlamGraph2D

///
/// File trailer for xdotc.cpp
///
/// [EOF]
///
