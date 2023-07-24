///
/// @file           : xrotg.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

/// @include file    : Include Files
#include "xrotg.h"
#include "rt_nonfinite.h"
#include <cmath>

/// Function Definitions
///
/// @fn             : xrotg
/// @brief          :
/// @param          : double &a
///                   double &b
///                   double &s
/// @return         : double
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
double xrotg(double &a, double &b, double &s)
{
  double absa;
  double absb;
  double c;
  double roe;
  double scale;
  roe = b;
  absa = std::abs(a);
  absb = std::abs(b);
  if (absa > absb) {
    roe = a;
  }
  scale = absa + absb;
  if (scale == 0.0) {
    s = 0.0;
    c = 1.0;
    a = 0.0;
    b = 0.0;
  } else {
    double ads;
    double bds;
    ads = absa / scale;
    bds = absb / scale;
    scale *= std::sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      scale = -scale;
    }
    c = a / scale;
    s = b / scale;
    if (absa > absb) {
      b = s;
    } else if (c != 0.0) {
      b = 1.0 / c;
    } else {
      b = 1.0;
    }
    a = scale;
  }
  return c;
}

} // namespace blas
} // namespace internal
} // namespace coder
} // namespace SlamGraph2D

///
/// File trailer for xrotg.cpp
///
/// [EOF]
///
