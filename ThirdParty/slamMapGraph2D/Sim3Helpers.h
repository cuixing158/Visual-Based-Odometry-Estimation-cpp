///
/// @file           : Sim3Helpers.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef SIM3HELPERS_H
#define SIM3HELPERS_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
class Sim3Helpers {
public:
  static void multiplyLogSim3(const double S1[16], const double S2[16],
                              const double S3[16], double e[7]);
  static void sim3ToSform(const double minVecSim3[7], double S[16]);
};

} // namespace internal
} // namespace core
} // namespace robotics
} // namespace coder
} // namespace SlamGraph2D

#endif
///
/// File trailer for Sim3Helpers.h
///
/// [EOF]
///
