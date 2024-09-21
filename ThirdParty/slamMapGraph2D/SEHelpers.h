///
/// @file           : SEHelpers.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef SEHELPERS_H
#define SEHELPERS_H

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
class SEHelpers {
   public:
    static void veelogmSE3(const double T[16], double vec[6]);
    static void expSE3hat(const double e[6], double T[16]);
};

}  // namespace internal
}  // namespace core
}  // namespace robotics
}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for SEHelpers.h
///
/// [EOF]
///
