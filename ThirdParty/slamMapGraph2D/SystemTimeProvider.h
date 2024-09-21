///
/// @file           : SystemTimeProvider.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef SYSTEMTIMEPROVIDER_H
#define SYSTEMTIMEPROVIDER_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_posix_time.h"
#include <cstddef>
#include <cstdlib>

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
class SystemTimeProvider {
   public:
    coderTimespec StartTime;
};

}  // namespace internal
}  // namespace core
}  // namespace robotics
}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for SystemTimeProvider.h
///
/// [EOF]
///
