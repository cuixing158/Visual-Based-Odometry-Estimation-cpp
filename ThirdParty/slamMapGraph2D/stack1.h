///
/// @file           : stack1.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef STACK1_H
#define STACK1_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "slamMapGraph2D_internal_types.h"
#include "coder_bounded_array.h"
#include <cstddef>
#include <cstdlib>

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace internal {
class stack {
   public:
    ::coder::bounded_array<struct_T, 120U, 1U> d;
    int n;
};

}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for stack1.h
///
/// [EOF]
///
