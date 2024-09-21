///
/// @file           : anonymous_function.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef ANONYMOUS_FUNCTION_H
#define ANONYMOUS_FUNCTION_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "slamMapGraph2D_internal_types.h"
#include <cstddef>
#include <cstdlib>

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
class anonymous_function {
   public:
    c_struct_T workspace;
};

}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for anonymous_function.h
///
/// [EOF]
///
