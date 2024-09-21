///
/// @file           : CXSparseAPI.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef CXSPARSEAPI_H
#define CXSPARSEAPI_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include "solve_from_qr.h"
#include <cstddef>
#include <cstdlib>

/// Type Declarations
namespace SlamGraph2D {
namespace coder {
class sparse;

}
}  // namespace SlamGraph2D

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace internal {
class CXSparseAPI {
   public:
    static void iteratedQR(const sparse &A, const ::coder::array<double, 1U> &b,
                           int n, ::coder::array<double, 1U> &out);
};

}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for CXSparseAPI.h
///
/// [EOF]
///
