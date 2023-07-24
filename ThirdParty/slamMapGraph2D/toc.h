///
/// @file           : toc.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef TOC_H
#define TOC_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

/// Type Declarations
namespace SlamGraph2D {
class myGraph;

}

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
double toc(myGraph *aInstancePtr, double tstart_tv_sec, double tstart_tv_nsec);

}
} // namespace SlamGraph2D

#endif
///
/// File trailer for toc.h
///
/// [EOF]
///
