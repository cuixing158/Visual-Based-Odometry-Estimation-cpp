///
/// @file           : toc.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "toc.h"
#include "myGraph.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_types.h"
#include "coder_posix_time.h"

/// Function Definitions
///
/// @fn             : toc
/// @brief          :
/// @param          : myGraph *aInstancePtr
///                   double tstart_tv_sec
///                   double tstart_tv_nsec
/// @return         : double
///
namespace SlamGraph2D {
namespace coder {
double toc(myGraph *aInstancePtr, double tstart_tv_sec, double tstart_tv_nsec) {
    coderTimespec b_timespec;
    slamMapGraph2DStackData *localSD;
    localSD = aInstancePtr->getStackData();
    if (!localSD->pd->freq_not_empty) {
        localSD->pd->freq_not_empty = true;
        coderInitTimeFunctions(&localSD->pd->freq);
    }
    coderTimeClockGettimeMonotonic(&b_timespec, localSD->pd->freq);
    return (b_timespec.tv_sec - tstart_tv_sec) +
           (b_timespec.tv_nsec - tstart_tv_nsec) / 1.0E+9;
}

}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for toc.cpp
///
/// [EOF]
///
