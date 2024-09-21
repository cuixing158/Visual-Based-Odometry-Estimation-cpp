///
/// @file           : slamMapGraph2D_types.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef SLAMMAPGRAPH2D_TYPES_H
#define SLAMMAPGRAPH2D_TYPES_H

/// @include file    : Include Files
#include "rtwtypes.h"

/// Type Definitions
namespace SlamGraph2D {
struct struct0_T {
    double MaxNumEdges;
    double MaxNumNodes;
};

struct slamMapGraph2DPersistentData {
    double freq;
    bool freq_not_empty;
};

struct slamMapGraph2DStackData {
    slamMapGraph2DPersistentData *pd;
};

}  // namespace SlamGraph2D

#endif
///
/// File trailer for slamMapGraph2D_types.h
///
/// [EOF]
///
