///
/// @file           : slamMapGraph2D_internal_types.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef SLAMMAPGRAPH2D_INTERNAL_TYPES_H
#define SLAMMAPGRAPH2D_INTERNAL_TYPES_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "slamMapGraph2D_types.h"
#include "coder_array.h"

/// Type Definitions
namespace SlamGraph2D {
struct struct_T {
  int xstart;
  int xend;
  int depth;
};

struct b_struct_T {
  ::coder::array<double, 2U> edgeNodePairs;
  ::coder::array<double, 2U> edgeMeasurements;
  ::coder::array<double, 2U> edgeInfoMats;
  double tformSize[2];
  double infoMatSize[2];
  double poseDeltaLength;
  ::coder::array<double, 1U> nodeMap;
  ::coder::array<double, 1U> nodeDims;
  ::coder::array<bool, 1U> IsLandmarkNode;
};

struct c_struct_T {
  ::coder::array<int, 1U> a;
  ::coder::array<int, 1U> b;
};

} // namespace SlamGraph2D

#endif
///
/// File trailer for slamMapGraph2D_internal_types.h
///
/// [EOF]
///
