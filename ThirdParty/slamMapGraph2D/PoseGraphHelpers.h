///
/// @file           : PoseGraphHelpers.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef POSEGRAPHHELPERS_H
#define POSEGRAPHHELPERS_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <algorithm>
#include <cstddef>
#include <cstdlib>

/// Type Declarations
namespace SlamGraph2D {
namespace coder {
class sparse;

}
} // namespace SlamGraph2D

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace nav {
namespace algs {
namespace internal {
class PoseGraphHelpers {
public:
  static double
  poseGraphCost(const ::coder::array<double, 2U> &posesMat,
                const ::coder::array<double, 2U> &args_edgeNodePairs,
                const ::coder::array<double, 2U> &args_edgeMeasurements,
                const ::coder::array<double, 2U> &args_edgeInfoMats,
                const double args_tformSize[2],
                const double args_infoMatSize[2], double args_poseDeltaLength,
                const ::coder::array<double, 1U> &args_nodeMap,
                const ::coder::array<double, 1U> &args_nodeDims,
                const ::coder::array<bool, 1U> &args_IsLandmarkNode,
                ::coder::array<double, 1U> &gradient, sparse &hessian);
  static double costBetweenTwoNodes(
      const ::coder::array<double, 2U> &Toi,
      const ::coder::array<double, 2U> &Toj,
      const ::coder::array<double, 2U> &measurement,
      const ::coder::array<double, 2U> &Omega, bool nodejIsLandmark,
      double gradi_data[], int &gradi_size, double gradj_data[],
      int &gradj_size, double hessii_data[], int hessii_size[2],
      double hessij_data[], int hessij_size[2], double hessji_data[],
      int hessji_size[2], double hessjj_data[], int hessjj_size[2]);
};

} // namespace internal
} // namespace algs
} // namespace nav
} // namespace coder
} // namespace SlamGraph2D

#endif
///
/// File trailer for PoseGraphHelpers.h
///
/// [EOF]
///
