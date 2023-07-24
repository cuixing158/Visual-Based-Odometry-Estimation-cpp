///
/// @file           : TrustRegionIndefiniteDogLegSE2.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#ifndef TRUSTREGIONINDEFINITEDOGLEGSE2_H
#define TRUSTREGIONINDEFINITEDOGLEGSE2_H

/// @include file    : Include Files
#include "SystemTimeProvider.h"
#include "rtwtypes.h"
#include "slamMapGraph2D_internal_types.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Type Declarations
namespace SlamGraph2D {
class myGraph;

namespace coder {
namespace robotics {
namespace core {
namespace internal {
class BlockMatrix;

}
} // namespace core
} // namespace robotics
class sparse;

} // namespace coder
} // namespace SlamGraph2D

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
class TrustRegionIndefiniteDogLegSE2 {
public:
  double solve(myGraph *aInstancePtr, const ::coder::array<double, 2U> &seed,
               BlockMatrix &iobj_0, BlockMatrix **xSol, sparse &hess,
               double &solutionInfo_Error, double &solutionInfo_ExitFlag);
  bool computeBasicSteps(const ::coder::array<double, 1U> &grad,
                         const sparse &B, ::coder::array<double, 1U> &stepSD,
                         ::coder::array<double, 1U> &stepGN) const;

protected:
  void incrementX(const ::coder::array<double, 2U> &x,
                  const ::coder::array<double, 1U> &epsilons,
                  ::coder::array<double, 2U> &xNew) const;

public:
  b_struct_T ExtraArgs;
  double MaxNumIteration;
  double MaxTime;
  ::coder::array<double, 2U> SeedInternal;
  double MaxTimeInternal;
  double MaxNumIterationInternal;
  double StepTolerance;
  SystemTimeProvider TimeObj;
  double GradientTolerance;
  double FunctionTolerance;
  double InitialTrustRegionRadius;
  double TrustRegionRadiusTolerance;
};

} // namespace internal
} // namespace core
} // namespace robotics
} // namespace coder
} // namespace SlamGraph2D

#endif
///
/// File trailer for TrustRegionIndefiniteDogLegSE2.h
///
/// [EOF]
///
