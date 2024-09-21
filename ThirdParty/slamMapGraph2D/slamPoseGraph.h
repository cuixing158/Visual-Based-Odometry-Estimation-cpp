///
/// @file           : slamPoseGraph.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef SLAMPOSEGRAPH_H
#define SLAMPOSEGRAPH_H

/// @include file    : Include Files
#include "poseGraph.h"
#include "rtwtypes.h"
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
}  // namespace core
}  // namespace robotics
}  // namespace coder
}  // namespace SlamGraph2D

/// Type Definitions
namespace SlamGraph2D {
class slamPoseGraph {
   public:
    slamPoseGraph *init(double MaxNumEdges, double MaxNumNodes);
    void addRelativePose2D(const double measurement[3], double fromNodeID,
                           double toNodeID);
    void optimizePoseGraph2D(myGraph *aInstancePtr,
                             coder::robotics::core::internal::BlockMatrix &iobj_0,
                             coder::poseGraph &iobj_1);
    void nodeEstimates2D(::coder::array<double, 2U> &measurements) const;
    coder::poseGraph _pobj0;

   private:
    coder::poseGraph *graph;
};

}  // namespace SlamGraph2D

#endif
///
/// File trailer for slamPoseGraph.h
///
/// [EOF]
///
