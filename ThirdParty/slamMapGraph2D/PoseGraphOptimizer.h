///
/// @file           : PoseGraphOptimizer.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef POSEGRAPHOPTIMIZER_H
#define POSEGRAPHOPTIMIZER_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

/// Type Declarations
namespace SlamGraph2D {
class myGraph;

namespace coder {
class poseGraph;

namespace robotics {
namespace core {
namespace internal {
class BlockMatrix;

}
}  // namespace core
}  // namespace robotics
class sparse;

}  // namespace coder
}  // namespace SlamGraph2D

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace nav {
namespace algs {
namespace internal {
class PoseGraphOptimizer {
   public:
    static void optimize(myGraph *aInstancePtr, poseGraph *b_poseGraph,
                         const double paramStruct_FirstNodePose[3],
                         robotics::core::internal::BlockMatrix &iobj_0,
                         poseGraph &iobj_1, poseGraph **poseGraphUpdated,
                         sparse &hessian);
};

}  // namespace internal
}  // namespace algs
}  // namespace nav
}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for PoseGraphOptimizer.h
///
/// [EOF]
///
