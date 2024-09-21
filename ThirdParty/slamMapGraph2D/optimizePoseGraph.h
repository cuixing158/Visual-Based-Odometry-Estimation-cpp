///
/// @file           : optimizePoseGraph.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef OPTIMIZEPOSEGRAPH_H
#define OPTIMIZEPOSEGRAPH_H

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
}  // namespace coder
}  // namespace SlamGraph2D

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
poseGraph *optimizePoseGraph(myGraph *aInstancePtr, poseGraph *b_poseGraph,
                             robotics::core::internal::BlockMatrix &iobj_0,
                             poseGraph &iobj_1);

}
}  // namespace SlamGraph2D

#endif
///
/// File trailer for optimizePoseGraph.h
///
/// [EOF]
///
