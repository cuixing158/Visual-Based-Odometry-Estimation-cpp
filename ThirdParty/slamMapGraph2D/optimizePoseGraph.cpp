///
/// @file           : optimizePoseGraph.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

/// @include file    : Include Files
#include "optimizePoseGraph.h"
#include "BlockMatrix.h"
#include "PoseGraphOptimizer.h"
#include "myGraph.h"
#include "poseGraph.h"
#include "rt_nonfinite.h"
#include "sparse1.h"

/// Function Definitions
///
/// @fn             : optimizePoseGraph
/// @brief          :
/// @param          : myGraph *aInstancePtr
///                   poseGraph *b_poseGraph
///                   robotics::core::internal::BlockMatrix &iobj_0
///                   poseGraph &iobj_1
/// @return         : poseGraph *
///
namespace SlamGraph2D {
namespace coder {
poseGraph *optimizePoseGraph(myGraph *aInstancePtr, poseGraph *b_poseGraph,
                             robotics::core::internal::BlockMatrix &iobj_0,
                             poseGraph &iobj_1)
{
  poseGraph *b_poseGraphUpdated;
  sparse hessian;
  double paramStruct_FirstNodePose[3];
  paramStruct_FirstNodePose[0] = 0.0;
  paramStruct_FirstNodePose[1] = 0.0;
  paramStruct_FirstNodePose[2] = 0.0;
  nav::algs::internal::PoseGraphOptimizer::optimize(
      aInstancePtr, b_poseGraph, paramStruct_FirstNodePose, (&iobj_0)[0],
      iobj_1, &b_poseGraphUpdated, hessian);
  return b_poseGraphUpdated;
}

} // namespace coder
} // namespace SlamGraph2D

///
/// File trailer for optimizePoseGraph.cpp
///
/// [EOF]
///
