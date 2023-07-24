///
/// @file           : slamPoseGraph.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

/// @include file    : Include Files
#include "slamPoseGraph.h"
#include "BlockMatrix.h"
#include "myGraph.h"
#include "optimizePoseGraph.h"
#include "poseGraph.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

/// Function Definitions
///
/// @fn             : addRelativePose2D
/// @brief          :
/// @param          : const double measurement[3]
///                   double fromNodeID
///                   double toNodeID
/// @return         : void
///
namespace SlamGraph2D {
void slamPoseGraph::addRelativePose2D(const double measurement[3],
                                      double fromNodeID, double toNodeID)
{
  graph->addRelativePose(measurement, fromNodeID, toNodeID);
}

///
/// @fn             : init
/// @brief          :
/// @param          : double MaxNumEdges
///                   double MaxNumNodes
/// @return         : slamPoseGraph *
///
slamPoseGraph *slamPoseGraph::init(double MaxNumEdges, double MaxNumNodes)
{
  slamPoseGraph *obj;
  ::coder::array<double, 2U> b_obj;
  ::coder::array<double, 2U> c_obj;
  obj = this;
  obj->graph = obj->_pobj0.init(MaxNumEdges, MaxNumNodes);
  obj->graph->get_LoopClosureEdgeIDs(b_obj);
  obj->graph->get_LandmarkNodeIDs(c_obj);
  return obj;
}

///
/// @fn             : nodeEstimates2D
/// @brief          :
/// @param          : ::coder::array<double, 2U> &measurements
/// @return         : void
///
void slamPoseGraph::nodeEstimates2D(
    ::coder::array<double, 2U> &measurements) const
{
  graph->nodeEstimates(measurements);
}

///
/// @fn             : optimizePoseGraph2D
/// @brief          :
/// @param          : myGraph *aInstancePtr
///                   coder::robotics::core::internal::BlockMatrix &iobj_0
///                   coder::poseGraph &iobj_1
/// @return         : void
///
void slamPoseGraph::optimizePoseGraph2D(
    myGraph *aInstancePtr, coder::robotics::core::internal::BlockMatrix &iobj_0,
    coder::poseGraph &iobj_1)
{
  ::coder::array<double, 2U> b_obj;
  ::coder::array<double, 2U> obj;
  graph = coder::optimizePoseGraph(aInstancePtr, graph, (&iobj_0)[0], iobj_1);
  graph->get_LoopClosureEdgeIDs(obj);
  graph->get_LandmarkNodeIDs(b_obj);
}

} // namespace SlamGraph2D

///
/// File trailer for slamPoseGraph.cpp
///
/// [EOF]
///
