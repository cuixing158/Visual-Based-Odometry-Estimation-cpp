///
/// @file           : poseGraph.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef POSEGRAPH_H
#define POSEGRAPH_H

/// @include file    : Include Files
#include "BlockMatrix.h"
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
class poseGraph {
   public:
    poseGraph *init(double varargin_2, double varargin_4);
    void get_LoopClosureEdgeIDs(::coder::array<double, 2U> &lcIds) const;
    void get_LandmarkNodeIDs(::coder::array<double, 2U> &lmIds) const;
    void addRelativePose(const double varargin_1[3], double varargin_3,
                         double varargin_4);
    void findEdgeID(const double nodePair[2],
                    ::coder::array<double, 1U> &edgeID) const;
    void nodeEstimates(::coder::array<double, 2U> &nodeEsts) const;
    double NumNodes;
    double NumEdges;
    double NumLoopClosureEdges;
    robotics::core::internal::BlockMatrix *NodeEstimates;
    ::coder::array<double, 1U> NodeMap;
    ::coder::array<double, 1U> NodeDims;
    ::coder::array<bool, 1U> IsLandmarkNode;
    ::coder::array<double, 2U> EdgeNodePairs;
    ::coder::array<double, 2U> LoopClosureEdgeNodePairs;
    robotics::core::internal::BlockMatrix *EdgeMeasurements;
    robotics::core::internal::BlockMatrix *EdgeInfoMatrices;
    ::coder::array<double, 2U> LoopClosureEdgeIDsInternal;
    double MaxNumEdges;
    bool __MaxNumEdges_AssignmentSentinel;
    double MaxNumNodes;
    bool __MaxNumNodes_AssignmentSentinel;
    robotics::core::internal::BlockMatrix _pobj0[3];
};

}  // namespace coder
}  // namespace SlamGraph2D

#endif
///
/// File trailer for poseGraph.h
///
/// [EOF]
///
