///
/// @file           : myGraph.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "myGraph.h"
#include "BlockMatrix.h"
#include "poseGraph.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_types.h"
#include "slamPoseGraph.h"
#include "coder_array.h"

/// Function Definitions
///
/// @fn             : getStackData
/// @brief          :
/// @param          : void
/// @return         : slamMapGraph2DStackData *
///
namespace SlamGraph2D {
slamMapGraph2DStackData *myGraph::getStackData() {
    return &SD_;
}

///
/// @fn             : myGraph
/// @brief          :
/// @param          : void
/// @return         : void
///
myGraph::myGraph() {
    SD_.pd = &pd_;
    pd_.freq_not_empty = false;
}

///
/// @fn             : myGraph
/// @brief          :
/// @param          : void
/// @return         : void
///
myGraph::~myGraph() = default;

///
/// @fn             : slamMapGraph2D
/// @brief          : pose graph
///
/// @param          : const struct0_T *poseParams
///                   const double measurement[3]
///                   double fromID
///                   double toID
///                   const double guesspose[3]
///                   ::coder::array<double, 2U> &poses1
/// @return         : void
///
void myGraph::slamMapGraph2D(const struct0_T *poseParams,
                             const double measurement[3], double fromID,
                             double toID, const double[3],
                             ::coder::array<double, 2U> &poses1) {
    coder::poseGraph lobj_2;
    coder::robotics::core::internal::BlockMatrix lobj_1[3];
    slamPoseGraph pg;
    pg.init(poseParams->MaxNumEdges, poseParams->MaxNumNodes);
    pg.addRelativePose2D(measurement, fromID, toID);
    pg.optimizePoseGraph2D(this, lobj_1[0], lobj_2);
    pg.nodeEstimates2D(poses1);
    //  factor graph
    //  fg = SlamGraph2D.slamFactorGraph();
    //  fg.addRelFactor2D([fromID,toID],measurement,guesspose);
    //  fg.optimizeFactorGraph2D();
    //  poses2 = fg.getNodeState2D();
}

}  // namespace SlamGraph2D

///
/// File trailer for myGraph.cpp
///
/// [EOF]
///
