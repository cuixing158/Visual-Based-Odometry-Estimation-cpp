///
/// @file           : PoseGraphOptimizer.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "PoseGraphOptimizer.h"
#include "BlockMatrix.h"
#include "SystemTimeProvider.h"
#include "TrustRegionIndefiniteDogLegSE2.h"
#include "myGraph.h"
#include "poseGraph.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_internal_types.h"
#include "sparse1.h"
#include "coder_array.h"
#include "coder_posix_time.h"
#include <cmath>

/// Function Definitions
///
/// @fn             : optimize
/// @brief          :
/// @param          : myGraph *aInstancePtr
///                   poseGraph *b_poseGraph
///                   const double paramStruct_FirstNodePose[3]
///                   robotics::core::internal::BlockMatrix &iobj_0
///                   poseGraph &iobj_1
///                   poseGraph **poseGraphUpdated
///                   sparse &hessian
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace nav {
namespace algs {
namespace internal {
void PoseGraphOptimizer::optimize(myGraph *aInstancePtr, poseGraph *b_poseGraph,
                                  const double paramStruct_FirstNodePose[3],
                                  robotics::core::internal::BlockMatrix &iobj_0,
                                  poseGraph &iobj_1,
                                  poseGraph **poseGraphUpdated, sparse &hessian) {
    robotics::core::internal::BlockMatrix lobj_1[2];
    robotics::core::internal::BlockMatrix *obj;
    robotics::core::internal::BlockMatrix *posesUpdated;
    robotics::core::internal::TrustRegionIndefiniteDogLegSE2 solver;
    ::coder::array<double, 2U> T1;
    ::coder::array<double, 2U> varargin_1;
    ::coder::array<double, 1U> c_poseGraph;
    ::coder::array<bool, 1U> d_poseGraph;
    double R[9];
    double T1Offset[9];
    double c_T_tmp[9];
    double R_idx_0;
    double T_tmp;
    double b_T_tmp;
    double colStart;
    double rowStart;
    int boffset;
    int coffset;
    int i;
    int loop_ub;
    rowStart = b_poseGraph->MaxNumEdges;
    colStart = b_poseGraph->MaxNumNodes;
    *poseGraphUpdated = iobj_1.init(rowStart, colStart);
    loop_ub = b_poseGraph->EdgeNodePairs.size(0) << 1;
    (*poseGraphUpdated)
        ->EdgeNodePairs.set_size(b_poseGraph->EdgeNodePairs.size(0), 2);
    c_poseGraph.set_size(loop_ub);
    for (i = 0; i < loop_ub; i++) {
        c_poseGraph[i] = b_poseGraph->EdgeNodePairs[i];
    }
    loop_ub = c_poseGraph.size(0);
    for (i = 0; i < loop_ub; i++) {
        (*poseGraphUpdated)->EdgeNodePairs[i] = c_poseGraph[i];
    }
    loop_ub = b_poseGraph->LoopClosureEdgeNodePairs.size(0) << 1;
    (*poseGraphUpdated)
        ->LoopClosureEdgeNodePairs.set_size(
            b_poseGraph->LoopClosureEdgeNodePairs.size(0), 2);
    c_poseGraph.set_size(loop_ub);
    for (i = 0; i < loop_ub; i++) {
        c_poseGraph[i] = b_poseGraph->LoopClosureEdgeNodePairs[i];
    }
    loop_ub = c_poseGraph.size(0);
    for (i = 0; i < loop_ub; i++) {
        (*poseGraphUpdated)->LoopClosureEdgeNodePairs[i] = c_poseGraph[i];
    }
    loop_ub = b_poseGraph->LoopClosureEdgeIDsInternal.size(1);
    (*poseGraphUpdated)
        ->LoopClosureEdgeIDsInternal.set_size(
            1, b_poseGraph->LoopClosureEdgeIDsInternal.size(1));
    c_poseGraph.set_size(loop_ub);
    for (i = 0; i < loop_ub; i++) {
        c_poseGraph[i] = b_poseGraph->LoopClosureEdgeIDsInternal[i];
    }
    loop_ub = c_poseGraph.size(0);
    for (i = 0; i < loop_ub; i++) {
        (*poseGraphUpdated)->LoopClosureEdgeIDsInternal[i] = c_poseGraph[i];
    }
    obj = b_poseGraph->NodeEstimates;
    varargin_1.set_size(obj->Matrix.size(0), 3);
    loop_ub = obj->Matrix.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        varargin_1[i] = obj->Matrix[i];
    }
    rowStart = obj->BlockSize[0];
    colStart = obj->BlockSize[1];
    (&(&iobj_0)[0])[0].Matrix.set_size(varargin_1.size(0), 3);
    loop_ub = varargin_1.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        (&(&iobj_0)[0])[0].Matrix[i] = varargin_1[i];
    }
    (&(&iobj_0)[0])[0].BlockSize[0] = rowStart;
    (&(&iobj_0)[0])[0].BlockSize[1] = colStart;
    (&(&iobj_0)[0])[0].NumRowBlocks =
        static_cast<double>(varargin_1.size(0)) / rowStart;
    (&(&iobj_0)[0])[0].NumColBlocks = 3.0 / colStart;
    (*poseGraphUpdated)->NodeEstimates = &(&(&iobj_0)[0])[0];
    obj = b_poseGraph->EdgeMeasurements;
    varargin_1.set_size(obj->Matrix.size(0), 3);
    loop_ub = obj->Matrix.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        varargin_1[i] = obj->Matrix[i];
    }
    rowStart = obj->BlockSize[0];
    colStart = obj->BlockSize[1];
    (&(&iobj_0)[0])[1].Matrix.set_size(varargin_1.size(0), 3);
    loop_ub = varargin_1.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        (&(&iobj_0)[0])[1].Matrix[i] = varargin_1[i];
    }
    (&(&iobj_0)[0])[1].BlockSize[0] = rowStart;
    (&(&iobj_0)[0])[1].BlockSize[1] = colStart;
    (&(&iobj_0)[0])[1].NumRowBlocks =
        static_cast<double>(varargin_1.size(0)) / rowStart;
    (&(&iobj_0)[0])[1].NumColBlocks = 3.0 / colStart;
    (*poseGraphUpdated)->EdgeMeasurements = &(&(&iobj_0)[0])[1];
    obj = b_poseGraph->EdgeInfoMatrices;
    varargin_1.set_size(obj->Matrix.size(0), 3);
    loop_ub = obj->Matrix.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        varargin_1[i] = obj->Matrix[i];
    }
    rowStart = obj->BlockSize[0];
    colStart = obj->BlockSize[1];
    (&(&iobj_0)[0])[2].Matrix.set_size(varargin_1.size(0), 3);
    loop_ub = varargin_1.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        (&(&iobj_0)[0])[2].Matrix[i] = varargin_1[i];
    }
    (&(&iobj_0)[0])[2].BlockSize[0] = rowStart;
    (&(&iobj_0)[0])[2].BlockSize[1] = colStart;
    (&(&iobj_0)[0])[2].NumRowBlocks =
        static_cast<double>(varargin_1.size(0)) / rowStart;
    (&(&iobj_0)[0])[2].NumColBlocks = 3.0 / colStart;
    (*poseGraphUpdated)->EdgeInfoMatrices = &(&(&iobj_0)[0])[2];
    (*poseGraphUpdated)->NumNodes = b_poseGraph->NumNodes;
    c_poseGraph.set_size(b_poseGraph->NodeMap.size(0));
    loop_ub = b_poseGraph->NodeMap.size(0);
    for (i = 0; i < loop_ub; i++) {
        c_poseGraph[i] = b_poseGraph->NodeMap[i];
    }
    (*poseGraphUpdated)->NodeMap.set_size(c_poseGraph.size(0));
    loop_ub = c_poseGraph.size(0);
    for (i = 0; i < loop_ub; i++) {
        (*poseGraphUpdated)->NodeMap[i] = c_poseGraph[i];
    }
    c_poseGraph.set_size(b_poseGraph->NodeDims.size(0));
    loop_ub = b_poseGraph->NodeDims.size(0);
    for (i = 0; i < loop_ub; i++) {
        c_poseGraph[i] = b_poseGraph->NodeDims[i];
    }
    (*poseGraphUpdated)->NodeDims.set_size(c_poseGraph.size(0));
    loop_ub = c_poseGraph.size(0);
    for (i = 0; i < loop_ub; i++) {
        (*poseGraphUpdated)->NodeDims[i] = c_poseGraph[i];
    }
    d_poseGraph.set_size(b_poseGraph->IsLandmarkNode.size(0));
    loop_ub = b_poseGraph->IsLandmarkNode.size(0);
    for (i = 0; i < loop_ub; i++) {
        d_poseGraph[i] = b_poseGraph->IsLandmarkNode[i];
    }
    (*poseGraphUpdated)->IsLandmarkNode.set_size(d_poseGraph.size(0));
    loop_ub = d_poseGraph.size(0);
    for (i = 0; i < loop_ub; i++) {
        (*poseGraphUpdated)->IsLandmarkNode[i] = d_poseGraph[i];
    }
    (*poseGraphUpdated)->NumEdges = b_poseGraph->NumEdges;
    (*poseGraphUpdated)->NumLoopClosureEdges = b_poseGraph->NumLoopClosureEdges;
    T_tmp = std::sin(paramStruct_FirstNodePose[2]);
    b_T_tmp = std::cos(paramStruct_FirstNodePose[2]);
    solver.TimeObj.StartTime.tv_sec = 0.0;
    solver.TimeObj.StartTime.tv_nsec = 0.0;
    solver.MaxNumIteration = 300.0;
    solver.MaxTime = 500.0;
    solver.GradientTolerance = 5.0E-9;
    solver.StepTolerance = 1.0E-12;
    solver.FunctionTolerance = 1.0E-8;
    solver.InitialTrustRegionRadius = 100.0;
    solver.TrustRegionRadiusTolerance = 1.0E-10;
    R_idx_0 = (*poseGraphUpdated)->NumEdges;
    if (R_idx_0 < 1.0) {
        loop_ub = 0;
    } else {
        loop_ub = static_cast<int>(R_idx_0);
    }
    solver.ExtraArgs.edgeNodePairs.set_size(loop_ub, 2);
    for (i = 0; i < 2; i++) {
        for (boffset = 0; boffset < loop_ub; boffset++) {
            solver.ExtraArgs
                .edgeNodePairs[boffset + solver.ExtraArgs.edgeNodePairs.size(0) * i] =
                (*poseGraphUpdated)
                    ->EdgeNodePairs[boffset +
                                    (*poseGraphUpdated)->EdgeNodePairs.size(0) * i];
        }
    }
    R_idx_0 = (*poseGraphUpdated)->NumNodes * 3.0;
    if (R_idx_0 < 1.0) {
        loop_ub = 0;
    } else {
        loop_ub = static_cast<int>(R_idx_0);
    }
    varargin_1.set_size(loop_ub, 3);
    for (i = 0; i < 3; i++) {
        for (boffset = 0; boffset < loop_ub; boffset++) {
            varargin_1[boffset + varargin_1.size(0) * i] =
                (*poseGraphUpdated)
                    ->NodeEstimates
                    ->Matrix[boffset +
                             (*poseGraphUpdated)->NodeEstimates->Matrix.size(0) * i];
        }
    }
    R_idx_0 = (*poseGraphUpdated)->NumEdges * 3.0;
    if (R_idx_0 < 1.0) {
        loop_ub = 0;
    } else {
        loop_ub = static_cast<int>(R_idx_0);
    }
    solver.ExtraArgs.edgeInfoMats.set_size(loop_ub, 3);
    for (i = 0; i < 3; i++) {
        for (boffset = 0; boffset < loop_ub; boffset++) {
            solver.ExtraArgs
                .edgeInfoMats[boffset + solver.ExtraArgs.edgeInfoMats.size(0) * i] =
                (*poseGraphUpdated)
                    ->EdgeInfoMatrices
                    ->Matrix[boffset +
                             (*poseGraphUpdated)->EdgeInfoMatrices->Matrix.size(0) *
                                 i];
        }
    }
    R_idx_0 = (*poseGraphUpdated)->NumEdges * 3.0;
    if (R_idx_0 < 1.0) {
        loop_ub = 0;
    } else {
        loop_ub = static_cast<int>(R_idx_0);
    }
    solver.ExtraArgs.edgeMeasurements.set_size(loop_ub, 3);
    for (i = 0; i < 3; i++) {
        for (boffset = 0; boffset < loop_ub; boffset++) {
            solver.ExtraArgs
                .edgeMeasurements[boffset +
                                  solver.ExtraArgs.edgeMeasurements.size(0) * i] =
                (*poseGraphUpdated)
                    ->EdgeMeasurements
                    ->Matrix[boffset +
                             (*poseGraphUpdated)->EdgeMeasurements->Matrix.size(0) *
                                 i];
        }
    }
    solver.ExtraArgs.tformSize[0] = 3.0;
    solver.ExtraArgs.infoMatSize[0] = 3.0;
    solver.ExtraArgs.tformSize[1] = 3.0;
    solver.ExtraArgs.infoMatSize[1] = 3.0;
    solver.ExtraArgs.poseDeltaLength = 3.0;
    solver.ExtraArgs.nodeMap.set_size((*poseGraphUpdated)->NodeMap.size(0));
    loop_ub = (*poseGraphUpdated)->NodeMap.size(0);
    for (i = 0; i < loop_ub; i++) {
        solver.ExtraArgs.nodeMap[i] = (*poseGraphUpdated)->NodeMap[i];
    }
    solver.ExtraArgs.nodeDims.set_size((*poseGraphUpdated)->NodeDims.size(0));
    loop_ub = (*poseGraphUpdated)->NodeDims.size(0);
    for (i = 0; i < loop_ub; i++) {
        solver.ExtraArgs.nodeDims[i] = (*poseGraphUpdated)->NodeDims[i];
    }
    solver.ExtraArgs.IsLandmarkNode.set_size(
        (*poseGraphUpdated)->IsLandmarkNode.size(0));
    loop_ub = (*poseGraphUpdated)->IsLandmarkNode.size(0);
    for (i = 0; i < loop_ub; i++) {
        solver.ExtraArgs.IsLandmarkNode[i] = (*poseGraphUpdated)->IsLandmarkNode[i];
    }
    solver.solve(aInstancePtr, varargin_1, lobj_1[0], &posesUpdated, hessian,
                 colStart, R_idx_0);
    rowStart = posesUpdated->BlockSize[0] * 0.0 + 1.0;
    colStart = posesUpdated->BlockSize[1] * 0.0 + 1.0;
    R_idx_0 = (rowStart + posesUpdated->BlockSize[0]) - 1.0;
    if (rowStart > R_idx_0) {
        loop_ub = 0;
    } else {
        loop_ub = static_cast<int>(R_idx_0);
    }
    R_idx_0 = (colStart + posesUpdated->BlockSize[1]) - 1.0;
    if (colStart > R_idx_0) {
        coffset = 0;
    } else {
        coffset = static_cast<int>(R_idx_0);
    }
    T1.set_size(loop_ub, coffset);
    for (i = 0; i < coffset; i++) {
        for (boffset = 0; boffset < loop_ub; boffset++) {
            T1[boffset + T1.size(0) * i] =
                posesUpdated->Matrix[boffset + posesUpdated->Matrix.size(0) * i];
        }
    }
    double R_idx_3;
    R_idx_0 = T1[0];
    rowStart = T1[T1.size(0)];
    colStart = T1[1];
    R_idx_3 = T1[T1.size(0) + 1];
    c_T_tmp[0] = b_T_tmp;
    c_T_tmp[3] = -T_tmp;
    c_T_tmp[6] = paramStruct_FirstNodePose[0];
    c_T_tmp[1] = T_tmp;
    c_T_tmp[4] = b_T_tmp;
    c_T_tmp[7] = paramStruct_FirstNodePose[1];
    c_T_tmp[2] = 0.0;
    c_T_tmp[5] = 0.0;
    c_T_tmp[8] = 1.0;
    R[0] = R_idx_0;
    R[1] = rowStart;
    R[6] = -R_idx_0 * T1[T1.size(0) * 2] + -colStart * T1[T1.size(0) * 2 + 1];
    R[3] = colStart;
    R[4] = R_idx_3;
    R[7] = -rowStart * T1[T1.size(0) * 2] + -R_idx_3 * T1[T1.size(0) * 2 + 1];
    R[2] = 0.0;
    R[5] = 0.0;
    R[8] = 1.0;
    for (i = 0; i < 3; i++) {
        R_idx_0 = c_T_tmp[i];
        rowStart = c_T_tmp[i + 3];
        colStart = c_T_tmp[i + 6];
        for (boffset = 0; boffset < 3; boffset++) {
            T1Offset[i + 3 * boffset] =
                (R_idx_0 * R[3 * boffset] + rowStart * R[3 * boffset + 1]) +
                colStart * R[3 * boffset + 2];
        }
    }
    R_idx_0 = (*poseGraphUpdated)->NumNodes;
    i = static_cast<int>(R_idx_0);
    for (int b_i{0}; b_i < i; b_i++) {
        int y_size_idx_1;
        posesUpdated->extractBlock(static_cast<double>(b_i) + 1.0, T1);
        loop_ub = T1.size(1);
        y_size_idx_1 = T1.size(1);
        for (int j{0}; j < loop_ub; j++) {
            coffset = j * 3;
            boffset = j * T1.size(0);
            for (int c_i{0}; c_i < 3; c_i++) {
                c_T_tmp[coffset + c_i] = (T1Offset[c_i] * T1[boffset] +
                                          T1Offset[c_i + 3] * T1[boffset + 1]) +
                                         T1Offset[c_i + 6] * T1[boffset + 2];
            }
        }
        rowStart =
            posesUpdated->BlockSize[0] * ((static_cast<double>(b_i) + 1.0) - 1.0) +
            1.0;
        R_idx_0 = (rowStart + posesUpdated->BlockSize[0]) - 1.0;
        if (rowStart > R_idx_0) {
            boffset = 1;
        } else {
            boffset = static_cast<int>(rowStart);
        }
        for (loop_ub = 0; loop_ub < y_size_idx_1; loop_ub++) {
            posesUpdated
                ->Matrix[(boffset + posesUpdated->Matrix.size(0) * loop_ub) - 1] =
                c_T_tmp[3 * loop_ub];
            posesUpdated->Matrix[boffset + posesUpdated->Matrix.size(0) * loop_ub] =
                c_T_tmp[3 * loop_ub + 1];
            posesUpdated
                ->Matrix[(boffset + posesUpdated->Matrix.size(0) * loop_ub) + 1] =
                c_T_tmp[3 * loop_ub + 2];
        }
    }
    varargin_1.set_size(posesUpdated->Matrix.size(0), 3);
    loop_ub = posesUpdated->Matrix.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        varargin_1[i] = posesUpdated->Matrix[i];
    }
    loop_ub = varargin_1.size(0);
    for (i = 0; i < 3; i++) {
        for (boffset = 0; boffset < loop_ub; boffset++) {
            (*poseGraphUpdated)
                ->NodeEstimates
                ->Matrix[boffset +
                         (*poseGraphUpdated)->NodeEstimates->Matrix.size(0) * i] =
                varargin_1[boffset + varargin_1.size(0) * i];
        }
    }
}

}  // namespace internal
}  // namespace algs
}  // namespace nav
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for PoseGraphOptimizer.cpp
///
/// [EOF]
///
