///
/// @file           : poseGraph.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

/// @include file    : Include Files
#include "poseGraph.h"
#include "BlockMatrix.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_data.h"
#include "slamMapGraph2D_rtwutil.h"
#include "coder_array.h"
#include <cmath>
#include <cstring>

/// Function Definitions
///
/// @fn             : findEdgeID
/// @brief          :
/// @param          : const double nodePair[2]
///                   ::coder::array<double, 1U> &edgeID
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
void poseGraph::findEdgeID(const double nodePair[2],
                           ::coder::array<double, 1U> &edgeID) const
{
  ::coder::array<int, 1U> ii;
  ::coder::array<bool, 1U> tf;
  int iRowA;
  int idx;
  int nx;
  bool exitg1;
  tf.set_size(EdgeNodePairs.size(0));
  idx = EdgeNodePairs.size(0);
  for (iRowA = 0; iRowA < idx; iRowA++) {
    bool p;
    tf[iRowA] = false;
    p = true;
    nx = 0;
    exitg1 = false;
    while ((!exitg1) && (nx < 2)) {
      if (EdgeNodePairs[iRowA + EdgeNodePairs.size(0) * nx] != nodePair[nx]) {
        p = false;
        exitg1 = true;
      } else {
        nx++;
      }
    }
    if (p) {
      tf[iRowA] = true;
    }
  }
  nx = tf.size(0);
  idx = 0;
  ii.set_size(tf.size(0));
  iRowA = 0;
  exitg1 = false;
  while ((!exitg1) && (iRowA <= nx - 1)) {
    if (tf[iRowA]) {
      idx++;
      ii[idx - 1] = iRowA + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        iRowA++;
      }
    } else {
      iRowA++;
    }
  }
  if (tf.size(0) == 1) {
    if (idx == 0) {
      ii.set_size(0);
    }
  } else {
    if (idx < 1) {
      idx = 0;
    }
    ii.set_size(idx);
  }
  edgeID.set_size(ii.size(0));
  idx = ii.size(0);
  for (nx = 0; nx < idx; nx++) {
    edgeID[nx] = ii[nx];
  }
}

///
/// @fn             : addRelativePose
/// @brief          :
/// @param          : const double varargin_1[3]
///                   double varargin_3
///                   double varargin_4
/// @return         : void
///
void poseGraph::addRelativePose(const double varargin_1[3], double varargin_3,
                                double varargin_4)
{
  robotics::core::internal::BlockMatrix *obj;
  ::coder::array<double, 2U> T;
  ::coder::array<double, 2U> a;
  ::coder::array<double, 2U> s;
  ::coder::array<double, 1U> id;
  double Omega[9];
  double RR[9];
  double Trel[9];
  double nodePair[2];
  double T_tmp;
  double d;
  double edgeId;
  double fromNodeId;
  double toNodeId;
  int coffset;
  int i;
  int iRowS;
  int inner;
  int j;
  int k;
  bool constraintNeedsInversion;
  bool exitg1;
  bool isLoopClosure;
  bool needNewPoseNode;
  bool tf;
  constraintNeedsInversion = false;
  if ((!(varargin_3 == varargin_4)) &&
      ((!(varargin_3 <= NumNodes)) ||
       (!IsLandmarkNode[static_cast<int>(varargin_3) - 1])) &&
      ((!(varargin_4 <= NumNodes)) ||
       (!IsLandmarkNode[static_cast<int>(varargin_4) - 1]))) {
    if ((varargin_3 <= NumNodes) && (varargin_4 <= NumNodes)) {
      bool guard1{false};
      nodePair[0] = varargin_3;
      nodePair[1] = varargin_4;
      findEdgeID(nodePair, id);
      guard1 = false;
      if (id.size(0) == 0) {
        guard1 = true;
      } else {
        nodePair[0] = varargin_3;
        nodePair[1] = varargin_4;
        tf = false;
        iRowS = 0;
        exitg1 = false;
        while ((!exitg1) && (iRowS <= LoopClosureEdgeNodePairs.size(0) - 1)) {
          bool exitg2;
          needNewPoseNode = true;
          j = 0;
          exitg2 = false;
          while ((!exitg2) && (j < 2)) {
            if (nodePair[j] !=
                LoopClosureEdgeNodePairs[iRowS +
                                         LoopClosureEdgeNodePairs.size(0) *
                                             j]) {
              needNewPoseNode = false;
              exitg2 = true;
            } else {
              j++;
            }
          }
          if (needNewPoseNode) {
            tf = true;
            exitg1 = true;
          } else {
            iRowS++;
          }
        }
        if (tf) {
          guard1 = true;
        } else {
          fromNodeId = std::fmin(varargin_3, varargin_4);
          toNodeId = std::fmax(varargin_3, varargin_4);
          needNewPoseNode = false;
          constraintNeedsInversion = (fromNodeId != varargin_3);
        }
      }
      if (guard1) {
        fromNodeId = varargin_3;
        toNodeId = varargin_4;
        needNewPoseNode = false;
      }
    } else {
      edgeId = std::fmin(varargin_3, varargin_4);
      if ((edgeId <= NumNodes) &&
          (std::fmax(varargin_3, varargin_4) - NumNodes == 1.0)) {
        fromNodeId = edgeId;
        toNodeId = NumNodes + 1.0;
        needNewPoseNode = true;
        constraintNeedsInversion = (edgeId != varargin_3);
      }
    }
  }
  nodePair[0] = fromNodeId;
  nodePair[1] = toNodeId;
  findEdgeID(nodePair, id);
  isLoopClosure = false;
  if (!needNewPoseNode) {
    if (id.size(0) == 0) {
      isLoopClosure = true;
    } else {
      s.set_size(1, LoopClosureEdgeIDsInternal.size(1));
      iRowS = LoopClosureEdgeIDsInternal.size(1);
      for (i = 0; i < iRowS; i++) {
        s[i] = LoopClosureEdgeIDsInternal[i];
      }
      tf = false;
      k = 0;
      exitg1 = false;
      while ((!exitg1) && (k <= s.size(1) - 1)) {
        if (id[0] == s[k]) {
          tf = true;
          exitg1 = true;
        } else {
          k++;
        }
      }
      if (tf) {
        isLoopClosure = true;
      }
    }
  }
  edgeId = std::sin(varargin_1[2]);
  T_tmp = std::cos(varargin_1[2]);
  Trel[0] = T_tmp;
  Trel[3] = -edgeId;
  Trel[6] = varargin_1[0];
  Trel[1] = edgeId;
  Trel[4] = T_tmp;
  Trel[7] = varargin_1[1];
  Trel[2] = 0.0;
  Trel[5] = 0.0;
  Trel[8] = 1.0;
  for (i = 0; i < 9; i++) {
    Omega[i] = iv[i];
  }
  iRowS = static_cast<int>(NumEdges + 1.0);
  EdgeNodePairs[iRowS - 1] = fromNodeId;
  EdgeNodePairs[(iRowS + EdgeNodePairs.size(0)) - 1] = toNodeId;
  if (constraintNeedsInversion) {
    double b_RR[9];
    Trel[0] = T_tmp;
    Trel[1] = -edgeId;
    Trel[6] = T_tmp * -varargin_1[0] + edgeId * -varargin_1[1];
    Trel[3] = edgeId;
    Trel[4] = T_tmp;
    Trel[7] = -edgeId * -varargin_1[0] + T_tmp * -varargin_1[1];
    Trel[2] = 0.0;
    Trel[5] = 0.0;
    Trel[8] = 1.0;
    std::memset(&RR[0], 0, 9U * sizeof(double));
    RR[0] = T_tmp;
    RR[1] = edgeId;
    RR[3] = -edgeId;
    RR[4] = T_tmp;
    RR[8] = 1.0;
    for (i = 0; i < 3; i++) {
      d = RR[i];
      edgeId = RR[i + 3];
      inner = static_cast<int>(RR[i + 6]);
      for (coffset = 0; coffset < 3; coffset++) {
        b_RR[i + 3 * coffset] =
            (d * static_cast<double>(iv[3 * coffset]) +
             edgeId * static_cast<double>(iv[3 * coffset + 1])) +
            static_cast<double>(inner * iv[3 * coffset + 2]);
      }
      d = b_RR[i];
      edgeId = b_RR[i + 3];
      T_tmp = b_RR[i + 6];
      for (inner = 0; inner < 3; inner++) {
        Omega[i + 3 * inner] =
            (d * RR[inner] + edgeId * RR[inner + 3]) + T_tmp * RR[inner + 6];
      }
    }
  }
  EdgeMeasurements->replaceBlock(NumEdges + 1.0, Trel);
  EdgeInfoMatrices->replaceBlock(NumEdges + 1.0, Omega);
  NumEdges++;
  edgeId = NumEdges;
  if (needNewPoseNode) {
    NodeEstimates->extractBlock(fromNodeId, a);
    iRowS = a.size(0);
    inner = a.size(1);
    T.set_size(a.size(0), 3);
    for (j = 0; j < 3; j++) {
      int boffset;
      coffset = j * iRowS;
      boffset = j * 3;
      for (int b_i{0}; b_i < iRowS; b_i++) {
        T[coffset + b_i] = 0.0;
      }
      for (k = 0; k < inner; k++) {
        int aoffset;
        aoffset = k * a.size(0);
        edgeId = Trel[boffset + k];
        for (int b_i{0}; b_i < iRowS; b_i++) {
          i = coffset + b_i;
          T[i] = T[i] + a[aoffset + b_i] * edgeId;
        }
      }
    }
    obj = NodeEstimates;
    edgeId = NumNodes + 1.0;
    edgeId = obj->BlockSize[0] * (edgeId - 1.0) + 1.0;
    d = (edgeId + obj->BlockSize[0]) - 1.0;
    if (edgeId > d) {
      i = 1;
    } else {
      i = static_cast<int>(edgeId);
    }
    iRowS = T.size(0);
    for (inner = 0; inner < 3; inner++) {
      for (coffset = 0; coffset < iRowS; coffset++) {
        obj->Matrix[((i + coffset) + obj->Matrix.size(0) * inner) - 1] =
            T[coffset + T.size(0) * inner];
      }
    }
    NumNodes++;
    NodeDims[static_cast<int>(NumNodes) - 1] = 3.0;
    NodeMap[static_cast<int>(NumNodes) - 1] =
        NodeMap[static_cast<int>(NumNodes - 1.0) - 1] +
        NodeDims[static_cast<int>(NumNodes - 1.0) - 1];
    IsLandmarkNode[static_cast<int>(NumNodes) - 1] = false;
  } else if (isLoopClosure) {
    iRowS = static_cast<int>(NumLoopClosureEdges + 1.0);
    LoopClosureEdgeNodePairs[iRowS - 1] = fromNodeId;
    LoopClosureEdgeNodePairs[(iRowS + LoopClosureEdgeNodePairs.size(0)) - 1] =
        toNodeId;
    LoopClosureEdgeIDsInternal[static_cast<int>(NumLoopClosureEdges + 1.0) -
                               1] = edgeId;
    NumLoopClosureEdges++;
  }
}

///
/// @fn             : get_LandmarkNodeIDs
/// @brief          :
/// @param          : ::coder::array<double, 2U> &lmIds
/// @return         : void
///
void poseGraph::get_LandmarkNodeIDs(::coder::array<double, 2U> &lmIds) const
{
  ::coder::array<double, 2U> nodeIds;
  double b;
  int end;
  int loop_ub;
  b = NumNodes;
  if (std::isnan(b)) {
    nodeIds.set_size(1, 1);
    nodeIds[0] = rtNaN;
  } else if (b < 1.0) {
    nodeIds.set_size(1, 0);
  } else {
    nodeIds.set_size(1, static_cast<int>(b - 1.0) + 1);
    loop_ub = static_cast<int>(b - 1.0);
    for (end = 0; end <= loop_ub; end++) {
      nodeIds[end] = static_cast<double>(end) + 1.0;
    }
  }
  end = IsLandmarkNode.size(0) - 1;
  loop_ub = 0;
  for (int i{0}; i <= end; i++) {
    if (IsLandmarkNode[i]) {
      loop_ub++;
    }
  }
  lmIds.set_size(1, loop_ub);
  loop_ub = 0;
  for (int i{0}; i <= end; i++) {
    if (IsLandmarkNode[i]) {
      lmIds[loop_ub] = nodeIds[i];
      loop_ub++;
    }
  }
}

///
/// @fn             : get_LoopClosureEdgeIDs
/// @brief          :
/// @param          : ::coder::array<double, 2U> &lcIds
/// @return         : void
///
void poseGraph::get_LoopClosureEdgeIDs(::coder::array<double, 2U> &lcIds) const
{
  double d;
  int loop_ub;
  d = NumLoopClosureEdges;
  if (d < 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(d);
  }
  lcIds.set_size(1, loop_ub);
  for (int i{0}; i < loop_ub; i++) {
    lcIds[i] = LoopClosureEdgeIDsInternal[i];
  }
}

///
/// @fn             : init
/// @brief          :
/// @param          : double varargin_2
///                   double varargin_4
/// @return         : poseGraph *
///
poseGraph *poseGraph::init(double varargin_2, double varargin_4)
{
  poseGraph *obj;
  robotics::core::internal::BlockMatrix *b_obj;
  double colStart;
  double d;
  double rowStart;
  int i;
  int loop_ub;
  int obj_tmp;
  obj = this;
  obj->MaxNumEdges = varargin_2;
  obj->MaxNumNodes = varargin_4;
  obj->NumEdges = 0.0;
  obj->NumNodes = 1.0;
  obj->NumLoopClosureEdges = 0.0;
  rowStart = obj->MaxNumNodes;
  i = static_cast<int>(rowStart * 3.0);
  obj->_pobj0[0].Matrix.set_size(i, 3);
  loop_ub = i * 3;
  for (i = 0; i < loop_ub; i++) {
    obj->_pobj0[0].Matrix[i] = 0.0;
  }
  obj->_pobj0[0].BlockSize[0] = 3.0;
  obj->_pobj0[0].BlockSize[1] = 3.0;
  obj->_pobj0[0].NumRowBlocks = rowStart;
  obj->_pobj0[0].NumColBlocks = 1.0;
  obj->NodeEstimates = &obj->_pobj0[0];
  b_obj = obj->NodeEstimates;
  rowStart = b_obj->BlockSize[0] * 0.0 + 1.0;
  colStart = b_obj->BlockSize[1] * 0.0 + 1.0;
  d = (rowStart + b_obj->BlockSize[0]) - 1.0;
  if (rowStart > d) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(d);
  }
  d = (colStart + b_obj->BlockSize[1]) - 1.0;
  if (colStart > d) {
    obj_tmp = 0;
  } else {
    obj_tmp = static_cast<int>(d);
  }
  for (i = 0; i < obj_tmp; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      b_obj->Matrix[i1 + b_obj->Matrix.size(0) * i] = iv[i1 + loop_ub * i];
    }
  }
  loop_ub = static_cast<int>(obj->MaxNumNodes);
  obj->NodeMap.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    obj->NodeMap[i] = 0.0;
  }
  loop_ub = static_cast<int>(obj->MaxNumNodes);
  obj->NodeDims.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    obj->NodeDims[i] = 0.0;
  }
  loop_ub = static_cast<int>(obj->MaxNumNodes);
  obj->IsLandmarkNode.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    obj->IsLandmarkNode[i] = false;
  }
  obj->NodeMap[0] = 1.0;
  obj->NodeDims[0] = 3.0;
  loop_ub = static_cast<int>(obj->MaxNumEdges);
  obj->EdgeNodePairs.set_size(loop_ub, 2);
  loop_ub <<= 1;
  for (i = 0; i < loop_ub; i++) {
    obj->EdgeNodePairs[i] = 0.0;
  }
  loop_ub = static_cast<int>(obj->MaxNumEdges);
  obj->LoopClosureEdgeNodePairs.set_size(loop_ub, 2);
  loop_ub <<= 1;
  for (i = 0; i < loop_ub; i++) {
    obj->LoopClosureEdgeNodePairs[i] = 0.0;
  }
  loop_ub = static_cast<int>(obj->MaxNumEdges);
  obj->LoopClosureEdgeIDsInternal.set_size(1, loop_ub);
  for (i = 0; i < loop_ub; i++) {
    obj->LoopClosureEdgeIDsInternal[i] = 0.0;
  }
  rowStart = obj->MaxNumEdges;
  i = static_cast<int>(rowStart * 3.0);
  obj->_pobj0[1].Matrix.set_size(i, 3);
  loop_ub = i * 3;
  for (i = 0; i < loop_ub; i++) {
    obj->_pobj0[1].Matrix[i] = 0.0;
  }
  obj->_pobj0[1].BlockSize[0] = 3.0;
  obj->_pobj0[1].BlockSize[1] = 3.0;
  obj->_pobj0[1].NumRowBlocks = rowStart;
  obj->_pobj0[1].NumColBlocks = 1.0;
  obj->EdgeMeasurements = &obj->_pobj0[1];
  rowStart = obj->MaxNumEdges;
  i = static_cast<int>(rowStart * 3.0);
  obj->_pobj0[2].Matrix.set_size(i, 3);
  loop_ub = i * 3;
  for (i = 0; i < loop_ub; i++) {
    obj->_pobj0[2].Matrix[i] = 0.0;
  }
  obj->_pobj0[2].BlockSize[0] = 3.0;
  obj->_pobj0[2].BlockSize[1] = 3.0;
  obj->_pobj0[2].NumRowBlocks = rowStart;
  obj->_pobj0[2].NumColBlocks = 1.0;
  obj->EdgeInfoMatrices = &obj->_pobj0[2];
  return obj;
}

///
/// @fn             : nodeEstimates
/// @brief          :
/// @param          : ::coder::array<double, 2U> &nodeEsts
/// @return         : void
///
void poseGraph::nodeEstimates(::coder::array<double, 2U> &nodeEsts) const
{
  ::coder::array<double, 2U> T;
  ::coder::array<double, 2U> nodeIds;
  double L;
  int i;
  int loop_ub;
  L = NumNodes;
  if (std::isnan(L)) {
    nodeIds.set_size(1, 1);
    nodeIds[0] = rtNaN;
  } else if (L < 1.0) {
    nodeIds.set_size(1, 0);
  } else {
    nodeIds.set_size(1, static_cast<int>(L - 1.0) + 1);
    loop_ub = static_cast<int>(L - 1.0);
    for (i = 0; i <= loop_ub; i++) {
      nodeIds[i] = static_cast<double>(i) + 1.0;
    }
  }
  i = static_cast<int>(L);
  nodeEsts.set_size(i, 3);
  for (loop_ub = 0; loop_ub < i; loop_ub++) {
    L = nodeIds[loop_ub];
    if (IsLandmarkNode[static_cast<int>(L) - 1]) {
      NodeEstimates->extractBlock(L, T);
      nodeEsts[loop_ub] = T[T.size(0) * 2];
      nodeEsts[loop_ub + nodeEsts.size(0)] = T[T.size(0) * 2 + 1];
      nodeEsts[loop_ub + nodeEsts.size(0) * 2] = rtNaN;
    } else {
      NodeEstimates->extractBlock(L, T);
      nodeEsts[loop_ub] = T[T.size(0) * 2];
      nodeEsts[loop_ub + nodeEsts.size(0)] = T[T.size(0) * 2 + 1];
      nodeEsts[loop_ub + nodeEsts.size(0) * 2] = rt_atan2d_snf(T[1], T[0]);
    }
  }
}

} // namespace coder
} // namespace SlamGraph2D

///
/// File trailer for poseGraph.cpp
///
/// [EOF]
///
