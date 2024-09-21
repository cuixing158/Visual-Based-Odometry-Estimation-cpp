///
/// @file           : TrustRegionIndefiniteDogLegSE2.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "TrustRegionIndefiniteDogLegSE2.h"
#include "BlockMatrix.h"
#include "CXSparseAPI.h"
#include "PoseGraphHelpers.h"
#include "SystemTimeProvider.h"
#include "TrustRegionIndefiniteDogLegInterface.h"
#include "myGraph.h"
#include "norm.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_internal_types.h"
#include "slamMapGraph2D_rtwutil.h"
#include "sparse1.h"
#include "tic.h"
#include "toc.h"
#include "coder_array.h"
#include "coder_posix_time.h"
#include "cs.h"
#include "makeCXSparseMatrix.h"
#include "solve_from_lu.h"
#include <cmath>

/// Type Definitions
namespace SlamGraph2D {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
enum class NLPSolverExitFlags : int {
    LocalMinimumFound = 1,  // Default value
    IterationLimitExceeded,
    TimeLimitExceeded,
    StepSizeBelowMinimum,
    ChangeInErrorBelowMinimum,
    SearchDirectionInvalid,
    HessianNotPositiveSemidefinite,
    TrustRegionRadiusBelowMinimum
};

}
}  // namespace core
}  // namespace robotics
}  // namespace coder
}  // namespace SlamGraph2D

/// Function Declarations
namespace SlamGraph2D {
static void binary_expand_op(::coder::array<double, 2U> &in1,
                             const ::coder::array<double, 1U> &in2, int in3,
                             int in4);

static void binary_expand_op(double in1[3],
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 1U> &in3, int in4,
                             int in5);

}  // namespace SlamGraph2D

/// Function Definitions
///
/// @fn             : computeBasicSteps
/// @brief          :
/// @param          : const ::coder::array<double, 1U> &grad
///                   const sparse &B
///                   ::coder::array<double, 1U> &stepSD
///                   ::coder::array<double, 1U> &stepGN
/// @return         : bool
///
namespace SlamGraph2D {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
bool TrustRegionIndefiniteDogLegSE2::computeBasicSteps(
    const ::coder::array<double, 1U> &grad, const sparse &B,
    ::coder::array<double, 1U> &stepSD,
    ::coder::array<double, 1U> &stepGN) const {
    sparse in;
    ::coder::array<double, 2U> a;
    ::coder::array<double, 1U> b_stepGN;
    double b_grad;
    double cd;
    int apend;
    int i;
    int loop_ub_tmp;
    a.set_size(1, B.n);
    loop_ub_tmp = B.n;
    for (i = 0; i < loop_ub_tmp; i++) {
        a[i] = 0.0;
    }
    if ((grad.size(0) != 0) && (B.n != 0) &&
        (B.colidx[B.colidx.size(0) - 1] - 1 != 0)) {
        for (int k{0}; k < loop_ub_tmp; k++) {
            cd = 0.0;
            apend = B.colidx[k + 1] - 1;
            i = B.colidx[k];
            for (int ap{i}; ap <= apend; ap++) {
                cd += B.d[ap - 1] * grad[B.rowidx[ap - 1] - 1];
            }
            a[k] = cd;
        }
    }
    b_grad = 0.0;
    apend = grad.size(0);
    for (i = 0; i < apend; i++) {
        cd = grad[i];
        b_grad += cd * cd;
    }
    cd = 0.0;
    apend = a.size(1);
    for (i = 0; i < apend; i++) {
        cd += a[i] * grad[i];
    }
    cd = b_grad / cd;
    stepGN.set_size(grad.size(0));
    apend = grad.size(0);
    for (i = 0; i < apend; i++) {
        stepGN[i] = -grad[i];
    }
    stepSD.set_size(stepGN.size(0));
    apend = stepGN.size(0);
    for (i = 0; i < apend; i++) {
        stepSD[i] = cd * stepGN[i];
    }
    if ((B.m == 0) || (B.n == 0) || (stepGN.size(0) == 0)) {
        stepGN.set_size(B.n);
        for (i = 0; i < loop_ub_tmp; i++) {
            stepGN[i] = 0.0;
        }
    } else if (stepGN.size(0) == B.n) {
        cs_di *cxA;
        cs_din *N;
        cs_dis *S;
        if (B.m < B.n) {
            B.ctranspose(in);
            cxA = makeCXSparseMatrix(in.colidx[in.colidx.size(0) - 1] - 1, in.n, in.m,
                                     &(in.colidx.data())[0], &(in.rowidx.data())[0],
                                     &(in.d.data())[0]);
        } else {
            cxA = makeCXSparseMatrix(
                B.colidx[B.colidx.size(0) - 1] - 1, B.n, B.m,
                &(((::coder::array<int, 1U> *)&B.colidx)->data())[0],
                &(((::coder::array<int, 1U> *)&B.rowidx)->data())[0],
                &(((::coder::array<double, 1U> *)&B.d)->data())[0]);
        }
        S = cs_di_sqr(2, cxA, 0);
        N = cs_di_lu(cxA, S, 1);
        cs_di_spfree(cxA);
        if (N == nullptr) {
            cs_di_sfree(S);
            cs_di_nfree(N);
            b_stepGN.set_size(stepGN.size(0));
            apend = stepGN.size(0) - 1;
            for (i = 0; i <= apend; i++) {
                b_stepGN[i] = stepGN[i];
            }
            ::SlamGraph2D::coder::internal::CXSparseAPI::iteratedQR(B, b_stepGN, B.n,
                                                                    stepGN);
        } else {
            solve_from_lu_di(N, S, (double *)&(stepGN.data())[0], stepGN.size(0));
            cs_di_sfree(S);
            cs_di_nfree(N);
        }
    } else {
        b_stepGN.set_size(stepGN.size(0));
        apend = stepGN.size(0) - 1;
        for (i = 0; i <= apend; i++) {
            b_stepGN[i] = stepGN[i];
        }
        ::SlamGraph2D::coder::internal::CXSparseAPI::iteratedQR(B, b_stepGN, B.n,
                                                                stepGN);
    }
    return b_norm(grad) < GradientTolerance;
}

///
/// @fn             : incrementX
/// @brief          :
/// @param          : const ::coder::array<double, 2U> &x
///                   const ::coder::array<double, 1U> &epsilons
///                   ::coder::array<double, 2U> &xNew
/// @return         : void
///
void TrustRegionIndefiniteDogLegSE2::incrementX(
    const ::coder::array<double, 2U> &x,
    const ::coder::array<double, 1U> &epsilons,
    ::coder::array<double, 2U> &xNew) const {
    BlockMatrix xBlk;
    b_BlockMatrix xBlkNew;
    ::coder::array<double, 2U> T;
    double b_T[3];
    int i;
    int i1;
    int loop_ub;
    xBlk.Matrix.set_size(x.size(0), 3);
    loop_ub = x.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        xBlk.Matrix[i] = x[i];
    }
    xBlk.BlockSize[0] = 3.0;
    xBlk.BlockSize[1] = 3.0;
    xBlk.NumRowBlocks = static_cast<double>(x.size(0)) / 3.0;
    xBlk.NumColBlocks = 1.0;
    i = static_cast<int>(xBlk.NumRowBlocks * 3.0);
    xBlkNew.Matrix.set_size(i, 3);
    loop_ub = i * 3;
    for (i = 0; i < loop_ub; i++) {
        xBlkNew.Matrix[i] = 0.0;
    }
    i = static_cast<int>(xBlk.NumRowBlocks);
    for (int b_i{0}; b_i < i; b_i++) {
        double idx;
        double len;
        xBlk.extractBlock(static_cast<double>(b_i) + 1.0, T);
        idx = ExtraArgs.nodeMap[b_i];
        len = ExtraArgs.nodeDims[b_i];
        if (len == 3.0) {
            double T_tmp[9];
            int i2;
            int rowStart;
            if (idx > (idx + 3.0) - 1.0) {
                i1 = 0;
                i2 = -2;
            } else {
                i1 = static_cast<int>(idx) - 1;
                i2 = static_cast<int>(static_cast<unsigned int>(idx));
            }
            if ((i2 - i1) + 2 == 3) {
                b_T[0] = T[T.size(0) * 2];
                b_T[1] = T[T.size(0) * 2 + 1];
                b_T[2] = rt_atan2d_snf(T[1], T[0]);
                b_T[0] += epsilons[i1];
                b_T[1] += epsilons[i1 + 1];
                b_T[2] += epsilons[i1 + 2];
            } else {
                binary_expand_op(b_T, T, epsilons, i1, i2 + 1);
            }
            len = std::sin(b_T[2]);
            idx = std::cos(b_T[2]);
            rowStart = 3 * b_i + 1;
            if (static_cast<unsigned int>(rowStart) >
                static_cast<unsigned int>(rowStart) + 2U) {
                i1 = 0;
                i2 = 0;
            } else {
                i1 = rowStart - 1;
                i2 = static_cast<int>(static_cast<unsigned int>(rowStart) + 2U);
            }
            T_tmp[0] = idx;
            T_tmp[3] = -len;
            T_tmp[6] = b_T[0];
            T_tmp[1] = len;
            T_tmp[4] = idx;
            T_tmp[7] = b_T[1];
            T_tmp[2] = 0.0;
            T_tmp[5] = 0.0;
            T_tmp[8] = 1.0;
            loop_ub = i2 - i1;
            for (i2 = 0; i2 < 3; i2++) {
                for (rowStart = 0; rowStart < loop_ub; rowStart++) {
                    xBlkNew.Matrix[(i1 + rowStart) + xBlkNew.Matrix.size(0) * i2] =
                        T_tmp[rowStart + loop_ub * i2];
                }
            }
        } else {
            int i2;
            int rowStart;
            len = (idx + len) - 1.0;
            if (idx > len) {
                i1 = 0;
                i2 = 0;
            } else {
                i1 = static_cast<int>(idx) - 1;
                i2 = static_cast<int>(len);
            }
            if (i2 - i1 == 2) {
                T[T.size(0) * 2] = T[T.size(0) * 2] + epsilons[i1];
                T[T.size(0) * 2 + 1] = T[T.size(0) * 2 + 1] + epsilons[i1 + 1];
            } else {
                binary_expand_op(T, epsilons, i1, i2 - 1);
            }
            rowStart = 3 * b_i;
            if (static_cast<unsigned int>(rowStart + 1) >
                static_cast<unsigned int>(rowStart) + 3U) {
                rowStart = 0;
            }
            loop_ub = T.size(1);
            for (i1 = 0; i1 < loop_ub; i1++) {
                int b_loop_ub;
                b_loop_ub = T.size(0);
                for (i2 = 0; i2 < b_loop_ub; i2++) {
                    xBlkNew.Matrix[(rowStart + i2) + xBlkNew.Matrix.size(0) * i1] =
                        T[i2 + T.size(0) * i1];
                }
            }
        }
    }
    xNew.set_size(xBlkNew.Matrix.size(0), 3);
    loop_ub = xBlkNew.Matrix.size(0);
    for (i = 0; i < 3; i++) {
        for (i1 = 0; i1 < loop_ub; i1++) {
            xNew[i1 + xNew.size(0) * i] =
                xBlkNew.Matrix[i1 + xBlkNew.Matrix.size(0) * i];
        }
    }
}

///
/// @fn             : binary_expand_op
/// @brief          :
/// @param          : ::coder::array<double, 2U> &in1
///                   const ::coder::array<double, 1U> &in2
///                   int in3
///                   int in4
/// @return         : void
///
}  // namespace internal
}  // namespace core
}  // namespace robotics
}  // namespace coder
static void binary_expand_op(::coder::array<double, 2U> &in1,
                             const ::coder::array<double, 1U> &in2, int in3,
                             int in4) {
    in1[in1.size(0) * 2] = in1[in1.size(0) * 2] + in2[in3];
    in1[in1.size(0) * 2 + 1] =
        in1[in1.size(0) * 2 + 1] + in2[in3 + ((in4 - in3) + 1 != 1)];
}

///
/// @fn             : binary_expand_op
/// @brief          :
/// @param          : double in1[3]
///                   const ::coder::array<double, 2U> &in2
///                   const ::coder::array<double, 1U> &in3
///                   int in4
///                   int in5
/// @return         : void
///
static void binary_expand_op(double in1[3],
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 1U> &in3, int in4,
                             int in5) {
    int stride_0_0;
    in1[0] = in2[in2.size(0) * 2];
    in1[1] = in2[in2.size(0) * 2 + 1];
    in1[2] = rt_atan2d_snf(in2[1], in2[0]);
    stride_0_0 = ((in5 - in4) + 1 != 1);
    in1[0] += in3[in4];
    in1[1] += in3[in4 + stride_0_0];
    in1[2] += in3[in4 + (stride_0_0 << 1)];
}

///
/// @fn             : solve
/// @brief          :
/// @param          : myGraph *aInstancePtr
///                   const ::coder::array<double, 2U> &seed
///                   BlockMatrix &iobj_0
///                   BlockMatrix **xSol
///                   sparse &hess
///                   double &solutionInfo_Error
///                   double &solutionInfo_ExitFlag
/// @return         : double
///
namespace coder {
namespace robotics {
namespace core {
namespace internal {
double TrustRegionIndefiniteDogLegSE2::solve(
    myGraph *aInstancePtr, const ::coder::array<double, 2U> &seed,
    BlockMatrix &iobj_0, BlockMatrix **xSol, sparse &hess,
    double &solutionInfo_Error, double &solutionInfo_ExitFlag) {
    BlockMatrix *b_xSol;
    sparse HessianTrial;
    sparse c;
    ::coder::array<double, 2U> x;
    ::coder::array<double, 2U> xn;
    ::coder::array<double, 1U> grad;
    ::coder::array<double, 1U> gradTrial;
    ::coder::array<double, 1U> stepDL_;
    ::coder::array<double, 1U> stepGN;
    ::coder::array<double, 1U> stepSD;
    ::coder::array<double, 1U> tmpd;
    ::coder::array<bool, 1U> b_x;
    double delta;
    double solutionInfo_Iterations;
    int i;
    int nx;
    bool localMin;
    NLPSolverExitFlags exitFlag;
    MaxNumIterationInternal = MaxNumIteration;
    MaxTimeInternal = MaxTime;
    SeedInternal.set_size(seed.size(0), 3);
    nx = seed.size(0) * 3;
    for (i = 0; i < nx; i++) {
        SeedInternal[i] = seed[i];
    }
    TimeObj.StartTime.tv_sec = tic(aInstancePtr, TimeObj.StartTime.tv_nsec);
    x.set_size(SeedInternal.size(0), 3);
    nx = SeedInternal.size(0) * 3;
    for (i = 0; i < nx; i++) {
        x[i] = SeedInternal[i];
    }
    (&(&iobj_0)[0])[0].Matrix.set_size(SeedInternal.size(0), 3);
    nx = SeedInternal.size(0) * 3;
    for (i = 0; i < nx; i++) {
        (&(&iobj_0)[0])[0].Matrix[i] = SeedInternal[i];
    }
    (&(&iobj_0)[0])[0].BlockSize[0] = 3.0;
    (&(&iobj_0)[0])[0].BlockSize[1] = 3.0;
    (&(&iobj_0)[0])[0].NumRowBlocks =
        static_cast<double>(SeedInternal.size(0)) / 3.0;
    (&(&iobj_0)[0])[0].NumColBlocks = 1.0;
    b_xSol = &(&(&iobj_0)[0])[0];
    solutionInfo_Iterations = 0.0;
    exitFlag = NLPSolverExitFlags::IterationLimitExceeded;
    solutionInfo_Error = nav::algs::internal::PoseGraphHelpers::poseGraphCost(
        SeedInternal, ExtraArgs.edgeNodePairs, ExtraArgs.edgeMeasurements,
        ExtraArgs.edgeInfoMats, ExtraArgs.tformSize, ExtraArgs.infoMatSize,
        ExtraArgs.poseDeltaLength, ExtraArgs.nodeMap, ExtraArgs.nodeDims,
        ExtraArgs.IsLandmarkNode, grad, hess);
    delta = InitialTrustRegionRadius;
    tmpd.set_size(grad.size(0));
    nx = grad.size(0);
    for (i = 0; i < nx; i++) {
        tmpd[i] = grad[i];
    }
    localMin = computeBasicSteps(tmpd, hess, stepSD, stepGN);
    if (localMin) {
        exitFlag = NLPSolverExitFlags::LocalMinimumFound;
    } else {
        double d;
        int b_i;
        bool exitg1;
        bool terminated;
        terminated = false;
        d = MaxNumIterationInternal;
        b_i = 0;
        exitg1 = false;
        while ((!exitg1) && (b_i <= static_cast<int>(d) - 1)) {
            double val;
            val = toc(aInstancePtr, TimeObj.StartTime.tv_sec,
                      TimeObj.StartTime.tv_nsec);
            if (val > MaxTimeInternal) {
                exitFlag = NLPSolverExitFlags::TimeLimitExceeded;
                terminated = true;
                exitg1 = true;
            } else {
                double b_stepDL_;
                double bc;
                int currRowIdx;
                bool exitg2;
                val = b_norm(stepSD);
                if (val >= delta) {
                    val = delta / val;
                    stepDL_.set_size(stepSD.size(0));
                    nx = stepSD.size(0);
                    for (i = 0; i < nx; i++) {
                        stepDL_[i] = val * stepSD[i];
                    }
                } else if (b_norm(stepGN) <= delta) {
                    stepDL_.set_size(stepGN.size(0));
                    nx = stepGN.size(0);
                    for (i = 0; i < nx; i++) {
                        stepDL_[i] = stepGN[i];
                    }
                } else {
                    if (stepGN.size(0) == stepSD.size(0)) {
                        stepDL_.set_size(stepGN.size(0));
                        nx = stepGN.size(0);
                        for (i = 0; i < nx; i++) {
                            stepDL_[i] = stepGN[i] - stepSD[i];
                        }
                    } else {
                        minus(stepDL_, stepGN, stepSD);
                    }
                    bc = 0.0;
                    nx = stepSD.size(0);
                    for (i = 0; i < nx; i++) {
                        bc += stepSD[i] * stepDL_[i];
                    }
                    b_stepDL_ = 0.0;
                    nx = stepDL_.size(0);
                    for (i = 0; i < nx; i++) {
                        b_stepDL_ += stepDL_[i] * stepDL_[i];
                    }
                    val = (-bc +
                           std::sqrt(bc * bc + b_stepDL_ * (delta * delta - val * val))) /
                          b_stepDL_;
                    if (stepSD.size(0) == stepDL_.size(0)) {
                        stepDL_.set_size(stepSD.size(0));
                        nx = stepSD.size(0);
                        for (i = 0; i < nx; i++) {
                            stepDL_[i] = stepSD[i] + val * stepDL_[i];
                        }
                    } else {
                        binary_expand_op(stepDL_, stepSD, val);
                    }
                }
                nx = stepDL_.size(0);
                tmpd.set_size(stepDL_.size(0));
                for (currRowIdx = 0; currRowIdx < nx; currRowIdx++) {
                    tmpd[currRowIdx] = std::abs(stepDL_[currRowIdx]);
                }
                b_x.set_size(tmpd.size(0));
                nx = tmpd.size(0);
                for (i = 0; i < nx; i++) {
                    b_x[i] = (tmpd[i] < StepTolerance);
                }
                localMin = true;
                nx = 1;
                exitg2 = false;
                while ((!exitg2) && (nx <= b_x.size(0))) {
                    if (!b_x[nx - 1]) {
                        localMin = false;
                        exitg2 = true;
                    } else {
                        nx++;
                    }
                }
                if (localMin) {
                    exitFlag = NLPSolverExitFlags::StepSizeBelowMinimum;
                    terminated = true;
                    exitg1 = true;
                } else {
                    double costTrial;
                    double d1;
                    incrementX(x, stepDL_, xn);
                    costTrial = nav::algs::internal::PoseGraphHelpers::poseGraphCost(
                        xn, ExtraArgs.edgeNodePairs, ExtraArgs.edgeMeasurements,
                        ExtraArgs.edgeInfoMats, ExtraArgs.tformSize,
                        ExtraArgs.infoMatSize, ExtraArgs.poseDeltaLength,
                        ExtraArgs.nodeMap, ExtraArgs.nodeDims, ExtraArgs.IsLandmarkNode,
                        gradTrial, HessianTrial);
                    d1 = solutionInfo_Error - costTrial;
                    if (std::abs(d1) < FunctionTolerance) {
                        exitFlag = NLPSolverExitFlags::ChangeInErrorBelowMinimum;
                        terminated = true;
                        exitg1 = true;
                    } else {
                        int b_c;
                        int ridx;
                        bool guard1{false};
                        ridx = hess.colidx[hess.colidx.size(0) - 1];
                        if (ridx - 1 < 1) {
                            nx = 0;
                        } else {
                            nx = ridx - 1;
                        }
                        tmpd.set_size(nx);
                        for (i = 0; i < nx; i++) {
                            tmpd[i] = 0.5 * hess.d[i];
                        }
                        if (ridx - 1 >= 1) {
                            nx = ridx - 2;
                        } else {
                            nx = 0;
                        }
                        c.d.set_size(nx + 1);
                        c.rowidx.set_size(nx + 1);
                        for (i = 0; i <= nx; i++) {
                            c.d[i] = 0.0;
                            c.rowidx[i] = 0;
                        }
                        if (ridx - 1 < 1) {
                            nx = 1;
                        } else {
                            nx = ridx;
                        }
                        for (i = 0; i <= nx - 2; i++) {
                            c.rowidx[i] = hess.rowidx[i];
                        }
                        c.colidx.set_size(hess.colidx.size(0));
                        nx = hess.colidx.size(0);
                        for (i = 0; i < nx; i++) {
                            c.colidx[i] = hess.colidx[i];
                        }
                        for (currRowIdx = 0; currRowIdx <= ridx - 2; currRowIdx++) {
                            c.d[currRowIdx] = tmpd[currRowIdx];
                        }
                        nx = 1;
                        i = hess.colidx.size(0);
                        for (b_c = 0; b_c <= i - 2; b_c++) {
                            ridx = c.colidx[b_c];
                            c.colidx[b_c] = nx;
                            while (ridx < c.colidx[b_c + 1]) {
                                currRowIdx = c.rowidx[ridx - 1];
                                val = c.d[ridx - 1];
                                ridx++;
                                if (val != 0.0) {
                                    c.d[nx - 1] = val;
                                    c.rowidx[nx - 1] = currRowIdx;
                                    nx++;
                                }
                            }
                        }
                        c.colidx[c.colidx.size(0) - 1] = nx;
                        tmpd.set_size(hess.m);
                        nx = hess.m;
                        for (i = 0; i < nx; i++) {
                            tmpd[i] = 0.0;
                        }
                        if ((hess.n != 0) && (hess.m != 0) &&
                            (c.colidx[c.colidx.size(0) - 1] - 1 != 0)) {
                            i = hess.n;
                            for (int acol{0}; acol < i; acol++) {
                                bc = stepDL_[acol];
                                b_c = c.colidx[acol];
                                ridx = c.colidx[acol + 1];
                                nx = ridx - c.colidx[acol];
                                if (nx >= 4) {
                                    int tmpd_tmp;
                                    currRowIdx = (ridx - nx) + ((nx / 4) << 2);
                                    for (int ap{b_c}; ap <= currRowIdx - 1; ap += 4) {
                                        tmpd_tmp = c.rowidx[ap - 1] - 1;
                                        tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap - 1] * bc;
                                        tmpd[c.rowidx[ap] - 1] =
                                            tmpd[c.rowidx[ap] - 1] + c.d[ap] * bc;
                                        tmpd_tmp = c.rowidx[ap + 1] - 1;
                                        tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap + 1] * bc;
                                        tmpd_tmp = c.rowidx[ap + 2] - 1;
                                        tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap + 2] * bc;
                                    }
                                    nx = ridx - 1;
                                    for (int ap{currRowIdx}; ap <= nx; ap++) {
                                        tmpd_tmp = c.rowidx[ap - 1] - 1;
                                        tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap - 1] * bc;
                                    }
                                } else {
                                    nx = ridx - 1;
                                    for (int ap{b_c}; ap <= nx; ap++) {
                                        int tmpd_tmp;
                                        tmpd_tmp = c.rowidx[ap - 1] - 1;
                                        tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap - 1] * bc;
                                    }
                                }
                            }
                        }
                        if (grad.size(0) == tmpd.size(0)) {
                            b_stepDL_ = 0.0;
                            nx = stepDL_.size(0);
                            for (i = 0; i < nx; i++) {
                                b_stepDL_ += -stepDL_[i] * (grad[i] + tmpd[i]);
                            }
                            val = d1 / b_stepDL_;
                        } else {
                            val = binary_expand_op(solutionInfo_Error, costTrial, stepDL_,
                                                   grad, tmpd);
                        }
                        guard1 = false;
                        if (val > 0.0) {
                            xn.set_size(x.size(0), 3);
                            nx = x.size(0) * x.size(1) - 1;
                            for (i = 0; i <= nx; i++) {
                                xn[i] = x[i];
                            }
                            incrementX(xn, stepDL_, x);
                            solutionInfo_Iterations = static_cast<double>(b_i) + 1.0;
                            solutionInfo_Error = costTrial;
                            grad.set_size(gradTrial.size(0));
                            nx = gradTrial.size(0);
                            for (i = 0; i < nx; i++) {
                                grad[i] = gradTrial[i];
                            }
                            hess = HessianTrial;
                            localMin =
                                computeBasicSteps(gradTrial, HessianTrial, stepSD, stepGN);
                            if (localMin) {
                                exitFlag = NLPSolverExitFlags::LocalMinimumFound;
                                terminated = true;
                                exitg1 = true;
                            } else {
                                guard1 = true;
                            }
                        } else {
                            guard1 = true;
                        }
                        if (guard1) {
                            bool b_guard1{false};
                            localMin = false;
                            b_guard1 = false;
                            if (val > 0.75) {
                                d1 = b_norm(stepDL_);
                                if (d1 > 0.9 * delta) {
                                    delta = 2.0 * d1;
                                } else {
                                    b_guard1 = true;
                                }
                            } else {
                                b_guard1 = true;
                            }
                            if (b_guard1 && (val < 0.25)) {
                                delta /= 4.0;
                                if (delta < TrustRegionRadiusTolerance) {
                                    localMin = true;
                                }
                            }
                            if (localMin) {
                                exitFlag = NLPSolverExitFlags::TrustRegionRadiusBelowMinimum;
                                terminated = true;
                                exitg1 = true;
                            } else {
                                b_i++;
                            }
                        }
                    }
                }
            }
        }
        b_xSol = &(&(&iobj_0)[0])[1];
        (&(&iobj_0)[0])[1].Matrix.set_size(x.size(0), 3);
        nx = x.size(0) * 3;
        for (i = 0; i < nx; i++) {
            (&(&iobj_0)[0])[1].Matrix[i] = x[i];
        }
        (&(&iobj_0)[0])[1].BlockSize[0] = 3.0;
        (&(&iobj_0)[0])[1].BlockSize[1] = 3.0;
        (&(&iobj_0)[0])[1].NumRowBlocks = static_cast<double>(x.size(0)) / 3.0;
        (&(&iobj_0)[0])[1].NumColBlocks = 1.0;
        if (!terminated) {
            solutionInfo_Iterations = d;
        }
    }
    *xSol = b_xSol;
    solutionInfo_ExitFlag = static_cast<double>(exitFlag);
    return solutionInfo_Iterations;
}

}  // namespace internal
}  // namespace core
}  // namespace robotics
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for TrustRegionIndefiniteDogLegSE2.cpp
///
/// [EOF]
///
