///
/// @file           : CXSparseAPI.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "CXSparseAPI.h"
#include "rt_nonfinite.h"
#include "sparse1.h"
#include "coder_array.h"
#include "cs.h"
#include "makeCXSparseMatrix.h"
#include "solve_from_qr.h"

/// Function Definitions
///
/// @fn             : iteratedQR
/// @brief          :
/// @param          : const sparse &A
///                   const ::coder::array<double, 1U> &b
///                   int n
///                   ::coder::array<double, 1U> &out
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
void CXSparseAPI::iteratedQR(const sparse &A,
                             const ::coder::array<double, 1U> &b, int n,
                             ::coder::array<double, 1U> &out) {
    cs_di *cxA;
    cs_din *N;
    cs_dis *S;
    sparse in;
    ::coder::array<double, 1U> outBuff;
    double tol;
    int loop_ub;
    if (A.m < A.n) {
        A.ctranspose(in);
        cxA = makeCXSparseMatrix(in.colidx[in.colidx.size(0) - 1] - 1, in.n, in.m,
                                 &(in.colidx.data())[0], &(in.rowidx.data())[0],
                                 &(in.d.data())[0]);
    } else {
        cxA =
            makeCXSparseMatrix(A.colidx[A.colidx.size(0) - 1] - 1, A.n, A.m,
                               &(((::coder::array<int, 1U> *)&A.colidx)->data())[0],
                               &(((::coder::array<int, 1U> *)&A.rowidx)->data())[0],
                               &(((::coder::array<double, 1U> *)&A.d)->data())[0]);
    }
    S = cs_di_sqr(2, cxA, 1);
    N = cs_di_qr(cxA, S);
    cs_di_spfree(cxA);
    qr_rank_di(N, &tol);
    out.set_size(n);
    if (b.size(0) < n) {
        outBuff.set_size(n);
    } else {
        outBuff.set_size(b.size(0));
    }
    loop_ub = b.size(0);
    for (int i{0}; i < loop_ub; i++) {
        outBuff[i] = b[i];
    }
    solve_from_qr_di(N, S, (double *)&(outBuff.data())[0], b.size(0), n);
    if (n < 1) {
        loop_ub = 0;
    } else {
        loop_ub = n;
    }
    for (int i{0}; i < loop_ub; i++) {
        out[i] = outBuff[i];
    }
    cs_di_sfree(S);
    cs_di_nfree(N);
}

}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for CXSparseAPI.cpp
///
/// [EOF]
///
