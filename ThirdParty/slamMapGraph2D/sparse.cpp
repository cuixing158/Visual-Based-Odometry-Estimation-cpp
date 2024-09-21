///
/// @file           : sparse.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "sparse.h"
#include "anonymous_function.h"
#include "introsort.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_internal_types.h"
#include "sparse1.h"
#include "coder_array.h"

/// Function Definitions
///
/// @fn             : b_sparse
/// @brief          :
/// @param          : const ::coder::array<double, 1U> &varargin_1
///                   const ::coder::array<double, 1U> &varargin_2
///                   const ::coder::array<double, 1U> &varargin_3
///                   sparse &y
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
void b_sparse(const ::coder::array<double, 1U> &varargin_1,
              const ::coder::array<double, 1U> &varargin_2,
              const ::coder::array<double, 1U> &varargin_3, sparse &y) {
    anonymous_function b_this;
    ::coder::array<int, 1U> sortedIndices;
    ::coder::array<int, 1U> t;
    int i;
    int nc;
    int ns;
    int ny;
    nc = varargin_2.size(0);
    ns = varargin_1.size(0);
    b_this.workspace.b.set_size(varargin_1.size(0));
    for (int k{0}; k < ns; k++) {
        b_this.workspace.b[k] = static_cast<int>(varargin_1[k]);
    }
    ns = varargin_2.size(0);
    b_this.workspace.a.set_size(varargin_2.size(0));
    sortedIndices.set_size(varargin_2.size(0));
    for (int k{0}; k < ns; k++) {
        b_this.workspace.a[k] = static_cast<int>(varargin_2[k]);
        sortedIndices[k] = k + 1;
    }
    internal::introsort(sortedIndices, b_this.workspace.a.size(0), b_this);
    ny = b_this.workspace.a.size(0);
    t.set_size(b_this.workspace.a.size(0));
    ns = b_this.workspace.a.size(0);
    for (i = 0; i < ns; i++) {
        t[i] = b_this.workspace.a[i];
    }
    for (int k{0}; k < ny; k++) {
        b_this.workspace.a[k] = t[sortedIndices[k] - 1];
    }
    ny = b_this.workspace.b.size(0);
    t.set_size(b_this.workspace.b.size(0));
    ns = b_this.workspace.b.size(0);
    for (i = 0; i < ns; i++) {
        t[i] = b_this.workspace.b[i];
    }
    for (int k{0}; k < ny; k++) {
        b_this.workspace.b[k] = t[sortedIndices[k] - 1];
    }
    if ((b_this.workspace.b.size(0) == 0) || (b_this.workspace.a.size(0) == 0)) {
        ny = 0;
        y.n = 0;
    } else {
        ns = b_this.workspace.b.size(0);
        ny = b_this.workspace.b[0];
        for (int k{2}; k <= ns; k++) {
            i = b_this.workspace.b[k - 1];
            if (ny < i) {
                ny = i;
            }
        }
        y.n = b_this.workspace.a[b_this.workspace.a.size(0) - 1];
    }
    y.m = ny;
    ns = varargin_2.size(0);
    if (ns < 1) {
        ns = 1;
    }
    y.d.set_size(ns);
    y.colidx.set_size(y.n + 1);
    y.colidx[0] = 1;
    y.rowidx.set_size(ns);
    for (i = 0; i < ns; i++) {
        y.d[i] = 0.0;
        y.rowidx[i] = 0;
    }
    ns = 0;
    i = y.n;
    for (ny = 0; ny < i; ny++) {
        while ((ns + 1 <= nc) && (b_this.workspace.a[ns] == ny + 1)) {
            y.rowidx[ns] = b_this.workspace.b[ns];
            ns++;
        }
        y.colidx[ny + 1] = ns + 1;
    }
    for (int k{0}; k < nc; k++) {
        y.d[k] = varargin_3[sortedIndices[k] - 1];
    }
    y.fillIn();
}

}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for sparse.cpp
///
/// [EOF]
///
