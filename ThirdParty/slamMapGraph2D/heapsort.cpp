///
/// @file           : heapsort.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "heapsort.h"
#include "anonymous_function.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_internal_types.h"
#include "coder_array.h"

/// Function Declarations
namespace SlamGraph2D {
namespace coder {
namespace internal {
static void heapify(::coder::array<int, 1U> &x, int idx, int xstart, int xend,
                    const anonymous_function &cmp);

}
}  // namespace coder
}  // namespace SlamGraph2D

/// Function Definitions
///
/// @fn             : heapify
/// @brief          :
/// @param          : ::coder::array<int, 1U> &x
///                   int idx
///                   int xstart
///                   int xend
///                   const anonymous_function &cmp
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
static void heapify(::coder::array<int, 1U> &x, int idx, int xstart, int xend,
                    const anonymous_function &cmp) {
    int extremum;
    int extremumIdx;
    int i;
    int i1;
    int leftIdx;
    bool changed;
    bool exitg1;
    bool varargout_1;
    changed = true;
    extremumIdx = (idx + xstart) - 2;
    leftIdx = ((idx << 1) + xstart) - 2;
    exitg1 = false;
    while ((!exitg1) && (leftIdx + 1 < xend)) {
        int cmpIdx;
        int i2;
        int xcmp;
        changed = false;
        extremum = x[extremumIdx];
        cmpIdx = leftIdx;
        xcmp = x[leftIdx];
        i = cmp.workspace.a[x[leftIdx] - 1];
        i1 = x[leftIdx + 1] - 1;
        i2 = cmp.workspace.a[i1];
        if (i < i2) {
            varargout_1 = true;
        } else if (i == i2) {
            varargout_1 = (cmp.workspace.b[x[leftIdx] - 1] < cmp.workspace.b[i1]);
        } else {
            varargout_1 = false;
        }
        if (varargout_1) {
            cmpIdx = leftIdx + 1;
            xcmp = x[leftIdx + 1];
        }
        i = cmp.workspace.a[x[extremumIdx] - 1];
        i1 = cmp.workspace.a[xcmp - 1];
        if (i < i1) {
            varargout_1 = true;
        } else if (i == i1) {
            varargout_1 =
                (cmp.workspace.b[x[extremumIdx] - 1] < cmp.workspace.b[xcmp - 1]);
        } else {
            varargout_1 = false;
        }
        if (varargout_1) {
            x[extremumIdx] = xcmp;
            x[cmpIdx] = extremum;
            extremumIdx = cmpIdx;
            leftIdx = ((((cmpIdx - xstart) + 2) << 1) + xstart) - 2;
            changed = true;
        } else {
            exitg1 = true;
        }
    }
    if (changed && (leftIdx + 1 <= xend)) {
        extremum = x[extremumIdx];
        i = cmp.workspace.a[x[extremumIdx] - 1];
        i1 = cmp.workspace.a[x[leftIdx] - 1];
        if (i < i1) {
            varargout_1 = true;
        } else if (i == i1) {
            varargout_1 = (cmp.workspace.b[x[extremumIdx] - 1] <
                           cmp.workspace.b[x[leftIdx] - 1]);
        } else {
            varargout_1 = false;
        }
        if (varargout_1) {
            x[extremumIdx] = x[leftIdx];
            x[leftIdx] = extremum;
        }
    }
}

///
/// @fn             : b_heapsort
/// @brief          :
/// @param          : ::coder::array<int, 1U> &x
///                   int xstart
///                   int xend
///                   const anonymous_function &cmp
/// @return         : void
///
void b_heapsort(::coder::array<int, 1U> &x, int xstart, int xend,
                const anonymous_function &cmp) {
    int idx;
    int n;
    n = (xend - xstart) - 1;
    for (idx = n + 2; idx >= 1; idx--) {
        heapify(x, idx, xstart, xend, cmp);
    }
    for (int k{0}; k <= n; k++) {
        int t;
        idx = (xend - k) - 1;
        t = x[idx];
        x[idx] = x[xstart - 1];
        x[xstart - 1] = t;
        heapify(x, 1, xstart, idx, cmp);
    }
}

}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for heapsort.cpp
///
/// [EOF]
///
