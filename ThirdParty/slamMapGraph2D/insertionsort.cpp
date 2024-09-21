///
/// @file           : insertionsort.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "insertionsort.h"
#include "anonymous_function.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_internal_types.h"
#include "coder_array.h"

/// Function Definitions
///
/// @fn             : insertionsort
/// @brief          :
/// @param          : ::coder::array<int, 1U> &x
///                   int xstart
///                   int xend
///                   const anonymous_function &cmp
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
void insertionsort(::coder::array<int, 1U> &x, int xstart, int xend,
                   const anonymous_function &cmp) {
    int i;
    i = xstart + 1;
    for (int k{i}; k <= xend; k++) {
        int idx;
        int xc;
        bool exitg1;
        xc = x[k - 1] - 1;
        idx = k - 2;
        exitg1 = false;
        while ((!exitg1) && (idx + 1 >= xstart)) {
            int i1;
            bool varargout_1;
            i1 = cmp.workspace.a[x[idx] - 1];
            if (cmp.workspace.a[xc] < i1) {
                varargout_1 = true;
            } else if (cmp.workspace.a[xc] == i1) {
                varargout_1 = (cmp.workspace.b[xc] < cmp.workspace.b[x[idx] - 1]);
            } else {
                varargout_1 = false;
            }
            if (varargout_1) {
                x[idx + 1] = x[idx];
                idx--;
            } else {
                exitg1 = true;
            }
        }
        x[idx + 1] = xc + 1;
    }
}

}  // namespace internal
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for insertionsort.cpp
///
/// [EOF]
///
