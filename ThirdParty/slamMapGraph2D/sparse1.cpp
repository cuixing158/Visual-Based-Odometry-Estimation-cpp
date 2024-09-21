///
/// @file           : sparse1.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "sparse1.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

/// Function Definitions
///
/// @fn             : ctranspose
/// @brief          :
/// @param          : sparse &y
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
void sparse::ctranspose(sparse &y) const {
    ::coder::array<int, 1U> counts;
    int i;
    int loop_ub;
    int nl;
    int numalloc;
    int outridx;
    nl = n;
    y.m = n;
    y.n = m;
    outridx = colidx[colidx.size(0) - 1];
    if (outridx - 1 >= 1) {
        numalloc = outridx - 2;
    } else {
        numalloc = 0;
    }
    y.d.set_size(numalloc + 1);
    i = m + 1;
    y.colidx.set_size(i);
    y.colidx[0] = 1;
    y.rowidx.set_size(numalloc + 1);
    for (int c{0}; c <= numalloc; c++) {
        y.d[c] = 0.0;
        y.rowidx[c] = 0;
    }
    loop_ub = m;
    for (int c{0}; c < loop_ub; c++) {
        y.colidx[c + 1] = 1;
    }
    y.fillIn();
    if ((m != 0) && (n != 0)) {
        numalloc = y.colidx.size(0);
        for (int c{0}; c < numalloc; c++) {
            y.colidx[c] = 0;
        }
        for (numalloc = 0; numalloc <= outridx - 2; numalloc++) {
            y.colidx[rowidx[numalloc]] = y.colidx[rowidx[numalloc]] + 1;
        }
        y.colidx[0] = 1;
        for (numalloc = 2; numalloc <= i; numalloc++) {
            y.colidx[numalloc - 1] = y.colidx[numalloc - 1] + y.colidx[numalloc - 2];
        }
        counts.set_size(m);
        for (i = 0; i < loop_ub; i++) {
            counts[i] = 0;
        }
        for (int c{0}; c < nl; c++) {
            for (numalloc = colidx[c] - 1; numalloc + 1 < colidx[c + 1]; numalloc++) {
                loop_ub = counts[rowidx[numalloc] - 1];
                outridx = (loop_ub + y.colidx[rowidx[numalloc] - 1]) - 1;
                y.d[outridx] = d[numalloc];
                y.rowidx[outridx] = c + 1;
                counts[rowidx[numalloc] - 1] = loop_ub + 1;
            }
        }
    }
}

///
/// @fn             : fillIn
/// @brief          :
/// @param          : void
/// @return         : void
///
void sparse::fillIn() {
    int i;
    int idx;
    idx = 1;
    i = colidx.size(0);
    for (int c{0}; c <= i - 2; c++) {
        int ridx;
        ridx = colidx[c];
        colidx[c] = idx;
        int exitg1;
        int i1;
        do {
            exitg1 = 0;
            i1 = colidx[c + 1];
            if (ridx < i1) {
                double val;
                int currRowIdx;
                val = 0.0;
                currRowIdx = rowidx[ridx - 1];
                while ((ridx < i1) && (rowidx[ridx - 1] == currRowIdx)) {
                    val += d[ridx - 1];
                    ridx++;
                }
                if (val != 0.0) {
                    d[idx - 1] = val;
                    rowidx[idx - 1] = currRowIdx;
                    idx++;
                }
            } else {
                exitg1 = 1;
            }
        } while (exitg1 == 0);
    }
    colidx[colidx.size(0) - 1] = idx;
}

}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for sparse1.cpp
///
/// [EOF]
///
