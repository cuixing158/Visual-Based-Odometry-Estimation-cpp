///
/// @file           : introsort.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

/// @include file    : Include Files
#include "introsort.h"
#include "anonymous_function.h"
#include "heapsort.h"
#include "insertionsort.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_internal_types.h"
#include "stack1.h"
#include "coder_array.h"
#include "coder_bounded_array.h"

/// Function Definitions
///
/// @fn             : introsort
/// @brief          :
/// @param          : ::coder::array<int, 1U> &x
///                   int xend
///                   const anonymous_function &cmp
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
void introsort(::coder::array<int, 1U> &x, int xend,
               const anonymous_function &cmp)
{
  struct_T frame;
  if (xend > 1) {
    if (xend <= 32) {
      insertionsort(x, 1, xend, cmp);
    } else {
      stack st;
      int MAXDEPTH;
      int i;
      int pmax;
      int pmin;
      int pow2p;
      int t;
      bool exitg1;
      pmax = 31;
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        t = (pmin + pmax) >> 1;
        pow2p = 1 << t;
        if (pow2p == xend) {
          pmax = t;
          exitg1 = true;
        } else if (pow2p > xend) {
          pmax = t;
        } else {
          pmin = t;
        }
      }
      MAXDEPTH = (pmax - 1) << 1;
      frame.xstart = 1;
      frame.xend = xend;
      frame.depth = 0;
      pmax = MAXDEPTH << 1;
      st.d.size[0] = pmax;
      for (i = 0; i < pmax; i++) {
        st.d.data[i] = frame;
      }
      st.d.data[0] = frame;
      st.n = 1;
      while (st.n > 0) {
        int frame_tmp_tmp;
        frame_tmp_tmp = st.n - 1;
        frame = st.d.data[st.n - 1];
        st.n--;
        i = frame.xend - frame.xstart;
        if (i + 1 <= 32) {
          insertionsort(x, frame.xstart, frame.xend, cmp);
        } else if (frame.depth == MAXDEPTH) {
          b_heapsort(x, frame.xstart, frame.xend, cmp);
        } else {
          int xmid;
          bool varargout_1;
          xmid = (frame.xstart + i / 2) - 1;
          i = cmp.workspace.a[x[xmid] - 1];
          pmax = x[frame.xstart - 1];
          pmin = cmp.workspace.a[pmax - 1];
          if (i < pmin) {
            varargout_1 = true;
          } else if (i == pmin) {
            varargout_1 =
                (cmp.workspace.b[x[xmid] - 1] < cmp.workspace.b[pmax - 1]);
          } else {
            varargout_1 = false;
          }
          if (varargout_1) {
            x[frame.xstart - 1] = x[xmid];
            x[xmid] = pmax;
          }
          i = x[frame.xend - 1];
          pmax = cmp.workspace.a[i - 1];
          pmin = x[frame.xstart - 1];
          t = cmp.workspace.a[pmin - 1];
          if (pmax < t) {
            varargout_1 = true;
          } else if (pmax == t) {
            varargout_1 = (cmp.workspace.b[i - 1] < cmp.workspace.b[pmin - 1]);
          } else {
            varargout_1 = false;
          }
          if (varargout_1) {
            x[frame.xstart - 1] = i;
            x[frame.xend - 1] = pmin;
          }
          i = x[frame.xend - 1];
          pmax = cmp.workspace.a[i - 1];
          pmin = cmp.workspace.a[x[xmid] - 1];
          if (pmax < pmin) {
            varargout_1 = true;
          } else if (pmax == pmin) {
            varargout_1 =
                (cmp.workspace.b[i - 1] < cmp.workspace.b[x[xmid] - 1]);
          } else {
            varargout_1 = false;
          }
          if (varargout_1) {
            t = x[xmid];
            x[xmid] = i;
            x[frame.xend - 1] = t;
          }
          pow2p = x[xmid] - 1;
          x[xmid] = x[frame.xend - 2];
          x[frame.xend - 2] = pow2p + 1;
          pmax = frame.xstart - 1;
          pmin = frame.xend - 2;
          int exitg2;
          do {
            int exitg3;
            exitg2 = 0;
            pmax++;
            do {
              exitg3 = 0;
              i = cmp.workspace.a[x[pmax] - 1];
              if (i < cmp.workspace.a[pow2p]) {
                varargout_1 = true;
              } else if (i == cmp.workspace.a[pow2p]) {
                varargout_1 =
                    (cmp.workspace.b[x[pmax] - 1] < cmp.workspace.b[pow2p]);
              } else {
                varargout_1 = false;
              }
              if (varargout_1) {
                pmax++;
              } else {
                exitg3 = 1;
              }
            } while (exitg3 == 0);
            pmin--;
            do {
              exitg3 = 0;
              i = cmp.workspace.a[x[pmin] - 1];
              if (cmp.workspace.a[pow2p] < i) {
                varargout_1 = true;
              } else if (cmp.workspace.a[pow2p] == i) {
                varargout_1 =
                    (cmp.workspace.b[pow2p] < cmp.workspace.b[x[pmin] - 1]);
              } else {
                varargout_1 = false;
              }
              if (varargout_1) {
                pmin--;
              } else {
                exitg3 = 1;
              }
            } while (exitg3 == 0);
            if (pmax + 1 >= pmin + 1) {
              exitg2 = 1;
            } else {
              t = x[pmax];
              x[pmax] = x[pmin];
              x[pmin] = t;
            }
          } while (exitg2 == 0);
          x[frame.xend - 2] = x[pmax];
          x[pmax] = pow2p + 1;
          if (pmax + 2 < frame.xend) {
            st.d.data[frame_tmp_tmp].xstart = pmax + 2;
            st.d.data[frame_tmp_tmp].xend = frame.xend;
            st.d.data[frame_tmp_tmp].depth = frame.depth + 1;
            st.n = frame_tmp_tmp + 1;
          }
          if (frame.xstart < pmax + 1) {
            st.d.data[st.n].xstart = frame.xstart;
            st.d.data[st.n].xend = pmax + 1;
            st.d.data[st.n].depth = frame.depth + 1;
            st.n++;
          }
        }
      }
    }
  }
}

} // namespace internal
} // namespace coder
} // namespace SlamGraph2D

///
/// File trailer for introsort.cpp
///
/// [EOF]
///
