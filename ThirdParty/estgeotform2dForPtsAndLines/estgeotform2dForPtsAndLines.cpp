#include "estgeotform2dForPtsAndLines.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "rt_defines.h"
#include <cfloat>
#include <cmath>
#include <cstring>

namespace estgeotform2dForPtsAndLines {
static unsigned int state[625];

static bool isInitialized_estgeotform2dForPtsAndLines{false};

} // namespace estgeotform2dForPtsAndLines

namespace estgeotform2dForPtsAndLines {
static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 1U> &in3,
                             const ::coder::array<double, 2U> &in4);

static void binary_expand_op(::coder::array<double, 2U> &in1,
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 2U> &in3);

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 1U> &in2,
                             const ::coder::array<double, 2U> &in3,
                             const ::coder::array<double, 1U> &in4,
                             const ::coder::array<double, 1U> &in5,
                             const ::coder::array<double, 1U> &in6);

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 1U> &in3,
                             const ::coder::array<double, 1U> &in4,
                             const ::coder::array<double, 1U> &in5);

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 2U> &in3,
                             const ::coder::array<double, 1U> &in4);

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 1U> &in2,
                             const ::coder::array<double, 1U> &in3);

namespace coder {
static bool all(const bool x[2]);

static bool any(const ::coder::array<bool, 1U> &x);

static bool b_all(const bool x[6]);

static void b_cosd(double &x);

static double b_rand();

static void b_sind(double &x);

static double combineVectorElements(const ::coder::array<double, 1U> &x);

namespace internal {
namespace blas {
static void mtimes(const double A[6], const ::coder::array<double, 2U> &B,
                   ::coder::array<double, 2U> &C);

static double xnrm2(int n, const ::coder::array<double, 2U> &x, int ix0);

static double xnrm2(const double x[4]);

static double xrotg(double &a, double &b, double &s);

static void xswap(double x[4]);

} // namespace blas
static void sort(double x[2]);

static void svd(const double A[4], double U[4], double s[2], double V[4]);

} // namespace internal
static void mean(const ::coder::array<double, 2U> &x, double y[2]);

static void mldivide(const ::coder::array<double, 2U> &A,
                     const ::coder::array<double, 1U> &B, double Y[2]);

static void randperm(double n, double p[2]);

} // namespace coder
static void eml_rand_mt19937ar_stateful_init();

static double evaluateModel(const double modelIn[6],
                            const ::coder::array<double, 2U> &matchedPoints1,
                            const ::coder::array<double, 2U> &matchedPoints2,
                            const ::coder::array<double, 2U> &matchedLines1,
                            const ::coder::array<double, 2U> &matchedLines2,
                            ::coder::array<double, 1U> &distances);

static bool judgeLinesValid(const ::coder::array<double, 2U> &lines);

static double rt_atan2d_snf(double u0, double u1);

static double rt_hypotd_snf(double u0, double u1);

static double rt_powd_snf(double u0, double u1);

static double rt_remd_snf(double u0, double u1);

static void
svdFitForInliersPointsAndLines(const ::coder::array<double, 2U> &inliersPts1,
                               const ::coder::array<double, 2U> &inliersPts2,
                               const ::coder::array<double, 2U> &inliersL1,
                               const ::coder::array<double, 2U> &inliersL2,
                               double tform[6]);

} // namespace estgeotform2dForPtsAndLines

namespace estgeotform2dForPtsAndLines {
static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 1U> &in3,
                             const ::coder::array<double, 2U> &in4)
{
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in4.size(0) == 1) {
    loop_ub = in3.size(0);
  } else {
    loop_ub = in4.size(0);
  }
  in1.set_size(loop_ub);
  stride_0_0 = (in3.size(0) != 1);
  stride_1_0 = (in4.size(0) != 1);
  for (int i{0}; i < loop_ub; i++) {
    double b_varargin_1;
    double varargin_1;
    int i1;
    i1 = i * stride_1_0;
    varargin_1 = in3[i * stride_0_0];
    b_varargin_1 = in4[i1 + in4.size(0) * 2] - in4[i1];
    in1[i] = varargin_1 * varargin_1 + b_varargin_1 * b_varargin_1;
  }
}

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 1U> &in2,
                             const ::coder::array<double, 2U> &in3,
                             const ::coder::array<double, 1U> &in4,
                             const ::coder::array<double, 1U> &in5,
                             const ::coder::array<double, 1U> &in6)
{
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0_tmp;
  int stride_2_0;
  int stride_4_0;
  int stride_5_0;
  if (in3.size(0) == 1) {
    i = in4.size(0);
  } else {
    i = in3.size(0);
  }
  if (in6.size(0) == 1) {
    if (in5.size(0) == 1) {
      if (i == 1) {
        if (in3.size(0) == 1) {
          loop_ub = in2.size(0);
        } else {
          loop_ub = in3.size(0);
        }
      } else {
        loop_ub = i;
      }
    } else {
      loop_ub = in5.size(0);
    }
  } else {
    loop_ub = in6.size(0);
  }
  in1.set_size(loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0_tmp = (in3.size(0) != 1);
  stride_2_0 = (in4.size(0) != 1);
  stride_4_0 = (in5.size(0) != 1);
  stride_5_0 = (in6.size(0) != 1);
  for (i = 0; i < loop_ub; i++) {
    int i1;
    i1 = i * stride_1_0_tmp;
    in1[i] = ((in2[i * stride_0_0] * in3[i1] +
               in4[i * stride_2_0] * in3[i1 + in3.size(0)]) +
              in5[i * stride_4_0]) -
             in6[i * stride_5_0];
  }
}

static void binary_expand_op(::coder::array<double, 2U> &in1,
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 2U> &in3)
{
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in3.size(0) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in3.size(0);
  }
  in1.set_size(loop_ub, 2);
  stride_0_0 = (in2.size(1) != 1);
  stride_1_0 = (in3.size(0) != 1);
  for (int i{0}; i < 2; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      in1[i1 + in1.size(0) * i] = in2[i + 2 * (i1 * stride_0_0)] -
                                  in3[i1 * stride_1_0 + in3.size(0) * i];
    }
  }
}

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 1U> &in3,
                             const ::coder::array<double, 1U> &in4,
                             const ::coder::array<double, 1U> &in5)
{
  ::coder::array<double, 1U> b_in1;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0_tmp;
  int stride_2_0;
  int stride_4_0;
  int stride_5_0;
  if (in2.size(0) == 1) {
    i = in3.size(0);
  } else {
    i = in2.size(0);
  }
  if (in5.size(0) == 1) {
    if (in4.size(0) == 1) {
      if (i == 1) {
        if (in2.size(0) == 1) {
          loop_ub = in1.size(0);
        } else {
          loop_ub = in2.size(0);
        }
      } else {
        loop_ub = i;
      }
    } else {
      loop_ub = in4.size(0);
    }
  } else {
    loop_ub = in5.size(0);
  }
  b_in1.set_size(loop_ub);
  stride_0_0 = (in1.size(0) != 1);
  stride_1_0_tmp = (in2.size(0) != 1);
  stride_2_0 = (in3.size(0) != 1);
  stride_4_0 = (in4.size(0) != 1);
  stride_5_0 = (in5.size(0) != 1);
  for (i = 0; i < loop_ub; i++) {
    int in1_tmp;
    in1_tmp = i * stride_1_0_tmp;
    b_in1[i] = ((in1[i * stride_0_0] * in2[in1_tmp + in2.size(0) * 2] +
                 in3[i * stride_2_0] * in2[in1_tmp + in2.size(0) * 3]) +
                in4[i * stride_4_0]) -
               in5[i * stride_5_0];
  }
  in1.set_size(b_in1.size(0));
  loop_ub = b_in1.size(0);
  for (i = 0; i < loop_ub; i++) {
    in1[i] = b_in1[i];
  }
}

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 2U> &in3,
                             const ::coder::array<double, 1U> &in4)
{
  ::coder::array<double, 1U> b_in2;
  int i;
  int loop_ub;
  int stride_0_0_tmp;
  int stride_1_0_tmp;
  int stride_4_0;
  int stride_5_0;
  if (in4.size(0) == 1) {
    i = in1.size(0);
  } else {
    i = in4.size(0);
  }
  if (i == 1) {
    if (in3.size(0) == 1) {
      loop_ub = in2.size(0);
    } else {
      loop_ub = in3.size(0);
    }
  } else {
    loop_ub = i;
  }
  b_in2.set_size(loop_ub);
  stride_0_0_tmp = (in2.size(0) != 1);
  stride_1_0_tmp = (in3.size(0) != 1);
  stride_4_0 = (in1.size(0) != 1);
  stride_5_0 = (in4.size(0) != 1);
  for (i = 0; i < loop_ub; i++) {
    int b_in2_tmp;
    int in2_tmp;
    in2_tmp = i * stride_0_0_tmp;
    b_in2_tmp = i * stride_1_0_tmp;
    b_in2[i] = (in2[in2_tmp] * in3[b_in2_tmp] +
                in2[in2_tmp + in2.size(0)] * in3[b_in2_tmp + in3.size(0)]) /
               (in1[i * stride_4_0] * in4[i * stride_5_0]);
  }
  in1.set_size(b_in2.size(0));
  loop_ub = b_in2.size(0);
  for (i = 0; i < loop_ub; i++) {
    in1[i] = b_in2[i];
  }
}

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 1U> &in2,
                             const ::coder::array<double, 1U> &in3)
{
  ::coder::array<double, 1U> b_in1;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0_tmp;
  int stride_2_0;
  if (in2.size(0) == 1) {
    i = in3.size(0);
  } else {
    i = in2.size(0);
  }
  if (i == 1) {
    if (in2.size(0) == 1) {
      loop_ub = in1.size(0);
    } else {
      loop_ub = in2.size(0);
    }
  } else {
    loop_ub = i;
  }
  b_in1.set_size(loop_ub);
  stride_0_0 = (in1.size(0) != 1);
  stride_1_0_tmp = (in2.size(0) != 1);
  stride_2_0 = (in3.size(0) != 1);
  for (i = 0; i < loop_ub; i++) {
    double in1_tmp;
    in1_tmp = in2[i * stride_1_0_tmp];
    b_in1[i] = (in1[i * stride_0_0] / in1_tmp + in3[i * stride_2_0] / in1_tmp) /
               2.0 * 0.5;
  }
  in1.set_size(b_in1.size(0));
  loop_ub = b_in1.size(0);
  for (i = 0; i < loop_ub; i++) {
    in1[i] = b_in1[i];
  }
}

namespace coder {
static bool all(const bool x[2])
{
  int k;
  bool exitg1;
  bool y;
  y = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 2)) {
    if (!x[k]) {
      y = false;
      exitg1 = true;
    } else {
      k++;
    }
  }
  return y;
}

static bool any(const ::coder::array<bool, 1U> &x)
{
  int ix;
  bool exitg1;
  bool y;
  y = false;
  ix = 1;
  exitg1 = false;
  while ((!exitg1) && (ix <= x.size(0))) {
    if (x[ix - 1]) {
      y = true;
      exitg1 = true;
    } else {
      ix++;
    }
  }
  return y;
}

static bool b_all(const bool x[6])
{
  int k;
  bool exitg1;
  bool y;
  y = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 6)) {
    if (!x[k]) {
      y = false;
      exitg1 = true;
    } else {
      k++;
    }
  }
  return y;
}

static void b_cosd(double &x)
{
  if (std::isinf(x) || std::isnan(x)) {
    x = rtNaN;
  } else {
    double absx;
    signed char n;
    x = rt_remd_snf(x, 360.0);
    absx = std::abs(x);
    if (absx > 180.0) {
      if (x > 0.0) {
        x -= 360.0;
      } else {
        x += 360.0;
      }
      absx = std::abs(x);
    }
    if (absx <= 45.0) {
      x *= 0.017453292519943295;
      n = 0;
    } else if (absx <= 135.0) {
      if (x > 0.0) {
        x = 0.017453292519943295 * (x - 90.0);
        n = 1;
      } else {
        x = 0.017453292519943295 * (x + 90.0);
        n = -1;
      }
    } else if (x > 0.0) {
      x = 0.017453292519943295 * (x - 180.0);
      n = 2;
    } else {
      x = 0.017453292519943295 * (x + 180.0);
      n = -2;
    }
    if (n == 0) {
      x = std::cos(x);
    } else if (n == 1) {
      x = -std::sin(x);
    } else if (n == -1) {
      x = std::sin(x);
    } else {
      x = -std::cos(x);
    }
  }
}

static double b_rand()
{
  double r;
  /* ========================= COPYRIGHT NOTICE ============================ */
  /*  This is a uniform (0,1) pseudorandom number generator based on:        */
  /*                                                                         */
  /*  A C-program for MT19937, with initialization improved 2002/1/26.       */
  /*  Coded by Takuji Nishimura and Makoto Matsumoto.                        */
  /*                                                                         */
  /*  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,      */
  /*  All rights reserved.                                                   */
  /*                                                                         */
  /*  Redistribution and use in source and binary forms, with or without     */
  /*  modification, are permitted provided that the following conditions     */
  /*  are met:                                                               */
  /*                                                                         */
  /*    1. Redistributions of source code must retain the above copyright    */
  /*       notice, this list of conditions and the following disclaimer.     */
  /*                                                                         */
  /*    2. Redistributions in binary form must reproduce the above copyright */
  /*       notice, this list of conditions and the following disclaimer      */
  /*       in the documentation and/or other materials provided with the     */
  /*       distribution.                                                     */
  /*                                                                         */
  /*    3. The names of its contributors may not be used to endorse or       */
  /*       promote products derived from this software without specific      */
  /*       prior written permission.                                         */
  /*                                                                         */
  /*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS    */
  /*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT      */
  /*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  */
  /*  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT  */
  /*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,  */
  /*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT       */
  /*  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,  */
  /*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY  */
  /*  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT    */
  /*  (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE */
  /*  OF THIS  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */
  /*                                                                         */
  /* =============================   END   ================================= */
  unsigned int u[2];
  do {
    for (int j{0}; j < 2; j++) {
      unsigned int mti;
      unsigned int y;
      mti = state[624] + 1U;
      if (state[624] + 1U >= 625U) {
        for (int kk{0}; kk < 227; kk++) {
          y = (state[kk] & 2147483648U) | (state[kk + 1] & 2147483647U);
          if ((y & 1U) == 0U) {
            y >>= 1U;
          } else {
            y = y >> 1U ^ 2567483615U;
          }
          state[kk] = state[kk + 397] ^ y;
        }
        for (int kk{0}; kk < 396; kk++) {
          y = (state[kk + 227] & 2147483648U) | (state[kk + 228] & 2147483647U);
          if ((y & 1U) == 0U) {
            y >>= 1U;
          } else {
            y = y >> 1U ^ 2567483615U;
          }
          state[kk + 227] = state[kk] ^ y;
        }
        y = (state[623] & 2147483648U) | (state[0] & 2147483647U);
        if ((y & 1U) == 0U) {
          y >>= 1U;
        } else {
          y = y >> 1U ^ 2567483615U;
        }
        state[623] = state[396] ^ y;
        mti = 1U;
      }
      y = state[static_cast<int>(mti) - 1];
      state[624] = mti;
      y ^= y >> 11U;
      y ^= y << 7U & 2636928640U;
      y ^= y << 15U & 4022730752U;
      u[j] = y ^ y >> 18U;
    }
    u[0] >>= 5U;
    u[1] >>= 6U;
    r = 1.1102230246251565E-16 *
        (static_cast<double>(u[0]) * 6.7108864E+7 + static_cast<double>(u[1]));
  } while (r == 0.0);
  return r;
}

static void b_sind(double &x)
{
  if (std::isinf(x) || std::isnan(x)) {
    x = rtNaN;
  } else {
    double absx;
    signed char n;
    x = rt_remd_snf(x, 360.0);
    absx = std::abs(x);
    if (absx > 180.0) {
      if (x > 0.0) {
        x -= 360.0;
      } else {
        x += 360.0;
      }
      absx = std::abs(x);
    }
    if (absx <= 45.0) {
      x *= 0.017453292519943295;
      n = 0;
    } else if (absx <= 135.0) {
      if (x > 0.0) {
        x = 0.017453292519943295 * (x - 90.0);
        n = 1;
      } else {
        x = 0.017453292519943295 * (x + 90.0);
        n = -1;
      }
    } else if (x > 0.0) {
      x = 0.017453292519943295 * (x - 180.0);
      n = 2;
    } else {
      x = 0.017453292519943295 * (x + 180.0);
      n = -2;
    }
    if (n == 0) {
      x = std::sin(x);
    } else if (n == 1) {
      x = std::cos(x);
    } else if (n == -1) {
      x = -std::cos(x);
    } else {
      x = -std::sin(x);
    }
  }
}

static double combineVectorElements(const ::coder::array<double, 1U> &x)
{
  double y;
  if (x.size(0) == 0) {
    y = 0.0;
  } else {
    int firstBlockLength;
    int lastBlockLength;
    int nblocks;
    if (x.size(0) <= 1024) {
      firstBlockLength = x.size(0);
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = static_cast<int>(static_cast<unsigned int>(x.size(0)) >> 10);
      lastBlockLength = x.size(0) - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    y = x[0];
    for (int k{2}; k <= firstBlockLength; k++) {
      y += x[k - 1];
    }
    for (int ib{2}; ib <= nblocks; ib++) {
      double bsum;
      int hi;
      firstBlockLength = (ib - 1) << 10;
      bsum = x[firstBlockLength];
      if (ib == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (int k{2}; k <= hi; k++) {
        bsum += x[(firstBlockLength + k) - 1];
      }
      y += bsum;
    }
  }
  return y;
}

namespace internal {
namespace blas {
static void mtimes(const double A[6], const ::coder::array<double, 2U> &B,
                   ::coder::array<double, 2U> &C)
{
  int n;
  n = B.size(1);
  C.set_size(2, B.size(1));
  for (int j{0}; j < n; j++) {
    int boffset;
    int coffset;
    coffset = j << 1;
    boffset = j * B.size(0);
    for (int i{0}; i < 2; i++) {
      C[coffset + i] = (A[i] * B[boffset] + A[i + 2] * B[boffset + 1]) +
                       A[i + 4] * B[boffset + 2];
    }
  }
}

static double xnrm2(int n, const ::coder::array<double, 2U> &x, int ix0)
{
  double y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (int k{ix0}; k <= kend; k++) {
        double absxk;
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * std::sqrt(y);
    }
  }
  return y;
}

static double xnrm2(const double x[4])
{
  double absxk;
  double scale;
  double t;
  double y;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }
  absxk = std::abs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  return scale * std::sqrt(y);
}

static double xrotg(double &a, double &b, double &s)
{
  double absa;
  double absb;
  double c;
  double roe;
  double scale;
  roe = b;
  absa = std::abs(a);
  absb = std::abs(b);
  if (absa > absb) {
    roe = a;
  }
  scale = absa + absb;
  if (scale == 0.0) {
    s = 0.0;
    c = 1.0;
    a = 0.0;
    b = 0.0;
  } else {
    double ads;
    double bds;
    ads = absa / scale;
    bds = absb / scale;
    scale *= std::sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      scale = -scale;
    }
    c = a / scale;
    s = b / scale;
    if (absa > absb) {
      b = s;
    } else if (c != 0.0) {
      b = 1.0 / c;
    } else {
      b = 1.0;
    }
    a = scale;
  }
  return c;
}

static void xswap(double x[4])
{
  double temp;
  temp = x[0];
  x[0] = x[2];
  x[2] = temp;
  temp = x[1];
  x[1] = x[3];
  x[3] = temp;
}

} // namespace blas
static void sort(double x[2])
{
  if ((!(x[0] <= x[1])) && (!std::isnan(x[1]))) {
    double tmp;
    tmp = x[0];
    x[0] = x[1];
    x[1] = tmp;
  }
}

static void svd(const double A[4], double U[4], double s[2], double V[4])
{
  double b_s[2];
  double e[2];
  double A_idx_3;
  double nrm;
  double rt;
  double sm;
  double snorm;
  double sqds;
  double temp;
  int iter;
  int kase;
  int m;
  int q;
  int qs;
  temp = A[0];
  sm = A[1];
  sqds = A[2];
  A_idx_3 = A[3];
  nrm = blas::xnrm2(A);
  if (nrm > 0.0) {
    if (A[0] < 0.0) {
      b_s[0] = -nrm;
    } else {
      b_s[0] = nrm;
    }
    if (std::abs(b_s[0]) >= 1.0020841800044864E-292) {
      rt = 1.0 / b_s[0];
      temp = rt * A[0];
      sm = rt * A[1];
    } else {
      temp = A[0] / b_s[0];
      sm = A[1] / b_s[0];
    }
    temp++;
    b_s[0] = -b_s[0];
    rt = -((temp * A[2] + sm * A[3]) / temp);
    if (!(rt == 0.0)) {
      sqds = A[2] + rt * temp;
      A_idx_3 = A[3] + rt * sm;
    }
  } else {
    b_s[0] = 0.0;
  }
  m = 2;
  b_s[1] = A_idx_3;
  e[0] = sqds;
  e[1] = 0.0;
  U[2] = 0.0;
  U[3] = 1.0;
  if (b_s[0] != 0.0) {
    rt = -((temp * 0.0 + sm) / temp);
    if (!(rt == 0.0)) {
      U[2] = rt * temp;
      U[3] = rt * sm + 1.0;
    }
    U[1] = -sm;
    U[0] = -temp + 1.0;
  } else {
    U[1] = 0.0;
    U[0] = 1.0;
  }
  V[2] = 0.0;
  V[3] = 1.0;
  V[1] = 0.0;
  V[0] = 1.0;
  for (q = 0; q < 2; q++) {
    nrm = b_s[q];
    if (nrm != 0.0) {
      rt = std::abs(nrm);
      nrm /= rt;
      b_s[q] = rt;
      if (q + 1 < 2) {
        e[0] /= nrm;
      }
      kase = q << 1;
      qs = kase + 2;
      for (int k{kase + 1}; k <= qs; k++) {
        U[k - 1] *= nrm;
      }
    }
    if ((q + 1 < 2) && (e[0] != 0.0)) {
      rt = std::abs(e[0]);
      nrm = rt / e[0];
      e[0] = rt;
      b_s[1] *= nrm;
      V[2] *= nrm;
      V[3] *= nrm;
    }
  }
  iter = 0;
  snorm = std::fmax(std::fmax(0.0, std::fmax(std::abs(b_s[0]), std::abs(e[0]))),
                    std::fmax(std::abs(b_s[1]), 0.0));
  while ((m > 0) && (iter < 75)) {
    int ii_tmp_tmp;
    bool exitg1;
    ii_tmp_tmp = m - 1;
    q = m - 1;
    exitg1 = false;
    while (!(exitg1 || (q == 0))) {
      nrm = std::abs(e[0]);
      if ((nrm <=
           2.2204460492503131E-16 * (std::abs(b_s[0]) + std::abs(b_s[1]))) ||
          (nrm <= 1.0020841800044864E-292) ||
          ((iter > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
        e[0] = 0.0;
        exitg1 = true;
      } else {
        q = 0;
      }
    }
    if (q == m - 1) {
      kase = 4;
    } else {
      qs = m;
      kase = m;
      exitg1 = false;
      while ((!exitg1) && (kase >= q)) {
        qs = kase;
        if (kase == q) {
          exitg1 = true;
        } else {
          nrm = 0.0;
          if (kase < m) {
            nrm = std::abs(e[0]);
          }
          if (kase > q + 1) {
            nrm += std::abs(e[0]);
          }
          rt = std::abs(b_s[kase - 1]);
          if ((rt <= 2.2204460492503131E-16 * nrm) ||
              (rt <= 1.0020841800044864E-292)) {
            b_s[kase - 1] = 0.0;
            exitg1 = true;
          } else {
            kase--;
          }
        }
      }
      if (qs == q) {
        kase = 3;
      } else if (qs == m) {
        kase = 1;
      } else {
        kase = 2;
        q = qs;
      }
    }
    switch (kase) {
    case 1: {
      double f;
      f = e[0];
      e[0] = 0.0;
      for (int k{ii_tmp_tmp}; k >= q + 1; k--) {
        A_idx_3 = blas::xrotg(b_s[0], f, rt);
        kase = (m - 1) << 1;
        temp = A_idx_3 * V[0] + rt * V[kase];
        V[kase] = A_idx_3 * V[kase] - rt * V[0];
        V[0] = temp;
        nrm = V[kase + 1];
        temp = A_idx_3 * V[1] + rt * nrm;
        V[kase + 1] = A_idx_3 * nrm - rt * V[1];
        V[1] = temp;
      }
    } break;
    case 2: {
      double f;
      f = e[q - 1];
      e[q - 1] = 0.0;
      for (int k{q + 1}; k <= m; k++) {
        A_idx_3 = blas::xrotg(b_s[k - 1], f, sm);
        rt = e[k - 1];
        f = -sm * rt;
        e[k - 1] = rt * A_idx_3;
        qs = (k - 1) << 1;
        kase = (q - 1) << 1;
        temp = A_idx_3 * U[qs] + sm * U[kase];
        U[kase] = A_idx_3 * U[kase] - sm * U[qs];
        U[qs] = temp;
        nrm = U[kase + 1];
        rt = U[qs + 1];
        U[kase + 1] = A_idx_3 * nrm - sm * rt;
        U[qs + 1] = A_idx_3 * rt + sm * nrm;
      }
    } break;
    case 3: {
      double f;
      double scale;
      nrm = b_s[m - 1];
      scale = std::fmax(
          std::fmax(std::fmax(std::fmax(std::abs(nrm), std::abs(b_s[0])),
                              std::abs(e[0])),
                    std::abs(b_s[q])),
          std::abs(e[q]));
      sm = nrm / scale;
      rt = b_s[0] / scale;
      nrm = e[0] / scale;
      sqds = b_s[q] / scale;
      temp = ((rt + sm) * (rt - sm) + nrm * nrm) / 2.0;
      A_idx_3 = sm * nrm;
      A_idx_3 *= A_idx_3;
      if ((temp != 0.0) || (A_idx_3 != 0.0)) {
        rt = std::sqrt(temp * temp + A_idx_3);
        if (temp < 0.0) {
          rt = -rt;
        }
        rt = A_idx_3 / (temp + rt);
      } else {
        rt = 0.0;
      }
      f = (sqds + sm) * (sqds - sm) + rt;
      rt = sqds * (e[q] / scale);
      for (int k{q + 1}; k < 2; k++) {
        A_idx_3 = blas::xrotg(f, rt, sm);
        f = A_idx_3 * b_s[0] + sm * e[0];
        nrm = A_idx_3 * e[0] - sm * b_s[0];
        e[0] = nrm;
        rt = sm * b_s[1];
        b_s[1] *= A_idx_3;
        temp = A_idx_3 * V[0] + sm * V[2];
        V[2] = A_idx_3 * V[2] - sm * V[0];
        V[0] = temp;
        temp = A_idx_3 * V[1] + sm * V[3];
        V[3] = A_idx_3 * V[3] - sm * V[1];
        V[1] = temp;
        b_s[0] = f;
        A_idx_3 = blas::xrotg(b_s[0], rt, sm);
        f = A_idx_3 * nrm + sm * b_s[1];
        b_s[1] = -sm * nrm + A_idx_3 * b_s[1];
        rt = sm * e[1];
        e[1] *= A_idx_3;
        temp = A_idx_3 * U[0] + sm * U[2];
        U[2] = A_idx_3 * U[2] - sm * U[0];
        U[0] = temp;
        temp = A_idx_3 * U[1] + sm * U[3];
        U[3] = A_idx_3 * U[3] - sm * U[1];
        U[1] = temp;
      }
      e[0] = f;
      iter++;
    } break;
    default:
      if (b_s[q] < 0.0) {
        b_s[q] = -b_s[q];
        kase = q << 1;
        qs = kase + 2;
        for (int k{kase + 1}; k <= qs; k++) {
          V[k - 1] = -V[k - 1];
        }
      }
      while ((q + 1 < 2) && (b_s[0] < b_s[1])) {
        rt = b_s[0];
        b_s[0] = b_s[1];
        b_s[1] = rt;
        blas::xswap(V);
        blas::xswap(U);
        q = 1;
      }
      iter = 0;
      m--;
      break;
    }
  }
  s[0] = b_s[0];
  s[1] = b_s[1];
}

} // namespace internal
static void mean(const ::coder::array<double, 2U> &x, double y[2])
{
  if (x.size(0) == 0) {
    y[0] = 0.0;
    y[1] = 0.0;
  } else {
    int firstBlockLength;
    int lastBlockLength;
    int nblocks;
    if (x.size(0) <= 1024) {
      firstBlockLength = x.size(0);
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = static_cast<int>(static_cast<unsigned int>(x.size(0)) >> 10);
      lastBlockLength = x.size(0) - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    for (int xi{0}; xi < 2; xi++) {
      int xpageoffset;
      xpageoffset = xi * x.size(0);
      y[xi] = x[xpageoffset];
      for (int k{2}; k <= firstBlockLength; k++) {
        y[xi] += x[(xpageoffset + k) - 1];
      }
      for (int ib{2}; ib <= nblocks; ib++) {
        double bsum;
        int hi;
        int xblockoffset;
        xblockoffset = xpageoffset + ((ib - 1) << 10);
        bsum = x[xblockoffset];
        if (ib == nblocks) {
          hi = lastBlockLength;
        } else {
          hi = 1024;
        }
        for (int k{2}; k <= hi; k++) {
          bsum += x[(xblockoffset + k) - 1];
        }
        y[xi] += bsum;
      }
    }
  }
  y[0] /= static_cast<double>(x.size(0));
  y[1] /= static_cast<double>(x.size(0));
}

static void mldivide(const ::coder::array<double, 2U> &A,
                     const ::coder::array<double, 1U> &B, double Y[2])
{
  ::coder::array<double, 2U> b_A;
  ::coder::array<double, 1U> b_B;
  double tau_data[2];
  if ((A.size(0) == 0) || (B.size(0) == 0)) {
    Y[0] = 0.0;
    Y[1] = 0.0;
  } else if (A.size(0) == 2) {
    double temp;
    int ix;
    int knt;
    if (std::abs(A[1]) > std::abs(A[0])) {
      knt = 1;
      ix = 0;
    } else {
      knt = 0;
      ix = 1;
    }
    temp = A[ix] / A[knt];
    Y[1] = (B[ix] - B[knt] * temp) /
           (A[ix + A.size(0)] - temp * A[knt + A.size(0)]);
    Y[0] = (B[knt] - Y[1] * A[knt + A.size(0)]) / A[knt];
  } else {
    double vn1[2];
    double vn2[2];
    double work[2];
    double temp;
    int i;
    int ix;
    int knt;
    int m;
    int ma;
    int rankA;
    int u0;
    signed char jpvt[2];
    b_A.set_size(A.size(0), 2);
    knt = A.size(0) << 1;
    for (i = 0; i < knt; i++) {
      b_A[i] = A[i];
    }
    m = A.size(0);
    u0 = A.size(0);
    if (u0 > 2) {
      u0 = 2;
    }
    std::memset(&tau_data[0], 0,
                static_cast<unsigned int>(u0) * sizeof(double));
    ma = A.size(0);
    jpvt[0] = 1;
    work[0] = 0.0;
    temp = internal::blas::xnrm2(A.size(0), A, 1);
    vn1[0] = temp;
    vn2[0] = temp;
    jpvt[1] = 2;
    work[1] = 0.0;
    temp = internal::blas::xnrm2(A.size(0), A, A.size(0) + 1);
    vn1[1] = temp;
    vn2[1] = temp;
    for (int b_i{0}; b_i < u0; b_i++) {
      double atmp;
      double beta1;
      int ii;
      int ip1;
      int lastv;
      int mmi;
      int pvt;
      ip1 = b_i + 2;
      knt = b_i * ma;
      ii = knt + b_i;
      mmi = m - b_i;
      ix = 0;
      if ((2 - b_i > 1) && (std::abs(vn1[1]) > std::abs(vn1[b_i]))) {
        ix = 1;
      }
      pvt = b_i + ix;
      if (pvt != b_i) {
        ix = pvt * ma;
        for (lastv = 0; lastv < m; lastv++) {
          rankA = ix + lastv;
          temp = b_A[rankA];
          i = knt + lastv;
          b_A[rankA] = b_A[i];
          b_A[i] = temp;
        }
        ix = jpvt[pvt];
        jpvt[pvt] = jpvt[b_i];
        jpvt[b_i] = static_cast<signed char>(ix);
        vn1[pvt] = vn1[b_i];
        vn2[pvt] = vn2[b_i];
      }
      if (b_i + 1 < m) {
        atmp = b_A[ii];
        ix = ii + 2;
        tau_data[b_i] = 0.0;
        if (mmi > 0) {
          temp = internal::blas::xnrm2(mmi - 1, b_A, ii + 2);
          if (temp != 0.0) {
            beta1 = rt_hypotd_snf(b_A[ii], temp);
            if (b_A[ii] >= 0.0) {
              beta1 = -beta1;
            }
            if (std::abs(beta1) < 1.0020841800044864E-292) {
              knt = 0;
              i = ii + mmi;
              do {
                knt++;
                for (lastv = ix; lastv <= i; lastv++) {
                  b_A[lastv - 1] = 9.9792015476736E+291 * b_A[lastv - 1];
                }
                beta1 *= 9.9792015476736E+291;
                atmp *= 9.9792015476736E+291;
              } while ((std::abs(beta1) < 1.0020841800044864E-292) &&
                       (knt < 20));
              beta1 = rt_hypotd_snf(
                  atmp, internal::blas::xnrm2(mmi - 1, b_A, ii + 2));
              if (atmp >= 0.0) {
                beta1 = -beta1;
              }
              tau_data[b_i] = (beta1 - atmp) / beta1;
              temp = 1.0 / (atmp - beta1);
              for (lastv = ix; lastv <= i; lastv++) {
                b_A[lastv - 1] = temp * b_A[lastv - 1];
              }
              for (lastv = 0; lastv < knt; lastv++) {
                beta1 *= 1.0020841800044864E-292;
              }
              atmp = beta1;
            } else {
              tau_data[b_i] = (beta1 - b_A[ii]) / beta1;
              temp = 1.0 / (b_A[ii] - beta1);
              i = ii + mmi;
              for (lastv = ix; lastv <= i; lastv++) {
                b_A[lastv - 1] = temp * b_A[lastv - 1];
              }
              atmp = beta1;
            }
          }
        }
        b_A[ii] = atmp;
      } else {
        tau_data[b_i] = 0.0;
      }
      if (b_i + 1 < 2) {
        int jA;
        atmp = b_A[ii];
        b_A[ii] = 1.0;
        jA = (ii + ma) + 1;
        if (tau_data[0] != 0.0) {
          lastv = mmi - 1;
          knt = (ii + mmi) - 1;
          while ((lastv + 1 > 0) && (b_A[knt] == 0.0)) {
            lastv--;
            knt--;
          }
          knt = 1;
          ix = jA;
          int exitg1;
          do {
            exitg1 = 0;
            if (ix <= jA + lastv) {
              if (b_A[ix - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ix++;
              }
            } else {
              knt = 0;
              exitg1 = 1;
            }
          } while (exitg1 == 0);
        } else {
          lastv = -1;
          knt = 0;
        }
        if (lastv + 1 > 0) {
          if (knt != 0) {
            work[0] = 0.0;
            rankA = 0;
            for (pvt = jA; ma < 0 ? pvt >= jA : pvt <= jA; pvt += ma) {
              temp = 0.0;
              i = pvt + lastv;
              for (ix = pvt; ix <= i; ix++) {
                temp += b_A[ix - 1] * b_A[(ii + ix) - pvt];
              }
              work[rankA] += temp;
              rankA++;
            }
          }
          if (!(-tau_data[0] == 0.0)) {
            for (ix = 0; ix < knt; ix++) {
              if (work[0] != 0.0) {
                temp = work[0] * -tau_data[0];
                i = lastv + jA;
                for (rankA = jA; rankA <= i; rankA++) {
                  b_A[rankA - 1] =
                      b_A[rankA - 1] + b_A[(ii + rankA) - jA] * temp;
                }
              }
              jA += ma;
            }
          }
        }
        b_A[ii] = atmp;
      }
      for (ix = ip1; ix < 3; ix++) {
        knt = b_i + ma;
        if (vn1[1] != 0.0) {
          temp = std::abs(b_A[knt]) / vn1[1];
          temp = 1.0 - temp * temp;
          if (temp < 0.0) {
            temp = 0.0;
          }
          beta1 = vn1[1] / vn2[1];
          beta1 = temp * (beta1 * beta1);
          if (beta1 <= 1.4901161193847656E-8) {
            if (b_i + 1 < m) {
              temp = internal::blas::xnrm2(mmi - 1, b_A, knt + 2);
              vn1[1] = temp;
              vn2[1] = temp;
            } else {
              vn1[1] = 0.0;
              vn2[1] = 0.0;
            }
          } else {
            vn1[1] *= std::sqrt(temp);
          }
        }
      }
    }
    rankA = 0;
    if (b_A.size(0) < 2) {
      knt = 1;
      ix = 2;
    } else {
      knt = 2;
      ix = b_A.size(0);
    }
    temp = std::fmin(1.4901161193847656E-8,
                     2.2204460492503131E-15 * static_cast<double>(ix)) *
           std::abs(b_A[0]);
    while ((rankA < knt) &&
           (!(std::abs(b_A[rankA + b_A.size(0) * rankA]) <= temp))) {
      rankA++;
    }
    b_B.set_size(B.size(0));
    knt = B.size(0);
    for (i = 0; i < knt; i++) {
      b_B[i] = B[i];
    }
    Y[0] = 0.0;
    Y[1] = 0.0;
    m = b_A.size(0);
    i = (b_A.size(0) >= 2);
    for (ix = 0; ix <= i; ix++) {
      if (tau_data[ix] != 0.0) {
        temp = b_B[ix];
        knt = ix + 2;
        for (int b_i{knt}; b_i <= m; b_i++) {
          temp += b_A[(b_i + b_A.size(0) * ix) - 1] * b_B[b_i - 1];
        }
        temp *= tau_data[ix];
        if (temp != 0.0) {
          b_B[ix] = b_B[ix] - temp;
          for (int b_i{knt}; b_i <= m; b_i++) {
            b_B[b_i - 1] =
                b_B[b_i - 1] - b_A[(b_i + b_A.size(0) * ix) - 1] * temp;
          }
        }
      }
    }
    for (int b_i{0}; b_i < rankA; b_i++) {
      Y[jpvt[b_i] - 1] = b_B[b_i];
    }
    for (ix = rankA; ix >= 1; ix--) {
      knt = jpvt[ix - 1] - 1;
      Y[knt] /= b_A[(ix + b_A.size(0) * (ix - 1)) - 1];
      for (int b_i{0}; b_i <= ix - 2; b_i++) {
        Y[jpvt[0] - 1] -= Y[knt] * b_A[b_A.size(0) * (ix - 1)];
      }
    }
  }
}

static void randperm(double n, double p[2])
{
  p[1] = 0.0;
  if (n <= 2.0) {
    double j;
    p[0] = 1.0;
    j = b_rand() * 2.0;
    j = std::floor(j);
    p[1] = p[static_cast<int>(j + 1.0) - 1];
    p[static_cast<int>(j + 1.0) - 1] = 2.0;
  } else if (n / 4.0 <= 2.0) {
    double j;
    double loc_idx_0;
    double selectedLoc;
    double val_idx_0;
    val_idx_0 = 0.0;
    j = n;
    selectedLoc = 2.0 / n;
    loc_idx_0 = b_rand();
    while (loc_idx_0 > selectedLoc) {
      val_idx_0++;
      j--;
      selectedLoc += (1.0 - selectedLoc) * (2.0 / j);
    }
    val_idx_0++;
    j = b_rand();
    j = std::floor(j);
    p[0] = 0.0;
    p[static_cast<int>(j + 1.0) - 1] = val_idx_0;
    j = n - val_idx_0;
    selectedLoc = 1.0 / j;
    loc_idx_0 = b_rand();
    while (loc_idx_0 > selectedLoc) {
      val_idx_0++;
      j--;
      selectedLoc += (1.0 - selectedLoc) * (1.0 / j);
    }
    val_idx_0++;
    j = b_rand() * 2.0;
    j = std::floor(j);
    p[1] = p[static_cast<int>(j + 1.0) - 1];
    p[static_cast<int>(j + 1.0) - 1] = val_idx_0;
  } else {
    double j;
    double loc_idx_0;
    double selectedLoc;
    double val_idx_0;
    int jlast;
    signed char hashTbl[2];
    hashTbl[0] = 0;
    hashTbl[1] = 0;
    selectedLoc = b_rand() * ((n - 1.0) + 1.0);
    selectedLoc = std::floor(selectedLoc);
    if (selectedLoc == 0.0) {
      j = 0.0;
    } else {
      j = std::fmod(selectedLoc, 2.0);
      if (j == 0.0) {
        j = 0.0;
      }
    }
    p[0] = selectedLoc + 1.0;
    loc_idx_0 = selectedLoc;
    hashTbl[static_cast<int>(j + 1.0) - 1] = 1;
    jlast = hashTbl[static_cast<int>(std::fmod(n - 1.0, 2.0) + 1.0) - 1];
    while ((jlast > 0) && (selectedLoc != n - 1.0)) {
      jlast = 0;
    }
    if (jlast > 0) {
      val_idx_0 = 0.0;
    } else {
      val_idx_0 = n - 1.0;
    }
    selectedLoc = b_rand() * ((n - 2.0) + 1.0);
    selectedLoc = std::floor(selectedLoc);
    if (selectedLoc == 0.0) {
      j = 0.0;
    } else {
      j = std::fmod(selectedLoc, 2.0);
      if (j == 0.0) {
        j = 0.0;
      }
    }
    j = hashTbl[static_cast<int>(j + 1.0) - 1];
    while ((j > 0.0) && (loc_idx_0 != selectedLoc)) {
      j = 0.0;
    }
    if (j > 0.0) {
      p[1] = val_idx_0 + 1.0;
    } else {
      p[1] = selectedLoc + 1.0;
    }
  }
}

} // namespace coder
static void eml_rand_mt19937ar_stateful_init()
{
  unsigned int r;
  std::memset(&state[0], 0, 625U * sizeof(unsigned int));
  r = 5489U;
  state[0] = 5489U;
  for (int mti{0}; mti < 623; mti++) {
    r = ((r ^ r >> 30U) * 1812433253U + static_cast<unsigned int>(mti)) + 1U;
    state[mti + 1] = r;
  }
  state[624] = 624U;
}

static double evaluateModel(const double modelIn[6],
                            const ::coder::array<double, 2U> &matchedPoints1,
                            const ::coder::array<double, 2U> &matchedPoints2,
                            const ::coder::array<double, 2U> &matchedLines1,
                            const ::coder::array<double, 2U> &matchedLines2,
                            ::coder::array<double, 1U> &distances)
{
  ::coder::array<double, 2U> B;
  ::coder::array<double, 2U> b_matchedLines1;
  ::coder::array<double, 2U> calPt2;
  ::coder::array<double, 2U> delta;
  ::coder::array<double, 2U> varargin_1;
  ::coder::array<double, 1U> D;
  ::coder::array<double, 1U> b_x;
  ::coder::array<double, 1U> dist1;
  ::coder::array<double, 1U> dist2;
  ::coder::array<double, 1U> r;
  ::coder::array<double, 1U> r1;
  ::coder::array<double, 1U> r2;
  ::coder::array<double, 1U> x;
  int boffset;
  int coffset;
  int j;
  int nx;
  dist1.set_size(matchedPoints1.size(0));
  nx = matchedPoints1.size(0);
  for (boffset = 0; boffset < nx; boffset++) {
    dist1[boffset] = 0.0;
  }
  if (matchedPoints1.size(0) != 0) {
    varargin_1.set_size(2, matchedPoints1.size(0));
    nx = matchedPoints1.size(0);
    for (boffset = 0; boffset < nx; boffset++) {
      varargin_1[2 * boffset] = matchedPoints1[boffset];
      varargin_1[2 * boffset + 1] =
          matchedPoints1[boffset + matchedPoints1.size(0)];
    }
    B.set_size(3, varargin_1.size(1));
    nx = varargin_1.size(1);
    for (boffset = 0; boffset < nx; boffset++) {
      for (j = 0; j < 2; j++) {
        B[j + B.size(0) * boffset] = varargin_1[j + 2 * boffset];
      }
      B[B.size(0) * boffset + 2] = 1.0;
    }
    nx = B.size(1);
    varargin_1.set_size(2, B.size(1));
    for (j = 0; j < nx; j++) {
      coffset = j << 1;
      boffset = j * 3;
      for (int i{0}; i < 2; i++) {
        varargin_1[coffset + i] =
            (modelIn[i] * B[boffset] + modelIn[i + 2] * B[boffset + 1]) +
            modelIn[i + 4] * B[boffset + 2];
      }
    }
    if (matchedPoints2.size(0) == varargin_1.size(1)) {
      delta.set_size(varargin_1.size(1), 2);
      nx = varargin_1.size(1);
      for (boffset = 0; boffset < 2; boffset++) {
        for (j = 0; j < nx; j++) {
          delta[j + delta.size(0) * boffset] =
              varargin_1[boffset + 2 * j] -
              matchedPoints2[j + matchedPoints2.size(0) * boffset];
        }
      }
    } else {
      binary_expand_op(delta, varargin_1, matchedPoints2);
    }
    dist1.set_size(delta.size(0));
    nx = delta.size(0);
    for (coffset = 0; coffset < nx; coffset++) {
      dist1[coffset] =
          rt_hypotd_snf(delta[coffset], delta[coffset + delta.size(0)]);
    }
    nx = dist1.size(0) - 1;
    for (int i{0}; i <= nx; i++) {
      if (dist1[i] > 1.5) {
        dist1[i] = 1.5;
      }
    }
  }
  dist2.set_size(matchedLines1.size(0));
  nx = matchedLines1.size(0);
  for (boffset = 0; boffset < nx; boffset++) {
    dist2[boffset] = 0.0;
  }
  if (matchedLines1.size(0) != 0) {
    double b_varargin_1;
    double c_varargin_1;
    unsigned int varargin_2;
    signed char input_sizes_idx_0;
    varargin_2 = static_cast<unsigned int>(matchedLines1.size(0)) << 1;
    if (static_cast<int>(varargin_2) != 0) {
      coffset = static_cast<int>(varargin_2);
    } else {
      coffset = matchedLines1.size(0) << 1;
    }
    if ((coffset == 0) || (static_cast<int>(varargin_2) != 0)) {
      input_sizes_idx_0 = 2;
    } else {
      input_sizes_idx_0 = 0;
    }
    b_matchedLines1.set_size(4, matchedLines1.size(0));
    nx = matchedLines1.size(0);
    for (boffset = 0; boffset < nx; boffset++) {
      b_matchedLines1[4 * boffset] = matchedLines1[boffset];
      b_matchedLines1[4 * boffset + 1] =
          matchedLines1[boffset + matchedLines1.size(0)];
      b_matchedLines1[4 * boffset + 2] =
          matchedLines1[boffset + matchedLines1.size(0) * 2];
      b_matchedLines1[4 * boffset + 3] =
          matchedLines1[boffset + matchedLines1.size(0) * 3];
    }
    nx = input_sizes_idx_0;
    B.set_size(input_sizes_idx_0 + 1, coffset);
    for (boffset = 0; boffset < coffset; boffset++) {
      for (j = 0; j < nx; j++) {
        B[j + B.size(0) * boffset] =
            b_matchedLines1[j + input_sizes_idx_0 * boffset];
      }
      B[input_sizes_idx_0 + B.size(0) * boffset] = 1.0;
    }
    nx = B.size(1);
    varargin_1.set_size(2, B.size(1));
    for (j = 0; j < nx; j++) {
      coffset = j << 1;
      boffset = j * B.size(0);
      for (int i{0}; i < 2; i++) {
        varargin_1[coffset + i] =
            (modelIn[i] * B[boffset] + modelIn[i + 2] * B[boffset + 1]) +
            modelIn[i + 4] * B[boffset + 2];
      }
    }
    calPt2.set_size(matchedLines1.size(0), 4);
    nx = matchedLines1.size(0);
    for (boffset = 0; boffset < 4; boffset++) {
      for (j = 0; j < nx; j++) {
        calPt2[j + calPt2.size(0) * boffset] = varargin_1[boffset + 4 * j];
      }
    }
    x.set_size(matchedLines2.size(0));
    nx = matchedLines2.size(0);
    for (boffset = 0; boffset < nx; boffset++) {
      x[boffset] = matchedLines2[boffset + matchedLines2.size(0) * 3] -
                   matchedLines2[boffset + matchedLines2.size(0)];
    }
    if (x.size(0) == matchedLines2.size(0)) {
      D.set_size(x.size(0));
      nx = x.size(0);
      for (boffset = 0; boffset < nx; boffset++) {
        b_varargin_1 = x[boffset];
        c_varargin_1 = matchedLines2[boffset + matchedLines2.size(0) * 2] -
                       matchedLines2[boffset];
        D[boffset] = b_varargin_1 * b_varargin_1 + c_varargin_1 * c_varargin_1;
      }
    } else {
      binary_expand_op(D, x, matchedLines2);
    }
    nx = D.size(0);
    for (coffset = 0; coffset < nx; coffset++) {
      D[coffset] = std::sqrt(D[coffset]);
    }
    r.set_size(matchedLines2.size(0));
    nx = matchedLines2.size(0);
    r1.set_size(matchedLines2.size(0));
    r2.set_size(matchedLines2.size(0));
    for (boffset = 0; boffset < nx; boffset++) {
      b_varargin_1 = matchedLines2[boffset];
      c_varargin_1 = matchedLines2[boffset + matchedLines2.size(0) * 2];
      r[boffset] = b_varargin_1 - c_varargin_1;
      r1[boffset] =
          c_varargin_1 * matchedLines2[boffset + matchedLines2.size(0)];
      r2[boffset] =
          b_varargin_1 * matchedLines2[boffset + matchedLines2.size(0) * 3];
    }
    if (x.size(0) == 1) {
      boffset = calPt2.size(0);
    } else {
      boffset = x.size(0);
    }
    if (r.size(0) == 1) {
      j = calPt2.size(0);
    } else {
      j = r.size(0);
    }
    if (boffset == 1) {
      nx = j;
    } else {
      nx = boffset;
    }
    if (nx == 1) {
      coffset = r1.size(0);
    } else {
      coffset = nx;
    }
    if ((x.size(0) == calPt2.size(0)) && (r.size(0) == calPt2.size(0)) &&
        (boffset == j) && (nx == r1.size(0)) && (coffset == r2.size(0))) {
      b_x.set_size(x.size(0));
      nx = x.size(0);
      for (boffset = 0; boffset < nx; boffset++) {
        b_x[boffset] = ((x[boffset] * calPt2[boffset] +
                         r[boffset] * calPt2[boffset + calPt2.size(0)]) +
                        r1[boffset]) -
                       r2[boffset];
      }
    } else {
      binary_expand_op(b_x, x, calPt2, r, r1, r2);
    }
    nx = b_x.size(0);
    dist2.set_size(b_x.size(0));
    for (coffset = 0; coffset < nx; coffset++) {
      dist2[coffset] = std::abs(b_x[coffset]);
    }
    if (x.size(0) == 1) {
      boffset = calPt2.size(0);
    } else {
      boffset = x.size(0);
    }
    if (r.size(0) == 1) {
      j = calPt2.size(0);
    } else {
      j = r.size(0);
    }
    if (boffset == 1) {
      nx = j;
    } else {
      nx = boffset;
    }
    if (nx == 1) {
      coffset = r1.size(0);
    } else {
      coffset = nx;
    }
    if ((x.size(0) == calPt2.size(0)) && (r.size(0) == calPt2.size(0)) &&
        (boffset == j) && (nx == r1.size(0)) && (coffset == r2.size(0))) {
      nx = x.size(0);
      for (boffset = 0; boffset < nx; boffset++) {
        x[boffset] = ((x[boffset] * calPt2[boffset + calPt2.size(0) * 2] +
                       r[boffset] * calPt2[boffset + calPt2.size(0) * 3]) +
                      r1[boffset]) -
                     r2[boffset];
      }
    } else {
      binary_expand_op(x, calPt2, r, r1, r2);
    }
    nx = x.size(0);
    b_x.set_size(x.size(0));
    for (coffset = 0; coffset < nx; coffset++) {
      b_x[coffset] = std::abs(x[coffset]);
    }
    if (dist2.size(0) == 1) {
      coffset = D.size(0);
    } else {
      coffset = dist2.size(0);
    }
    if (b_x.size(0) == 1) {
      nx = D.size(0);
    } else {
      nx = b_x.size(0);
    }
    if ((dist2.size(0) == D.size(0)) && (b_x.size(0) == D.size(0)) &&
        (coffset == nx)) {
      nx = dist2.size(0);
      for (boffset = 0; boffset < nx; boffset++) {
        dist2[boffset] =
            (dist2[boffset] / D[boffset] + b_x[boffset] / D[boffset]) / 2.0 *
            0.5;
      }
    } else {
      binary_expand_op(dist2, D, b_x);
    }
    nx = dist2.size(0) - 1;
    for (int i{0}; i <= nx; i++) {
      if (dist2[i] > 1.5) {
        dist2[i] = 1.5;
      }
    }
  }
  distances.set_size(dist1.size(0) + dist2.size(0));
  nx = dist1.size(0);
  for (boffset = 0; boffset < nx; boffset++) {
    distances[boffset] = dist1[boffset];
  }
  nx = dist2.size(0);
  for (boffset = 0; boffset < nx; boffset++) {
    distances[boffset + dist1.size(0)] = dist2[boffset];
  }
  return coder::combineVectorElements(distances);
}

static bool judgeLinesValid(const ::coder::array<double, 2U> &lines)
{
  int i;
  bool flag;
  flag = false;
  i = lines.size(0);
  for (int b_i{0}; b_i <= i - 2; b_i++) {
    double a_idx_0;
    double a_idx_1;
    double b_idx_0;
    double b_idx_1;
    a_idx_0 = lines[b_i + lines.size(0) * 2] - lines[b_i];
    a_idx_1 = lines[b_i + lines.size(0) * 3] - lines[b_i + lines.size(0)];
    b_idx_0 = lines[(b_i + lines.size(0) * 2) + 1] - lines[b_i + 1];
    b_idx_1 =
        lines[(b_i + lines.size(0) * 3) + 1] - lines[(b_i + lines.size(0)) + 1];
    flag = ((57.295779513082323 *
                 std::acos((a_idx_0 * b_idx_0 + a_idx_1 * b_idx_1) /
                           (std::sqrt(a_idx_0 * a_idx_0 + a_idx_1 * a_idx_1) *
                            std::sqrt(b_idx_0 * b_idx_0 + b_idx_1 * b_idx_1))) >
             5.0) ||
            flag);
  }
  return flag;
}

static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaN;
  } else if (std::isinf(u0) && std::isinf(u1)) {
    int i;
    int i1;
    if (u0 > 0.0) {
      i = 1;
    } else {
      i = -1;
    }
    if (u1 > 0.0) {
      i1 = 1;
    } else {
      i1 = -1;
    }
    y = std::atan2(static_cast<double>(i), static_cast<double>(i1));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }
  return y;
}

static double rt_hypotd_snf(double u0, double u1)
{
  double a;
  double b;
  double y;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0);
  } else if (std::isnan(b)) {
    y = rtNaN;
  } else {
    y = a * 1.4142135623730951;
  }
  return y;
}

static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaN;
  } else {
    double d;
    double d1;
    d = std::abs(u0);
    d1 = std::abs(u1);
    if (std::isinf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = std::pow(u0, u1);
    }
  }
  return y;
}

static double rt_remd_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1) || std::isinf(u0)) {
    y = rtNaN;
  } else if (std::isinf(u1)) {
    y = u0;
  } else if ((u1 != 0.0) && (u1 != std::trunc(u1))) {
    double q;
    q = std::abs(u0 / u1);
    if (!(std::abs(q - std::floor(q + 0.5)) > DBL_EPSILON * q)) {
      y = 0.0 * u0;
    } else {
      y = std::fmod(u0, u1);
    }
  } else {
    y = std::fmod(u0, u1);
  }
  return y;
}

static void
svdFitForInliersPointsAndLines(const ::coder::array<double, 2U> &inliersPts1,
                               const ::coder::array<double, 2U> &inliersPts2,
                               const ::coder::array<double, 2U> &inliersL1,
                               const ::coder::array<double, 2U> &inliersL2,
                               double tform[6])
{
  ::coder::array<double, 2U> b_inliersPts1;
  ::coder::array<double, 2U> normPoints1;
  ::coder::array<double, 2U> normPoints2;
  ::coder::array<double, 2U> t1;
  ::coder::array<double, 1U> b_x;
  ::coder::array<double, 1U> x;
  double a[2];
  double centroid1[2];
  double centroid2[2];
  double bkj;
  double d;
  double d1;
  double d2;
  double d3;
  double theta1;
  double theta2;
  int acoef;
  int boffset;
  bool t1_flag;
  bool t2_flag;
  bool theta1_flag;
  bool theta2_flag;
  theta1 = 0.0;
  t1.set_size(2, 1);
  t1[0] = 0.0;
  t1[1] = 0.0;
  theta1_flag = false;
  t1_flag = false;
  if (inliersPts1.size(0) >= 2) {
    double C[4];
    double U[4];
    double V[4];
    int coffset;
    signed char csz[2];
    coder::mean(inliersPts1, centroid1);
    coder::mean(inliersPts2, centroid2);
    normPoints1.set_size(inliersPts1.size(0), 2);
    for (int k{0}; k < 2; k++) {
      boffset = normPoints1.size(0) - 1;
      for (coffset = 0; coffset <= boffset; coffset++) {
        normPoints1[coffset + normPoints1.size(0) * k] =
            inliersPts1[coffset + inliersPts1.size(0) * k] - centroid1[k];
      }
    }
    normPoints2.set_size(inliersPts2.size(0), 2);
    if (inliersPts2.size(0) != 0) {
      acoef = (inliersPts2.size(0) != 1);
      for (int k{0}; k < 2; k++) {
        boffset = normPoints2.size(0) - 1;
        for (coffset = 0; coffset <= boffset; coffset++) {
          normPoints2[coffset + normPoints2.size(0) * k] =
              inliersPts2[acoef * coffset + inliersPts2.size(0) * k] -
              centroid2[k];
        }
      }
    }
    acoef = normPoints1.size(0);
    for (int j{0}; j < 2; j++) {
      coffset = j << 1;
      boffset = j * normPoints2.size(0);
      C[coffset] = 0.0;
      C[coffset + 1] = 0.0;
      for (int k{0}; k < acoef; k++) {
        bkj = normPoints2[boffset + k];
        C[coffset] += normPoints1[k] * bkj;
        C[coffset + 1] += normPoints1[normPoints1.size(0) + k] * bkj;
      }
    }
    theta1_flag = true;
    if (std::isinf(C[0]) || std::isnan(C[0]) ||
        (std::isinf(C[1]) || std::isnan(C[1]))) {
      theta1_flag = false;
    }
    if ((!theta1_flag) || (std::isinf(C[2]) || std::isnan(C[2]))) {
      theta1_flag = false;
    }
    if ((!theta1_flag) || (std::isinf(C[3]) || std::isnan(C[3]))) {
      theta1_flag = false;
    }
    if (theta1_flag) {
      coder::internal::svd(C, U, a, V);
    } else {
      U[0] = rtNaN;
      V[0] = rtNaN;
      U[1] = rtNaN;
      V[1] = rtNaN;
      U[2] = rtNaN;
      V[2] = rtNaN;
      U[3] = rtNaN;
      V[3] = rtNaN;
    }
    d = V[0];
    d1 = V[2];
    d2 = V[1];
    d3 = V[3];
    for (boffset = 0; boffset < 2; boffset++) {
      bkj = U[boffset + 2];
      theta2 = U[boffset];
      C[boffset] = theta2 * d + bkj * d1;
      C[boffset + 2] = theta2 * d2 + bkj * d3;
      csz[boffset] = static_cast<signed char>(boffset + 1);
    }
    acoef = 0;
    if (std::abs(C[1]) > std::abs(C[0])) {
      acoef = 1;
    }
    if (C[acoef] != 0.0) {
      if (acoef != 0) {
        csz[0] = 2;
        bkj = C[0];
        C[0] = C[1];
        C[1] = bkj;
        bkj = C[2];
        C[2] = C[3];
        C[3] = bkj;
      }
      C[1] /= C[0];
    }
    if (C[2] != 0.0) {
      C[3] += C[1] * -C[2];
    }
    bkj = C[0] * C[3];
    if (csz[0] > 1) {
      bkj = -bkj;
    }
    if (std::isnan(bkj)) {
      C[3] = rtNaN;
    } else if (bkj < 0.0) {
      C[3] = -1.0;
    } else {
      C[3] = (bkj > 0.0);
    }
    t1.set_size(2, 1);
    d = U[0];
    d1 = U[2];
    d2 = U[1];
    d3 = U[3];
    bkj = centroid1[0];
    theta2 = centroid1[1];
    for (boffset = 0; boffset < 2; boffset++) {
      double d4;
      double d5;
      double d6;
      theta1 = V[boffset + 2];
      d4 = V[boffset];
      d5 = d4 + theta1 * 0.0;
      theta1 = d4 * 0.0 + theta1 * C[3];
      d4 = d5 * d + theta1 * d1;
      C[boffset] = d4;
      d6 = d4 * bkj;
      d4 = d5 * d2 + theta1 * d3;
      C[boffset + 2] = d4;
      d6 += d4 * theta2;
      t1[boffset] = centroid2[boffset] - d6;
    }
    bkj = rt_atan2d_snf(C[1], C[0]);
    if (C[0] * C[0] + C[1] * C[1] < 2.2204460492503131E-15) {
      bkj = 0.0;
    }
    theta1 = bkj * 180.0 / 3.1415926535897931;
    theta1_flag = true;
    t1_flag = true;
  } else if (inliersPts1.size(0) == 1) {
    b_inliersPts1.set_size(2, 1);
    b_inliersPts1[0] = inliersPts1[0];
    b_inliersPts1[1] = inliersPts1[inliersPts1.size(0)];
    t1.set_size(2, inliersPts2.size(0));
    acoef = inliersPts2.size(0);
    for (boffset = 0; boffset < acoef; boffset++) {
      t1[2 * boffset] = inliersPts2[boffset] - b_inliersPts1[0];
      t1[2 * boffset + 1] =
          inliersPts2[boffset + inliersPts2.size(0)] - b_inliersPts1[1];
    }
    t1_flag = true;
  }
  theta2 = 0.0;
  centroid2[0] = 0.0;
  centroid2[1] = 0.0;
  theta2_flag = false;
  t2_flag = false;
  if (inliersL1.size(0) >= 2) {
    normPoints1.set_size(inliersL1.size(0), 2);
    acoef = inliersL1.size(0);
    for (boffset = 0; boffset < acoef; boffset++) {
      normPoints1[boffset] =
          inliersL1[boffset + inliersL1.size(0) * 2] - inliersL1[boffset];
      normPoints1[boffset + normPoints1.size(0)] =
          inliersL1[boffset + inliersL1.size(0) * 3] -
          inliersL1[boffset + inliersL1.size(0)];
    }
    normPoints2.set_size(inliersL2.size(0), 2);
    acoef = inliersL2.size(0);
    for (boffset = 0; boffset < acoef; boffset++) {
      normPoints2[boffset] =
          inliersL2[boffset + inliersL2.size(0) * 2] - inliersL2[boffset];
      normPoints2[boffset + normPoints2.size(0)] =
          inliersL2[boffset + inliersL2.size(0) * 3] -
          inliersL2[boffset + inliersL2.size(0)];
    }
    x.set_size(normPoints1.size(0));
    acoef = normPoints1.size(0);
    for (boffset = 0; boffset < acoef; boffset++) {
      bkj = normPoints1[boffset];
      theta2 = normPoints1[boffset + normPoints1.size(0)];
      x[boffset] = bkj * bkj + theta2 * theta2;
    }
    acoef = x.size(0);
    for (int k{0}; k < acoef; k++) {
      x[k] = std::sqrt(x[k]);
    }
    b_x.set_size(normPoints2.size(0));
    acoef = normPoints2.size(0);
    for (boffset = 0; boffset < acoef; boffset++) {
      bkj = normPoints2[boffset];
      theta2 = normPoints2[boffset + normPoints2.size(0)];
      b_x[boffset] = bkj * bkj + theta2 * theta2;
    }
    acoef = b_x.size(0);
    for (int k{0}; k < acoef; k++) {
      b_x[k] = std::sqrt(b_x[k]);
    }
    if ((normPoints1.size(0) == normPoints2.size(0)) &&
        (x.size(0) == b_x.size(0)) && (normPoints1.size(0) == x.size(0))) {
      x.set_size(normPoints1.size(0));
      acoef = normPoints1.size(0);
      for (boffset = 0; boffset < acoef; boffset++) {
        x[boffset] = (normPoints1[boffset] * normPoints2[boffset] +
                      normPoints1[boffset + normPoints1.size(0)] *
                          normPoints2[boffset + normPoints2.size(0)]) /
                     (x[boffset] * b_x[boffset]);
      }
    } else {
      binary_expand_op(x, normPoints1, normPoints2, b_x);
    }
    acoef = x.size(0);
    for (int k{0}; k < acoef; k++) {
      x[k] = 57.295779513082323 * std::acos(x[k]);
    }
    theta2 = coder::combineVectorElements(x) / static_cast<double>(x.size(0));
    normPoints1.set_size(inliersL2.size(0), 2);
    acoef = inliersL2.size(0);
    x.set_size(inliersL2.size(0));
    for (boffset = 0; boffset < acoef; boffset++) {
      d = inliersL2[boffset + inliersL2.size(0) * 3];
      d1 = inliersL2[boffset + inliersL2.size(0)];
      normPoints1[boffset] = d - d1;
      d2 = inliersL2[boffset];
      d3 = inliersL2[boffset + inliersL2.size(0) * 2];
      normPoints1[boffset + normPoints1.size(0)] = d2 - d3;
      x[boffset] = d2 * d - d3 * d1;
    }
    coder::mldivide(normPoints1, x, centroid2);
    normPoints1.set_size(inliersL1.size(0), 2);
    acoef = inliersL1.size(0);
    x.set_size(inliersL1.size(0));
    for (boffset = 0; boffset < acoef; boffset++) {
      d = inliersL1[boffset + inliersL1.size(0) * 3];
      d1 = inliersL1[boffset + inliersL1.size(0)];
      normPoints1[boffset] = d - d1;
      d2 = inliersL1[boffset];
      d3 = inliersL1[boffset + inliersL1.size(0) * 2];
      normPoints1[boffset + normPoints1.size(0)] = d2 - d3;
      x[boffset] = d2 * d - d3 * d1;
    }
    coder::mldivide(normPoints1, x, a);
    centroid2[0] -= a[0];
    centroid2[1] -= a[1];
    theta2_flag = true;
    t2_flag = true;
  } else if (inliersL1.size(0) == 1) {
    a[0] = inliersL1[inliersL1.size(0) * 2] - inliersL1[0];
    a[1] = inliersL1[inliersL1.size(0) * 3] - inliersL1[inliersL1.size(0)];
    centroid1[0] = inliersL2[inliersL2.size(0) * 2] - inliersL2[0];
    centroid1[1] =
        inliersL2[inliersL2.size(0) * 3] - inliersL2[inliersL2.size(0)];
    theta2 = 57.295779513082323 *
             std::acos((a[0] * centroid1[0] + a[1] * centroid1[1]) /
                       (std::sqrt(a[0] * a[0] + a[1] * a[1]) *
                        std::sqrt(centroid1[0] * centroid1[0] +
                                  centroid1[1] * centroid1[1])));
    theta2_flag = true;
  }
  if (theta1_flag && theta2_flag) {
    theta2 = (theta1 + theta2) / 2.0;
  } else if (theta1_flag && (!theta2_flag)) {
    theta2 = theta1;
  } else if (theta1_flag || (!theta2_flag)) {
    theta2 = 0.0;
  }
  if (t1_flag && t2_flag) {
    b_inliersPts1.set_size(2, t1.size(1));
    acoef = t1.size(1);
    for (boffset = 0; boffset < acoef; boffset++) {
      b_inliersPts1[2 * boffset] = (t1[2 * boffset] + centroid2[0]) / 2.0;
      b_inliersPts1[2 * boffset + 1] =
          (t1[2 * boffset + 1] + centroid2[1]) / 2.0;
    }
    t1.set_size(2, b_inliersPts1.size(1));
    acoef = b_inliersPts1.size(1) << 1;
    for (boffset = 0; boffset < acoef; boffset++) {
      t1[boffset] = b_inliersPts1[boffset];
    }
  } else if ((!t1_flag) || t2_flag) {
    if ((!t1_flag) && t2_flag) {
      t1.set_size(2, 1);
      t1[0] = centroid2[0];
      t1[1] = centroid2[1];
    } else {
      t1.set_size(2, 1);
      t1[0] = 0.0;
      t1[1] = 0.0;
    }
  }
  bkj = theta2;
  coder::b_sind(bkj);
  coder::b_cosd(theta2);
  tform[0] = theta2;
  tform[2] = -bkj;
  tform[4] = t1[0];
  tform[1] = bkj;
  tform[3] = theta2;
  tform[5] = t1[1];
}

void estgeotform2dForPtsAndLines(
    const ::coder::array<double, 2U> &matchedPoints1,
    const ::coder::array<double, 2U> &matchedPoints2,
    const ::coder::array<double, 2U> &matchedLines1,
    const ::coder::array<double, 2U> &matchedLines2, double tform[6],
    ::coder::array<bool, 1U> &inlierIndex, double *status)
{
  static const signed char iv[6]{1, 0, 0, 1, 0, 0};
  ::coder::array<double, 2U> b_matchedLines1;
  ::coder::array<double, 2U> b_matchedLines2;
  ::coder::array<double, 2U> b_matchedPoints1;
  ::coder::array<double, 2U> b_matchedPoints2;
  ::coder::array<double, 2U> b_varargin_1;
  ::coder::array<double, 2U> calPt2;
  ::coder::array<double, 2U> delta;
  ::coder::array<double, 2U> varargin_1;
  ::coder::array<double, 1U> D;
  ::coder::array<double, 1U> b_delta;
  ::coder::array<double, 1U> dist1;
  ::coder::array<double, 1U> dist2;
  ::coder::array<double, 1U> r4;
  ::coder::array<double, 1U> r5;
  ::coder::array<double, 1U> r6;
  ::coder::array<double, 1U> x;
  ::coder::array<int, 1U> r;
  ::coder::array<int, 1U> r1;
  ::coder::array<int, 1U> r2;
  ::coder::array<int, 1U> r3;
  ::coder::array<bool, 1U> bestInliers;
  int b_status;
  int i;
  int loop_ub;
  int m;
  unsigned int u;
  bool cond2;
  bool cond3;
  if (!isInitialized_estgeotform2dForPtsAndLines) {
    estgeotform2dForPtsAndLines_initialize();
  }
  m = matchedPoints1.size(0);
  if ((matchedPoints1.size(0) == 1) && (matchedLines1.size(0) >= 1)) {
    cond2 = true;
  } else {
    cond2 = false;
  }
  cond3 = false;
  if ((matchedLines1.size(0) >= 2) &&
      (judgeLinesValid(matchedLines1) || judgeLinesValid(matchedLines2))) {
    cond3 = true;
  }
  for (i = 0; i < 6; i++) {
    tform[i] = iv[i];
  }
  u = static_cast<unsigned int>(matchedPoints1.size(0)) +
      static_cast<unsigned int>(matchedLines1.size(0));
  loop_ub = static_cast<int>(u);
  inlierIndex.set_size(static_cast<int>(u));
  for (i = 0; i < loop_ub; i++) {
    inlierIndex[i] = false;
  }
  b_status = 0;
  if ((matchedPoints1.size(0) >= 2) || cond2 || cond3) {
    double modelParams[6];
    double bestDis;
    double t_idx_0;
    double t_idx_1;
    int idxTrial;
    int k;
    int numTrials;
    int sizes_idx_1;
    int skipTrials;
    int vlen;
    bool bv[6];
    idxTrial = 1;
    numTrials = 1000;
    skipTrials = 0;
    bestDis = 1.5 * static_cast<double>(u);
    bestInliers.set_size(static_cast<int>(u));
    for (i = 0; i < loop_ub; i++) {
      bestInliers[i] = false;
    }
    while ((idxTrial <= numTrials) && (skipTrials < 10000)) {
      double indices[2];
      bool b_indices[2];
      coder::randperm(static_cast<double>(u), indices);
      b_indices[0] = (indices[0] <= m);
      b_indices[1] = (indices[1] <= m);
      if (coder::all(b_indices)) {
        double a;
        double a_idx_0;
        double a_idx_1;
        double b_a;
        double b_d_tmp;
        double b_idx_0;
        double b_idx_1;
        double c;
        double c_d_tmp;
        double d_d_tmp;
        double d_tmp;
        double rotateT;
        t_idx_1 = matchedPoints1[static_cast<int>(indices[0]) - 1];
        t_idx_0 = matchedPoints1[static_cast<int>(indices[1]) - 1];
        a_idx_0 = t_idx_0 - t_idx_1;
        d_tmp = matchedPoints1[(static_cast<int>(indices[0]) +
                                matchedPoints1.size(0)) -
                               1];
        b_d_tmp = matchedPoints1[(static_cast<int>(indices[1]) +
                                  matchedPoints1.size(0)) -
                                 1];
        a_idx_1 = b_d_tmp - d_tmp;
        c_d_tmp = matchedPoints2[static_cast<int>(indices[0]) - 1];
        d_d_tmp = matchedPoints2[static_cast<int>(indices[1]) - 1];
        b_idx_0 = d_d_tmp - c_d_tmp;
        a = matchedPoints2[(static_cast<int>(indices[0]) +
                            matchedPoints2.size(0)) -
                           1];
        b_a = matchedPoints2[(static_cast<int>(indices[1]) +
                              matchedPoints2.size(0)) -
                             1];
        b_idx_1 = b_a - a;
        rotateT = 57.295779513082323 *
                  std::acos((a_idx_0 * b_idx_0 + a_idx_1 * b_idx_1) /
                            (std::sqrt(a_idx_0 * a_idx_0 + a_idx_1 * a_idx_1) *
                             std::sqrt(b_idx_0 * b_idx_0 + b_idx_1 * b_idx_1)));
        c = rotateT;
        coder::b_sind(c);
        coder::b_cosd(rotateT);
        t_idx_0 = (t_idx_1 + t_idx_0) / 2.0;
        t_idx_1 = (d_tmp + b_d_tmp) / 2.0;
        modelParams[0] = rotateT;
        modelParams[1] = c;
        modelParams[4] =
            (c_d_tmp + d_d_tmp) / 2.0 - (rotateT * t_idx_0 + -c * t_idx_1);
        modelParams[2] = -c;
        modelParams[3] = rotateT;
        modelParams[5] = (a + b_a) / 2.0 - (c * t_idx_0 + rotateT * t_idx_1);
      } else {
        b_indices[0] = (indices[0] > m);
        b_indices[1] = (indices[1] > m);
        if (coder::all(b_indices)) {
          double a;
          double a_idx_0;
          double a_idx_1;
          double b_a;
          double b_d_tmp;
          double b_idx_0;
          double b_idx_1;
          double c;
          double c_d_tmp;
          double d;
          double d_d_tmp;
          double d_tmp;
          double e_d_tmp;
          double f_d_tmp;
          double g_d_tmp;
          double h_d_tmp;
          double i_d_tmp;
          double intesect_x1;
          double intesect_y1;
          double j_d_tmp;
          double k_d_tmp;
          double l_d_tmp;
          double m_d_tmp;
          double n_d_tmp;
          double o_d_tmp;
          double p_d_tmp;
          double q_d_tmp;
          double r_d_tmp;
          double rotateT;
          double s_d_tmp;
          vlen = static_cast<int>(indices[0] - static_cast<double>(m)) - 1;
          sizes_idx_1 =
              static_cast<int>(indices[1] - static_cast<double>(m)) - 1;
          rotateT = matchedLines1[sizes_idx_1];
          c_d_tmp = matchedLines1[vlen + matchedLines1.size(0)];
          e_d_tmp = matchedLines1[sizes_idx_1 + matchedLines1.size(0)];
          d_tmp = matchedLines1[vlen];
          d_d_tmp = matchedLines1[vlen + matchedLines1.size(0) * 3];
          f_d_tmp = matchedLines1[sizes_idx_1 + matchedLines1.size(0) * 2];
          b_d_tmp = matchedLines1[vlen + matchedLines1.size(0) * 2];
          g_d_tmp = matchedLines1[sizes_idx_1 + matchedLines1.size(0) * 3];
          h_d_tmp = d_tmp * e_d_tmp;
          i_d_tmp = b_d_tmp * e_d_tmp;
          j_d_tmp = d_tmp * g_d_tmp;
          k_d_tmp = b_d_tmp * g_d_tmp;
          d = ((((((c_d_tmp * rotateT - h_d_tmp) - d_d_tmp * rotateT) -
                  c_d_tmp * f_d_tmp) +
                 i_d_tmp) +
                j_d_tmp) +
               d_d_tmp * f_d_tmp) -
              k_d_tmp;
          t_idx_1 = rotateT * g_d_tmp;
          t_idx_0 = e_d_tmp * f_d_tmp;
          intesect_x1 =
              (((((((t_idx_1 * d_tmp - t_idx_0 * d_tmp) - t_idx_1 * b_d_tmp) +
                   t_idx_0 * b_d_tmp) -
                  d_tmp * rotateT * d_d_tmp) +
                 b_d_tmp * rotateT * c_d_tmp) +
                d_tmp * f_d_tmp * d_d_tmp) -
               b_d_tmp * f_d_tmp * c_d_tmp) /
              d;
          intesect_y1 = (((((((t_idx_1 * c_d_tmp - t_idx_0 * c_d_tmp) -
                              t_idx_1 * d_d_tmp) +
                             t_idx_0 * d_d_tmp) -
                            h_d_tmp * d_d_tmp) +
                           i_d_tmp * c_d_tmp) +
                          j_d_tmp * d_d_tmp) -
                         k_d_tmp * c_d_tmp) /
                        d;
          h_d_tmp = matchedLines2[sizes_idx_1];
          i_d_tmp = matchedLines2[vlen + matchedLines2.size(0)];
          j_d_tmp = matchedLines2[sizes_idx_1 + matchedLines2.size(0)];
          k_d_tmp = matchedLines2[vlen];
          l_d_tmp = matchedLines2[vlen + matchedLines2.size(0) * 3];
          m_d_tmp = matchedLines2[sizes_idx_1 + matchedLines2.size(0) * 2];
          n_d_tmp = matchedLines2[vlen + matchedLines2.size(0) * 2];
          o_d_tmp = matchedLines2[sizes_idx_1 + matchedLines2.size(0) * 3];
          p_d_tmp = k_d_tmp * j_d_tmp;
          q_d_tmp = n_d_tmp * j_d_tmp;
          r_d_tmp = k_d_tmp * o_d_tmp;
          s_d_tmp = n_d_tmp * o_d_tmp;
          d = ((((((i_d_tmp * h_d_tmp - p_d_tmp) - l_d_tmp * h_d_tmp) -
                  i_d_tmp * m_d_tmp) +
                 q_d_tmp) +
                r_d_tmp) +
               l_d_tmp * m_d_tmp) -
              s_d_tmp;
          a_idx_0 = b_d_tmp - d_tmp;
          a_idx_1 = d_d_tmp - c_d_tmp;
          b_idx_0 = n_d_tmp - k_d_tmp;
          b_idx_1 = l_d_tmp - i_d_tmp;
          t_idx_1 = a_idx_0;
          t_idx_0 = a_idx_1;
          a = b_idx_0;
          b_a = b_idx_1;
          c = a_idx_0 * b_idx_0 + a_idx_1 * b_idx_1;
          a_idx_0 = f_d_tmp - rotateT;
          a_idx_1 = g_d_tmp - e_d_tmp;
          b_idx_0 = m_d_tmp - h_d_tmp;
          b_idx_1 = o_d_tmp - j_d_tmp;
          rotateT =
              (57.295779513082323 *
                   std::acos(c /
                             (std::sqrt(t_idx_1 * t_idx_1 + t_idx_0 * t_idx_0) *
                              std::sqrt(a * a + b_a * b_a))) +
               57.295779513082323 *
                   std::acos(
                       (a_idx_0 * b_idx_0 + a_idx_1 * b_idx_1) /
                       (std::sqrt(a_idx_0 * a_idx_0 + a_idx_1 * a_idx_1) *
                        std::sqrt(b_idx_0 * b_idx_0 + b_idx_1 * b_idx_1)))) /
              2.0;
          c = rotateT;
          coder::b_sind(c);
          coder::b_cosd(rotateT);
          modelParams[0] = rotateT;
          modelParams[2] = -c;
          t_idx_0 = h_d_tmp * o_d_tmp;
          t_idx_1 = j_d_tmp * m_d_tmp;
          modelParams[4] = (((((((t_idx_0 * k_d_tmp - t_idx_1 * k_d_tmp) -
                                 t_idx_0 * n_d_tmp) +
                                t_idx_1 * n_d_tmp) -
                               k_d_tmp * h_d_tmp * l_d_tmp) +
                              n_d_tmp * h_d_tmp * i_d_tmp) +
                             k_d_tmp * m_d_tmp * l_d_tmp) -
                            n_d_tmp * m_d_tmp * i_d_tmp) /
                               d -
                           intesect_x1;
          modelParams[1] = c;
          modelParams[3] = rotateT;
          modelParams[5] = (((((((t_idx_0 * i_d_tmp - t_idx_1 * i_d_tmp) -
                                 t_idx_0 * l_d_tmp) +
                                t_idx_1 * l_d_tmp) -
                               p_d_tmp * l_d_tmp) +
                              q_d_tmp * i_d_tmp) +
                             r_d_tmp * l_d_tmp) -
                            s_d_tmp * i_d_tmp) /
                               d -
                           intesect_y1;
        } else {
          double a_idx_0;
          double a_idx_1;
          double b_idx_0;
          double b_idx_1;
          double rotateT;
          coder::internal::sort(indices);
          vlen = static_cast<int>(indices[1] - static_cast<double>(m)) - 1;
          a_idx_0 = matchedLines1[vlen + matchedLines1.size(0) * 2] -
                    matchedLines1[vlen];
          a_idx_1 = matchedLines1[vlen + matchedLines1.size(0) * 3] -
                    matchedLines1[vlen + matchedLines1.size(0)];
          b_idx_0 = matchedLines2[vlen + matchedLines2.size(0) * 2] -
                    matchedLines2[vlen];
          b_idx_1 = matchedLines2[vlen + matchedLines2.size(0) * 3] -
                    matchedLines2[vlen + matchedLines2.size(0)];
          rotateT =
              57.295779513082323 *
              std::acos((a_idx_0 * b_idx_0 + a_idx_1 * b_idx_1) /
                        (std::sqrt(a_idx_0 * a_idx_0 + a_idx_1 * a_idx_1) *
                         std::sqrt(b_idx_0 * b_idx_0 + b_idx_1 * b_idx_1)));
          t_idx_1 = rotateT;
          coder::b_sind(t_idx_1);
          coder::b_cosd(rotateT);
          modelParams[0] = rotateT;
          modelParams[2] = -t_idx_1;
          modelParams[1] = t_idx_1;
          modelParams[3] = rotateT;
          modelParams[4] = matchedPoints2[static_cast<int>(indices[0]) - 1] -
                           matchedPoints1[static_cast<int>(indices[0]) - 1];
          modelParams[5] = matchedPoints2[(static_cast<int>(indices[0]) +
                                           matchedPoints2.size(0)) -
                                          1] -
                           matchedPoints1[(static_cast<int>(indices[0]) +
                                           matchedPoints1.size(0)) -
                                          1];
        }
      }
      for (i = 0; i < 6; i++) {
        t_idx_1 = modelParams[i];
        bv[i] = ((!std::isinf(t_idx_1)) && (!std::isnan(t_idx_1)));
      }
      if (coder::b_all(bv)) {
        t_idx_0 = evaluateModel(modelParams, matchedPoints1, matchedPoints2,
                                matchedLines1, matchedLines2, dist1);
        if (t_idx_0 < bestDis) {
          bestDis = t_idx_0;
          bestInliers.set_size(dist1.size(0));
          k = dist1.size(0);
          for (i = 0; i < k; i++) {
            bestInliers[i] = (dist1[i] < 1.5);
          }
          for (i = 0; i < 6; i++) {
            tform[i] = modelParams[i];
          }
          vlen = bestInliers.size(0);
          if (bestInliers.size(0) == 0) {
            sizes_idx_1 = 0;
          } else {
            sizes_idx_1 = bestInliers[0];
            for (k = 2; k <= vlen; k++) {
              sizes_idx_1 += bestInliers[k - 1];
            }
          }
          t_idx_0 = rt_powd_snf(
              static_cast<double>(sizes_idx_1) / static_cast<double>(u), 2.0);
          if (t_idx_0 < 2.2204460492503131E-16) {
            sizes_idx_1 = MAX_int32_T;
          } else {
            t_idx_1 =
                std::ceil(-1.9999999999999996 / std::log10(1.0 - t_idx_0));
            if (t_idx_1 < 2.147483648E+9) {
              sizes_idx_1 = static_cast<int>(t_idx_1);
            } else if (t_idx_1 >= 2.147483648E+9) {
              sizes_idx_1 = MAX_int32_T;
            } else {
              sizes_idx_1 = 0;
            }
          }
          if (numTrials > sizes_idx_1) {
            numTrials = sizes_idx_1;
          }
        }
        idxTrial++;
      } else {
        skipTrials++;
      }
    }
    vlen = bestInliers.size(0);
    if (bestInliers.size(0) == 0) {
      sizes_idx_1 = 0;
    } else {
      sizes_idx_1 = bestInliers[0];
      for (k = 2; k <= vlen; k++) {
        sizes_idx_1 += bestInliers[k - 1];
      }
    }
    if (sizes_idx_1 >= 2) {
      if (static_cast<unsigned int>(matchedPoints1.size(0)) + 1U >
          static_cast<unsigned int>(bestInliers.size(0))) {
        i = 0;
        idxTrial = 0;
        sizes_idx_1 = 0;
      } else {
        i = matchedPoints1.size(0);
        idxTrial = bestInliers.size(0);
        sizes_idx_1 = matchedPoints1.size(0);
      }
      if (matchedPoints1.size(0) < 1) {
        numTrials = 0;
      } else {
        numTrials = matchedPoints1.size(0);
      }
      k = numTrials - 1;
      vlen = 0;
      for (m = 0; m <= k; m++) {
        if (bestInliers[m]) {
          vlen++;
        }
      }
      r.set_size(vlen);
      vlen = 0;
      for (m = 0; m <= k; m++) {
        if (bestInliers[m]) {
          r[vlen] = m;
          vlen++;
        }
      }
      vlen = 0;
      for (m = 0; m <= k; m++) {
        if (bestInliers[m]) {
          vlen++;
        }
      }
      r1.set_size(vlen);
      vlen = 0;
      for (m = 0; m <= k; m++) {
        if (bestInliers[m]) {
          r1[vlen] = m;
          vlen++;
        }
      }
      k = (idxTrial - i) - 1;
      vlen = 0;
      for (m = 0; m <= k; m++) {
        if (bestInliers[i + m]) {
          vlen++;
        }
      }
      r2.set_size(vlen);
      vlen = 0;
      for (m = 0; m <= k; m++) {
        if (bestInliers[i + m]) {
          r2[vlen] = m;
          vlen++;
        }
      }
      vlen = 0;
      for (m = 0; m <= k; m++) {
        if (bestInliers[sizes_idx_1 + m]) {
          vlen++;
        }
      }
      r3.set_size(vlen);
      vlen = 0;
      for (m = 0; m <= k; m++) {
        if (bestInliers[sizes_idx_1 + m]) {
          r3[vlen] = m;
          vlen++;
        }
      }
      b_matchedPoints1.set_size(r.size(0), 2);
      b_matchedPoints2.set_size(r1.size(0), 2);
      k = r.size(0);
      vlen = r1.size(0);
      for (i = 0; i < 2; i++) {
        for (idxTrial = 0; idxTrial < k; idxTrial++) {
          b_matchedPoints1[idxTrial + b_matchedPoints1.size(0) * i] =
              matchedPoints1[r[idxTrial] + matchedPoints1.size(0) * i];
        }
        for (idxTrial = 0; idxTrial < vlen; idxTrial++) {
          b_matchedPoints2[idxTrial + b_matchedPoints2.size(0) * i] =
              matchedPoints2[r1[idxTrial] + matchedPoints2.size(0) * i];
        }
      }
      calPt2.set_size(r2.size(0), 4);
      b_matchedLines2.set_size(r3.size(0), 4);
      k = r2.size(0);
      vlen = r3.size(0);
      for (i = 0; i < 4; i++) {
        for (idxTrial = 0; idxTrial < k; idxTrial++) {
          calPt2[idxTrial + calPt2.size(0) * i] =
              matchedLines1[r2[idxTrial] + matchedLines1.size(0) * i];
        }
        for (idxTrial = 0; idxTrial < vlen; idxTrial++) {
          b_matchedLines2[idxTrial + b_matchedLines2.size(0) * i] =
              matchedLines2[r3[idxTrial] + matchedLines2.size(0) * i];
        }
      }
      svdFitForInliersPointsAndLines(b_matchedPoints1, b_matchedPoints2, calPt2,
                                     b_matchedLines2, modelParams);
      dist1.set_size(matchedPoints1.size(0));
      k = matchedPoints1.size(0);
      for (i = 0; i < k; i++) {
        dist1[i] = 0.0;
      }
      if (matchedPoints1.size(0) != 0) {
        varargin_1.set_size(2, matchedPoints1.size(0));
        k = matchedPoints1.size(0);
        for (i = 0; i < k; i++) {
          varargin_1[2 * i] = matchedPoints1[i];
          varargin_1[2 * i + 1] = matchedPoints1[i + matchedPoints1.size(0)];
        }
        b_varargin_1.set_size(3, varargin_1.size(1));
        k = varargin_1.size(1);
        for (i = 0; i < k; i++) {
          for (idxTrial = 0; idxTrial < 2; idxTrial++) {
            b_varargin_1[idxTrial + b_varargin_1.size(0) * i] =
                varargin_1[idxTrial + 2 * i];
          }
          b_varargin_1[b_varargin_1.size(0) * i + 2] = 1.0;
        }
        coder::internal::blas::mtimes(modelParams, b_varargin_1, varargin_1);
        if (matchedPoints2.size(0) == varargin_1.size(1)) {
          delta.set_size(varargin_1.size(1), 2);
          k = varargin_1.size(1);
          for (i = 0; i < 2; i++) {
            for (idxTrial = 0; idxTrial < k; idxTrial++) {
              delta[idxTrial + delta.size(0) * i] =
                  varargin_1[i + 2 * idxTrial] -
                  matchedPoints2[idxTrial + matchedPoints2.size(0) * i];
            }
          }
        } else {
          binary_expand_op(delta, varargin_1, matchedPoints2);
        }
        dist2.set_size(delta.size(0));
        k = delta.size(0);
        b_delta.set_size(delta.size(0));
        for (i = 0; i < k; i++) {
          dist2[i] = delta[i];
          b_delta[i] = delta[i + delta.size(0)];
        }
        vlen = dist2.size(0);
        sizes_idx_1 = b_delta.size(0);
        if (vlen <= sizes_idx_1) {
          sizes_idx_1 = vlen;
        }
        dist1.set_size(sizes_idx_1);
        for (k = 0; k < sizes_idx_1; k++) {
          dist1[k] = rt_hypotd_snf(dist2[k], b_delta[k]);
        }
        vlen = dist1.size(0) - 1;
        for (m = 0; m <= vlen; m++) {
          if (dist1[m] > 1.5) {
            dist1[m] = 1.5;
          }
        }
      }
      dist2.set_size(matchedLines1.size(0));
      k = matchedLines1.size(0);
      for (i = 0; i < k; i++) {
        dist2[i] = 0.0;
      }
      if (matchedLines1.size(0) != 0) {
        unsigned int varargin_2;
        signed char input_sizes_idx_0;
        varargin_2 = static_cast<unsigned int>(matchedLines1.size(0)) << 1;
        if (static_cast<int>(varargin_2) != 0) {
          sizes_idx_1 = static_cast<int>(varargin_2);
        } else {
          sizes_idx_1 = matchedLines1.size(0) << 1;
        }
        if ((sizes_idx_1 == 0) || (static_cast<int>(varargin_2) != 0)) {
          input_sizes_idx_0 = 2;
        } else {
          input_sizes_idx_0 = 0;
        }
        b_matchedLines1.set_size(4, matchedLines1.size(0));
        k = matchedLines1.size(0);
        for (i = 0; i < k; i++) {
          b_matchedLines1[4 * i] = matchedLines1[i];
          b_matchedLines1[4 * i + 1] = matchedLines1[i + matchedLines1.size(0)];
          b_matchedLines1[4 * i + 2] =
              matchedLines1[i + matchedLines1.size(0) * 2];
          b_matchedLines1[4 * i + 3] =
              matchedLines1[i + matchedLines1.size(0) * 3];
        }
        vlen = input_sizes_idx_0;
        b_varargin_1.set_size(input_sizes_idx_0 + 1, sizes_idx_1);
        for (i = 0; i < sizes_idx_1; i++) {
          for (idxTrial = 0; idxTrial < vlen; idxTrial++) {
            b_varargin_1[idxTrial + b_varargin_1.size(0) * i] =
                b_matchedLines1[idxTrial + input_sizes_idx_0 * i];
          }
          b_varargin_1[input_sizes_idx_0 + b_varargin_1.size(0) * i] = 1.0;
        }
        coder::internal::blas::mtimes(modelParams, b_varargin_1, varargin_1);
        calPt2.set_size(matchedLines1.size(0), 4);
        k = matchedLines1.size(0);
        for (i = 0; i < 4; i++) {
          for (idxTrial = 0; idxTrial < k; idxTrial++) {
            calPt2[idxTrial + calPt2.size(0) * i] =
                varargin_1[i + 4 * idxTrial];
          }
        }
        b_delta.set_size(matchedLines2.size(0));
        k = matchedLines2.size(0);
        for (i = 0; i < k; i++) {
          b_delta[i] = matchedLines2[i + matchedLines2.size(0) * 3] -
                       matchedLines2[i + matchedLines2.size(0)];
        }
        if (b_delta.size(0) == matchedLines2.size(0)) {
          D.set_size(b_delta.size(0));
          k = b_delta.size(0);
          for (i = 0; i < k; i++) {
            t_idx_0 = b_delta[i];
            t_idx_1 =
                matchedLines2[i + matchedLines2.size(0) * 2] - matchedLines2[i];
            D[i] = t_idx_0 * t_idx_0 + t_idx_1 * t_idx_1;
          }
        } else {
          binary_expand_op(D, b_delta, matchedLines2);
        }
        vlen = D.size(0);
        for (k = 0; k < vlen; k++) {
          D[k] = std::sqrt(D[k]);
        }
        r4.set_size(matchedLines2.size(0));
        k = matchedLines2.size(0);
        r5.set_size(matchedLines2.size(0));
        r6.set_size(matchedLines2.size(0));
        for (i = 0; i < k; i++) {
          t_idx_1 = matchedLines2[i];
          t_idx_0 = matchedLines2[i + matchedLines2.size(0) * 2];
          r4[i] = t_idx_1 - t_idx_0;
          r5[i] = t_idx_0 * matchedLines2[i + matchedLines2.size(0)];
          r6[i] = t_idx_1 * matchedLines2[i + matchedLines2.size(0) * 3];
        }
        if (b_delta.size(0) == 1) {
          i = calPt2.size(0);
        } else {
          i = b_delta.size(0);
        }
        if (r4.size(0) == 1) {
          idxTrial = calPt2.size(0);
        } else {
          idxTrial = r4.size(0);
        }
        if (i == 1) {
          sizes_idx_1 = idxTrial;
        } else {
          sizes_idx_1 = i;
        }
        if (sizes_idx_1 == 1) {
          numTrials = r5.size(0);
        } else {
          numTrials = sizes_idx_1;
        }
        if ((b_delta.size(0) == calPt2.size(0)) &&
            (r4.size(0) == calPt2.size(0)) && (i == idxTrial) &&
            (sizes_idx_1 == r5.size(0)) && (numTrials == r6.size(0))) {
          x.set_size(b_delta.size(0));
          k = b_delta.size(0);
          for (i = 0; i < k; i++) {
            x[i] =
                ((b_delta[i] * calPt2[i] + r4[i] * calPt2[i + calPt2.size(0)]) +
                 r5[i]) -
                r6[i];
          }
        } else {
          binary_expand_op(x, b_delta, calPt2, r4, r5, r6);
        }
        vlen = x.size(0);
        dist2.set_size(x.size(0));
        for (k = 0; k < vlen; k++) {
          dist2[k] = std::abs(x[k]);
        }
        if (b_delta.size(0) == 1) {
          i = calPt2.size(0);
        } else {
          i = b_delta.size(0);
        }
        if (r4.size(0) == 1) {
          idxTrial = calPt2.size(0);
        } else {
          idxTrial = r4.size(0);
        }
        if (i == 1) {
          sizes_idx_1 = idxTrial;
        } else {
          sizes_idx_1 = i;
        }
        if (sizes_idx_1 == 1) {
          numTrials = r5.size(0);
        } else {
          numTrials = sizes_idx_1;
        }
        if ((b_delta.size(0) == calPt2.size(0)) &&
            (r4.size(0) == calPt2.size(0)) && (i == idxTrial) &&
            (sizes_idx_1 == r5.size(0)) && (numTrials == r6.size(0))) {
          k = b_delta.size(0);
          for (i = 0; i < k; i++) {
            b_delta[i] = ((b_delta[i] * calPt2[i + calPt2.size(0) * 2] +
                           r4[i] * calPt2[i + calPt2.size(0) * 3]) +
                          r5[i]) -
                         r6[i];
          }
        } else {
          binary_expand_op(b_delta, calPt2, r4, r5, r6);
        }
        vlen = b_delta.size(0);
        x.set_size(b_delta.size(0));
        for (k = 0; k < vlen; k++) {
          x[k] = std::abs(b_delta[k]);
        }
        if (dist2.size(0) == 1) {
          i = D.size(0);
        } else {
          i = dist2.size(0);
        }
        if (x.size(0) == 1) {
          numTrials = D.size(0);
        } else {
          numTrials = x.size(0);
        }
        if ((dist2.size(0) == D.size(0)) && (x.size(0) == D.size(0)) &&
            (i == numTrials)) {
          k = dist2.size(0);
          for (i = 0; i < k; i++) {
            dist2[i] = (dist2[i] / D[i] + x[i] / D[i]) / 2.0 * 0.5;
          }
        } else {
          binary_expand_op(dist2, D, x);
        }
        vlen = dist2.size(0) - 1;
        for (m = 0; m <= vlen; m++) {
          if (dist2[m] > 1.5) {
            dist2[m] = 1.5;
          }
        }
      }
      inlierIndex.set_size(dist1.size(0) + dist2.size(0));
      k = dist1.size(0);
      for (i = 0; i < k; i++) {
        inlierIndex[i] = (dist1[i] < 1.5);
      }
      k = dist2.size(0);
      for (i = 0; i < k; i++) {
        inlierIndex[i + dist1.size(0)] = (dist2[i] < 1.5);
      }
      for (i = 0; i < 6; i++) {
        t_idx_1 = modelParams[i];
        bv[i] = ((!std::isinf(t_idx_1)) && (!std::isnan(t_idx_1)));
      }
      if ((!coder::b_all(bv)) || (!coder::any(inlierIndex))) {
        inlierIndex.set_size(static_cast<int>(u));
        for (i = 0; i < loop_ub; i++) {
          inlierIndex[i] = false;
        }
      } else {
        for (i = 0; i < 6; i++) {
          tform[i] = modelParams[i];
        }
      }
    } else {
      b_status = 2;
    }
  } else {
    b_status = 1;
  }
  *status = b_status;
}

void estgeotform2dForPtsAndLines_initialize()
{
  eml_rand_mt19937ar_stateful_init();
  isInitialized_estgeotform2dForPtsAndLines = true;
}

void estgeotform2dForPtsAndLines_terminate()
{
  isInitialized_estgeotform2dForPtsAndLines = false;
}

} // namespace estgeotform2dForPtsAndLines
