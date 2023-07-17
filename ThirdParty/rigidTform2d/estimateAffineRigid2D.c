#include "estimateAffineRigid2D.h"
#include "estimateAffineRigid2D_emxutil.h"
#include "estimateAffineRigid2D_types.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include <float.h>
#include <math.h>
#include <string.h>

static unsigned int state[625];

static boolean_T isInitialized_estimateAffineRigid2D = false;

static void b_cosd(double *x);

static double b_rand(void);

static void b_sind(double *x);

static void binary_expand_op(emxArray_real_T *in1, const emxArray_real_T *in2,
                             const emxArray_real_T *in3);

static void c_eml_rand_mt19937ar_stateful_i(void);

static void computeRigid2d(const emxArray_real_T *points, double T[9]);

static double constrainToRotationMatrix2D(const double R[4], double Rc[4]);

static void evaluateTform2d(const double tform[9],
                            const emxArray_real_T *points,
                            emxArray_real_T *dis);

static void mean(const emxArray_real_T *x, double y[2]);

static boolean_T msac(const emxArray_real_T *allPoints,
                      double bestModelParams_data[],
                      int bestModelParams_size[2], emxArray_boolean_T *inliers);

static double rt_atan2d_snf(double u0, double u1);

static double rt_hypotd_snf(double u0, double u1);

static double rt_powd_snf(double u0, double u1);

static double rt_remd_snf(double u0, double u1);

static double rt_roundd_snf(double u);

static void svd(const double A[4], double U[4], double s[2], double V[4]);

static double xnrm2(const double x[4]);

static double xrotg(double *a, double *b, double *s);

static void xswap(double x[4]);

static void b_cosd(double *x)
{
  if (rtIsInf(*x) || rtIsNaN(*x)) {
    *x = rtNaN;
  } else {
    double absx;
    signed char n;
    *x = rt_remd_snf(*x, 360.0);
    absx = fabs(*x);
    if (absx > 180.0) {
      if (*x > 0.0) {
        *x -= 360.0;
      } else {
        *x += 360.0;
      }
      absx = fabs(*x);
    }
    if (absx <= 45.0) {
      *x *= 0.017453292519943295;
      n = 0;
    } else if (absx <= 135.0) {
      if (*x > 0.0) {
        *x = 0.017453292519943295 * (*x - 90.0);
        n = 1;
      } else {
        *x = 0.017453292519943295 * (*x + 90.0);
        n = -1;
      }
    } else if (*x > 0.0) {
      *x = 0.017453292519943295 * (*x - 180.0);
      n = 2;
    } else {
      *x = 0.017453292519943295 * (*x + 180.0);
      n = -2;
    }
    if (n == 0) {
      *x = cos(*x);
    } else if (n == 1) {
      *x = -sin(*x);
    } else if (n == -1) {
      *x = sin(*x);
    } else {
      *x = -cos(*x);
    }
  }
}

static double b_rand(void)
{
  double r;
  int j;
  int kk;
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
    for (j = 0; j < 2; j++) {
      unsigned int mti;
      unsigned int y;
      mti = state[624] + 1U;
      if (state[624] + 1U >= 625U) {
        for (kk = 0; kk < 227; kk++) {
          y = (state[kk] & 2147483648U) | (state[kk + 1] & 2147483647U);
          if ((y & 1U) == 0U) {
            y >>= 1U;
          } else {
            y = y >> 1U ^ 2567483615U;
          }
          state[kk] = state[kk + 397] ^ y;
        }
        for (kk = 0; kk < 396; kk++) {
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
      y = state[(int)mti - 1];
      state[624] = mti;
      y ^= y >> 11U;
      y ^= y << 7U & 2636928640U;
      y ^= y << 15U & 4022730752U;
      u[j] = y ^ y >> 18U;
    }
    u[0] >>= 5U;
    u[1] >>= 6U;
    r = 1.1102230246251565E-16 * ((double)u[0] * 6.7108864E+7 + (double)u[1]);
  } while (r == 0.0);
  return r;
}

static void b_sind(double *x)
{
  if (rtIsInf(*x) || rtIsNaN(*x)) {
    *x = rtNaN;
  } else {
    double absx;
    signed char n;
    *x = rt_remd_snf(*x, 360.0);
    absx = fabs(*x);
    if (absx > 180.0) {
      if (*x > 0.0) {
        *x -= 360.0;
      } else {
        *x += 360.0;
      }
      absx = fabs(*x);
    }
    if (absx <= 45.0) {
      *x *= 0.017453292519943295;
      n = 0;
    } else if (absx <= 135.0) {
      if (*x > 0.0) {
        *x = 0.017453292519943295 * (*x - 90.0);
        n = 1;
      } else {
        *x = 0.017453292519943295 * (*x + 90.0);
        n = -1;
      }
    } else if (*x > 0.0) {
      *x = 0.017453292519943295 * (*x - 180.0);
      n = 2;
    } else {
      *x = 0.017453292519943295 * (*x + 180.0);
      n = -2;
    }
    if (n == 0) {
      *x = sin(*x);
    } else if (n == 1) {
      *x = cos(*x);
    } else if (n == -1) {
      *x = -cos(*x);
    } else {
      *x = -sin(*x);
    }
  }
}

static void binary_expand_op(emxArray_real_T *in1, const emxArray_real_T *in2,
                             const emxArray_real_T *in3)
{
  emxArray_real_T *b_in1;
  const double *in2_data;
  const double *in3_data;
  double *b_in1_data;
  double *in1_data;
  int i;
  int i1;
  int in2_idx_0;
  int stride_0_0;
  int stride_1_0;
  in3_data = in3->data;
  in2_data = in2->data;
  in2_idx_0 = in2->size[0];
  i = in1->size[0] * in1->size[1];
  in1->size[0] = in2_idx_0;
  in1->size[1] = 2;
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  for (i = 0; i < in2_idx_0; i++) {
    double d;
    d = in2_data[i + in2->size[0] * 2];
    in1_data[i] = in2_data[i] / d;
    in1_data[i + in1->size[0]] = in2_data[i + in2->size[0]] / d;
  }
  emxInit_real_T(&b_in1, 2);
  if (in3->size[0] == 1) {
    in2_idx_0 = in1->size[0];
  } else {
    in2_idx_0 = in3->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = in2_idx_0;
  b_in1->size[1] = 2;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_1_0 = (in3->size[0] != 1);
  for (i = 0; i < 2; i++) {
    for (i1 = 0; i1 < in2_idx_0; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[i1 * stride_0_0 + in1->size[0] * i] -
          in3_data[(i1 * stride_1_0 + in3->size[0] * i) + in3->size[0] * 2];
    }
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = b_in1->size[0];
  in1->size[1] = 2;
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  in2_idx_0 = b_in1->size[0];
  for (i = 0; i < 2; i++) {
    for (i1 = 0; i1 < in2_idx_0; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real_T(&b_in1);
}

static void c_eml_rand_mt19937ar_stateful_i(void)
{
  int mti;
  unsigned int r;
  memset(&state[0], 0, 625U * sizeof(unsigned int));
  r = 5489U;
  state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = ((r ^ r >> 30U) * 1812433253U + (unsigned int)mti) + 1U;
    state[mti + 1] = r;
  }
  state[624] = 624U;
}

static void computeRigid2d(const emxArray_real_T *points, double T[9])
{
  emxArray_real_T *normPoints1;
  emxArray_real_T *normPoints2;
  double C[4];
  double U[4];
  double V[4];
  double centroid1[2];
  double centroid2[2];
  const double *points_data;
  double bkj;
  double d;
  double d1;
  double d2;
  double d3;
  double d4;
  double *normPoints1_data;
  double *normPoints2_data;
  int acoef;
  int boffset;
  int coffset;
  int j;
  int k;
  signed char csz[2];
  boolean_T p;
  points_data = points->data;
  emxInit_real_T(&normPoints1, 2);
  boffset = normPoints1->size[0] * normPoints1->size[1];
  normPoints1->size[0] = points->size[0];
  normPoints1->size[1] = 2;
  emxEnsureCapacity_real_T(normPoints1, boffset);
  normPoints1_data = normPoints1->data;
  acoef = points->size[0];
  for (boffset = 0; boffset < 2; boffset++) {
    for (coffset = 0; coffset < acoef; coffset++) {
      normPoints1_data[coffset + normPoints1->size[0] * boffset] =
          points_data[coffset + points->size[0] * boffset];
    }
  }
  mean(normPoints1, centroid1);
  boffset = normPoints1->size[0] * normPoints1->size[1];
  normPoints1->size[0] = points->size[0];
  normPoints1->size[1] = 2;
  emxEnsureCapacity_real_T(normPoints1, boffset);
  normPoints1_data = normPoints1->data;
  acoef = points->size[0];
  for (boffset = 0; boffset < 2; boffset++) {
    for (coffset = 0; coffset < acoef; coffset++) {
      normPoints1_data[coffset + normPoints1->size[0] * boffset] =
          points_data[(coffset + points->size[0] * boffset) +
                      points->size[0] * 2];
    }
  }
  mean(normPoints1, centroid2);
  boffset = normPoints1->size[0] * normPoints1->size[1];
  normPoints1->size[0] = points->size[0];
  normPoints1->size[1] = 2;
  emxEnsureCapacity_real_T(normPoints1, boffset);
  normPoints1_data = normPoints1->data;
  if (points->size[0] != 0) {
    acoef = (points->size[0] != 1);
    for (k = 0; k < 2; k++) {
      boffset = normPoints1->size[0] - 1;
      for (coffset = 0; coffset <= boffset; coffset++) {
        normPoints1_data[coffset + normPoints1->size[0] * k] =
            points_data[acoef * coffset + points->size[0] * k] - centroid1[k];
      }
    }
  }
  emxInit_real_T(&normPoints2, 2);
  boffset = normPoints2->size[0] * normPoints2->size[1];
  normPoints2->size[0] = points->size[0];
  normPoints2->size[1] = 2;
  emxEnsureCapacity_real_T(normPoints2, boffset);
  normPoints2_data = normPoints2->data;
  if (points->size[0] != 0) {
    acoef = (points->size[0] != 1);
    for (k = 0; k < 2; k++) {
      boffset = normPoints2->size[0] - 1;
      for (coffset = 0; coffset <= boffset; coffset++) {
        normPoints2_data[coffset + normPoints2->size[0] * k] =
            points_data[(acoef * coffset + points->size[0] * k) +
                        points->size[0] * 2] -
            centroid2[k];
      }
    }
  }
  acoef = normPoints1->size[0];
  for (j = 0; j < 2; j++) {
    coffset = j << 1;
    boffset = j * normPoints2->size[0];
    C[coffset] = 0.0;
    C[coffset + 1] = 0.0;
    for (k = 0; k < acoef; k++) {
      bkj = normPoints2_data[boffset + k];
      C[coffset] += normPoints1_data[k] * bkj;
      C[coffset + 1] += normPoints1_data[normPoints1->size[0] + k] * bkj;
    }
  }
  emxFree_real_T(&normPoints2);
  emxFree_real_T(&normPoints1);
  p = true;
  if (rtIsInf(C[0]) || rtIsNaN(C[0]) || (rtIsInf(C[1]) || rtIsNaN(C[1]))) {
    p = false;
  }
  if ((!p) || (rtIsInf(C[2]) || rtIsNaN(C[2]))) {
    p = false;
  }
  if ((!p) || (rtIsInf(C[3]) || rtIsNaN(C[3]))) {
    p = false;
  }
  if (p) {
    double s[2];
    svd(C, U, s, V);
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
  bkj = V[0];
  d = V[2];
  d1 = V[1];
  d2 = V[3];
  for (boffset = 0; boffset < 2; boffset++) {
    d3 = U[boffset + 2];
    d4 = U[boffset];
    C[boffset] = d4 * bkj + d3 * d;
    C[boffset + 2] = d4 * d1 + d3 * d2;
    csz[boffset] = (signed char)(boffset + 1);
  }
  acoef = 0;
  if (fabs(C[1]) > fabs(C[0])) {
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
  if (rtIsNaN(bkj)) {
    C[3] = rtNaN;
  } else if (bkj < 0.0) {
    C[3] = -1.0;
  } else {
    C[3] = (bkj > 0.0);
  }
  bkj = U[0];
  d = U[2];
  d1 = U[1];
  d2 = U[3];
  for (boffset = 0; boffset < 2; boffset++) {
    double d5;
    d3 = V[boffset + 2];
    d4 = V[boffset];
    d5 = d4 + d3 * 0.0;
    d3 = d4 * 0.0 + d3 * C[3];
    C[boffset] = d5 * bkj + d3 * d;
    C[boffset + 2] = d5 * d1 + d3 * d2;
  }
  memset(&T[0], 0, 9U * sizeof(double));
  T[8] = 1.0;
  T[0] = C[0];
  T[1] = C[2];
  T[2] = centroid2[0] - (C[0] * centroid1[0] + centroid1[1] * C[2]);
  T[3] = C[1];
  T[4] = C[3];
  T[5] = centroid2[1] - (centroid1[0] * C[1] + centroid1[1] * C[3]);
}

static double constrainToRotationMatrix2D(const double R[4], double Rc[4])
{
  double e[2];
  double s[2];
  double R_clamped_idx_0;
  double R_clamped_idx_1;
  double R_clamped_idx_2;
  double R_clamped_idx_3;
  double b;
  double d;
  double f;
  double minval_idx_2;
  double r;
  double sm;
  double snorm;
  double sqds;
  double wpr;
  int kase;
  boolean_T close_enough;
  R_clamped_idx_0 = fmax(fmin(R[0], 1.0), -1.0);
  R_clamped_idx_1 = fmax(fmin(R[1], 1.0), -1.0);
  R_clamped_idx_2 = fmax(fmin(R[2], 1.0), -1.0);
  R_clamped_idx_3 = fmax(fmin(R[3], 1.0), -1.0);
  wpr = 57.295779513082323 * rt_atan2d_snf(R_clamped_idx_1, R_clamped_idx_3) +
        180.0;
  if (rtIsNaN(wpr) || rtIsInf(wpr)) {
    r = rtNaN;
  } else if (wpr == 0.0) {
    r = 0.0;
  } else {
    r = fmod(wpr, 360.0);
    if (r == 0.0) {
      r = 0.0;
    } else if (wpr < 0.0) {
      r += 360.0;
    }
  }
  wpr = rt_roundd_snf(r - 180.0);
  if (r - 180.0 == wpr) {
    close_enough = true;
  } else {
    d = fabs((r - 180.0) - wpr);
    if ((r - 180.0 == 0.0) || (wpr == 0.0)) {
      close_enough = (d < 4.94065645841247E-324);
    } else {
      b = fabs(r - 180.0) + fabs(wpr);
      if (b < 2.2250738585072014E-308) {
        close_enough = (d < 4.94065645841247E-324);
      } else {
        close_enough =
            (d / fmin(b, 1.7976931348623157E+308) < 2.2204460492503131E-16);
      }
    }
  }
  r -= 180.0;
  if (close_enough) {
    r = wpr;
  }
  wpr = r;
  b_sind(&wpr);
  d = r;
  b_cosd(&d);
  Rc[0] = d;
  Rc[2] = -wpr;
  Rc[1] = wpr;
  Rc[3] = d;
  f = R_clamped_idx_0 - d;
  sqds = R_clamped_idx_1 - wpr;
  minval_idx_2 = R_clamped_idx_2 - (-wpr);
  snorm = R_clamped_idx_3 - d;
  wpr = 0.0;
  b = fabs(f);
  if (rtIsNaN(b) || (b > 0.0)) {
    wpr = b;
  }
  sm = fabs(sqds);
  if (rtIsNaN(sm) || (sm > wpr)) {
    wpr = sm;
  }
  d = fabs(minval_idx_2);
  if (rtIsNaN(d) || (d > wpr)) {
    wpr = d;
  }
  d = fabs(snorm);
  if (rtIsNaN(d) || (d > wpr)) {
    wpr = d;
  }
  if ((!rtIsInf(wpr)) && (!rtIsNaN(wpr))) {
    double scale;
    int iter;
    int m;
    scale = 3.3121686421112381E-170;
    if (b > 3.3121686421112381E-170) {
      d = 1.0;
      scale = b;
    } else {
      wpr = b / 3.3121686421112381E-170;
      d = wpr * wpr;
    }
    if (sm > scale) {
      wpr = scale / sm;
      d = d * wpr * wpr + 1.0;
      scale = sm;
    } else {
      wpr = sm / scale;
      d += wpr * wpr;
    }
    d = scale * sqrt(d);
    if (d > 0.0) {
      if (f < 0.0) {
        s[0] = -d;
      } else {
        s[0] = d;
      }
      if (fabs(s[0]) >= 1.0020841800044864E-292) {
        wpr = 1.0 / s[0];
        f *= wpr;
        sqds *= wpr;
      } else {
        f /= s[0];
        sqds /= s[0];
      }
      f++;
      s[0] = -s[0];
      wpr = -((f * minval_idx_2 + sqds * snorm) / f);
      if (!(wpr == 0.0)) {
        minval_idx_2 += wpr * f;
        snorm += wpr * sqds;
      }
    } else {
      s[0] = 0.0;
    }
    m = 2;
    s[1] = snorm;
    e[0] = minval_idx_2;
    e[1] = 0.0;
    if (s[0] != 0.0) {
      d = fabs(s[0]);
      wpr = s[0] / d;
      s[0] = d;
      e[0] = minval_idx_2 / wpr;
    }
    if (e[0] != 0.0) {
      d = fabs(e[0]);
      wpr = d / e[0];
      e[0] = d;
      s[1] = snorm * wpr;
    }
    if (s[1] != 0.0) {
      s[1] = fabs(s[1]);
    }
    iter = 0;
    snorm = fmax(fmax(s[0], e[0]), fmax(s[1], 0.0));
    while ((m > 0) && (iter < 75)) {
      int ii_tmp_tmp;
      int q;
      boolean_T exitg1;
      ii_tmp_tmp = m - 1;
      q = m - 1;
      exitg1 = false;
      while (!(exitg1 || (q == 0))) {
        wpr = fabs(e[0]);
        if ((wpr <= 2.2204460492503131E-16 * (fabs(s[0]) + fabs(s[1]))) ||
            (wpr <= 1.0020841800044864E-292) ||
            ((iter > 20) && (wpr <= 2.2204460492503131E-16 * snorm))) {
          e[0] = 0.0;
          exitg1 = true;
        } else {
          q = 0;
        }
      }
      if (q == m - 1) {
        kase = 4;
      } else {
        int qs;
        qs = m;
        kase = m;
        exitg1 = false;
        while ((!exitg1) && (kase >= q)) {
          qs = kase;
          if (kase == q) {
            exitg1 = true;
          } else {
            wpr = 0.0;
            if (kase < m) {
              wpr = fabs(e[0]);
            }
            if (kase > q + 1) {
              wpr += fabs(e[0]);
            }
            d = fabs(s[kase - 1]);
            if ((d <= 2.2204460492503131E-16 * wpr) ||
                (d <= 1.0020841800044864E-292)) {
              s[kase - 1] = 0.0;
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
      case 1:
        f = e[0];
        e[0] = 0.0;
        for (kase = ii_tmp_tmp; kase >= q + 1; kase--) {
          xrotg(&s[0], &f, &sm);
        }
        break;
      case 2:
        f = e[q - 1];
        e[q - 1] = 0.0;
        for (kase = q + 1; kase <= m; kase++) {
          d = xrotg(&s[kase - 1], &f, &sm);
          wpr = e[kase - 1];
          f = -sm * wpr;
          e[kase - 1] = wpr * d;
        }
        break;
      case 3:
        wpr = s[m - 1];
        scale = fmax(
            fmax(fmax(fmax(fabs(wpr), fabs(s[0])), fabs(e[0])), fabs(s[q])),
            fabs(e[q]));
        sm = wpr / scale;
        wpr = s[0] / scale;
        d = e[0] / scale;
        sqds = s[q] / scale;
        b = ((wpr + sm) * (wpr - sm) + d * d) / 2.0;
        wpr = sm * d;
        wpr *= wpr;
        if ((b != 0.0) || (wpr != 0.0)) {
          d = sqrt(b * b + wpr);
          if (b < 0.0) {
            d = -d;
          }
          d = wpr / (b + d);
        } else {
          d = 0.0;
        }
        f = (sqds + sm) * (sqds - sm) + d;
        wpr = sqds * (e[q] / scale);
        for (kase = q + 1; kase < 2; kase++) {
          d = xrotg(&f, &wpr, &sm);
          f = d * s[0] + sm * e[0];
          b = d * e[0] - sm * s[0];
          e[0] = b;
          wpr = sm * s[1];
          s[1] *= d;
          s[0] = f;
          d = xrotg(&s[0], &wpr, &sm);
          f = d * b + sm * s[1];
          s[1] = -sm * b + d * s[1];
          wpr = sm * e[1];
          e[1] *= d;
        }
        e[0] = f;
        iter++;
        break;
      default:
        if (s[q] < 0.0) {
          s[q] = -s[q];
        }
        while ((q + 1 < 2) && (s[0] < s[1])) {
          d = s[0];
          s[0] = s[1];
          s[1] = d;
          q = 1;
        }
        iter = 0;
        m--;
        break;
      }
    }
    wpr = s[0];
  }
  if (wpr / 2.2204460492503131E-16 < 10.0) {
    Rc[0] = R_clamped_idx_0;
    Rc[1] = R_clamped_idx_1;
    Rc[2] = R_clamped_idx_2;
    Rc[3] = R_clamped_idx_3;
  }
  return r;
}

static void evaluateTform2d(const double tform[9],
                            const emxArray_real_T *points, emxArray_real_T *dis)
{
  emxArray_real_T *b_points;
  emxArray_real_T *b_pt1h;
  emxArray_real_T *delta;
  emxArray_real_T *pt1h;
  emxArray_real_T *y;
  const double *points_data;
  double bkj;
  double *b_pt1h_data;
  double *delta_data;
  double *pt1h_data;
  int b_i;
  int coffset;
  int i;
  int j;
  int k;
  int loop_ub;
  int nx;
  signed char input_sizes_idx_1;
  signed char sizes_idx_1;
  boolean_T empty_non_axis_sizes;
  points_data = points->data;
  if (points->size[0] != 0) {
    nx = points->size[0];
  } else {
    nx = 0;
  }
  empty_non_axis_sizes = (nx == 0);
  if (empty_non_axis_sizes || (points->size[0] != 0)) {
    input_sizes_idx_1 = 2;
  } else {
    input_sizes_idx_1 = 0;
  }
  if (empty_non_axis_sizes || (points->size[0] != 0)) {
    sizes_idx_1 = 1;
  } else {
    sizes_idx_1 = 0;
  }
  emxInit_real_T(&b_points, 2);
  i = b_points->size[0] * b_points->size[1];
  b_points->size[0] = points->size[0];
  b_points->size[1] = 2;
  emxEnsureCapacity_real_T(b_points, i);
  delta_data = b_points->data;
  loop_ub = points->size[0];
  for (i = 0; i < 2; i++) {
    for (coffset = 0; coffset < loop_ub; coffset++) {
      delta_data[coffset + b_points->size[0] * i] =
          points_data[coffset + points->size[0] * i];
    }
  }
  emxInit_real_T(&pt1h, 2);
  i = pt1h->size[0] * pt1h->size[1];
  pt1h->size[0] = nx;
  pt1h->size[1] = input_sizes_idx_1 + sizes_idx_1;
  emxEnsureCapacity_real_T(pt1h, i);
  pt1h_data = pt1h->data;
  loop_ub = input_sizes_idx_1;
  for (i = 0; i < loop_ub; i++) {
    for (coffset = 0; coffset < nx; coffset++) {
      pt1h_data[coffset + pt1h->size[0] * i] = delta_data[coffset + nx * i];
    }
  }
  emxFree_real_T(&b_points);
  loop_ub = sizes_idx_1;
  for (i = 0; i < loop_ub; i++) {
    for (coffset = 0; coffset < nx; coffset++) {
      pt1h_data[coffset + pt1h->size[0] * input_sizes_idx_1] = 1.0;
    }
  }
  nx = pt1h->size[0];
  loop_ub = pt1h->size[1];
  emxInit_real_T(&b_pt1h, 2);
  i = b_pt1h->size[0] * b_pt1h->size[1];
  b_pt1h->size[0] = pt1h->size[0];
  b_pt1h->size[1] = 3;
  emxEnsureCapacity_real_T(b_pt1h, i);
  b_pt1h_data = b_pt1h->data;
  for (j = 0; j < 3; j++) {
    int boffset;
    coffset = j * nx;
    boffset = j * 3;
    for (b_i = 0; b_i < nx; b_i++) {
      b_pt1h_data[coffset + b_i] = 0.0;
    }
    for (k = 0; k < loop_ub; k++) {
      int aoffset;
      aoffset = k * pt1h->size[0];
      bkj = tform[boffset + k];
      for (b_i = 0; b_i < nx; b_i++) {
        i = coffset + b_i;
        b_pt1h_data[i] += pt1h_data[aoffset + b_i] * bkj;
      }
    }
  }
  emxFree_real_T(&pt1h);
  emxInit_real_T(&delta, 2);
  if (b_pt1h->size[0] == points->size[0]) {
    i = delta->size[0] * delta->size[1];
    delta->size[0] = b_pt1h->size[0];
    delta->size[1] = 2;
    emxEnsureCapacity_real_T(delta, i);
    delta_data = delta->data;
    loop_ub = b_pt1h->size[0];
    for (i = 0; i < loop_ub; i++) {
      bkj = b_pt1h_data[i + b_pt1h->size[0] * 2];
      delta_data[i] = b_pt1h_data[i] / bkj;
      delta_data[i + delta->size[0]] = b_pt1h_data[i + b_pt1h->size[0]] / bkj;
    }
    i = delta->size[0] * delta->size[1];
    delta->size[1] = 2;
    emxEnsureCapacity_real_T(delta, i);
    delta_data = delta->data;
    for (i = 0; i < 2; i++) {
      loop_ub = delta->size[0];
      for (coffset = 0; coffset < loop_ub; coffset++) {
        delta_data[coffset + delta->size[0] * i] -=
            points_data[(coffset + points->size[0] * i) + points->size[0] * 2];
      }
    }
  } else {
    binary_expand_op(delta, b_pt1h, points);
    delta_data = delta->data;
  }
  i = dis->size[0];
  dis->size[0] = delta->size[0];
  emxEnsureCapacity_real_T(dis, i);
  pt1h_data = dis->data;
  nx = delta->size[0];
  for (k = 0; k < nx; k++) {
    pt1h_data[k] = rt_hypotd_snf(delta_data[k], delta_data[k + delta->size[0]]);
  }
  emxFree_real_T(&delta);
  nx = b_pt1h->size[0];
  emxInit_real_T(&y, 1);
  i = y->size[0];
  y->size[0] = b_pt1h->size[0];
  emxEnsureCapacity_real_T(y, i);
  delta_data = y->data;
  for (k = 0; k < nx; k++) {
    delta_data[k] = fabs(b_pt1h_data[k + b_pt1h->size[0] * 2]);
  }
  emxFree_real_T(&b_pt1h);
  nx = y->size[0] - 1;
  for (b_i = 0; b_i <= nx; b_i++) {
    if (delta_data[b_i] < 2.2204460492503131E-16) {
      pt1h_data[b_i] = rtInf;
    }
  }
  emxFree_real_T(&y);
}

static void mean(const emxArray_real_T *x, double y[2])
{
  const double *x_data;
  int ib;
  int k;
  int xi;
  x_data = x->data;
  if (x->size[0] == 0) {
    y[0] = 0.0;
    y[1] = 0.0;
  } else {
    int firstBlockLength;
    int lastBlockLength;
    int nblocks;
    if (x->size[0] <= 1024) {
      firstBlockLength = x->size[0];
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = (int)((unsigned int)x->size[0] >> 10);
      lastBlockLength = x->size[0] - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    for (xi = 0; xi < 2; xi++) {
      int xpageoffset;
      xpageoffset = xi * x->size[0];
      y[xi] = x_data[xpageoffset];
      for (k = 2; k <= firstBlockLength; k++) {
        y[xi] += x_data[(xpageoffset + k) - 1];
      }
      for (ib = 2; ib <= nblocks; ib++) {
        double bsum;
        int hi;
        int xblockoffset;
        xblockoffset = xpageoffset + ((ib - 1) << 10);
        bsum = x_data[xblockoffset];
        if (ib == nblocks) {
          hi = lastBlockLength;
        } else {
          hi = 1024;
        }
        for (k = 2; k <= hi; k++) {
          bsum += x_data[(xblockoffset + k) - 1];
        }
        y[xi] += bsum;
      }
    }
  }
  y[0] /= (double)x->size[0];
  y[1] /= (double)x->size[0];
}

static boolean_T msac(const emxArray_real_T *allPoints,
                      double bestModelParams_data[],
                      int bestModelParams_size[2], emxArray_boolean_T *inliers)
{
  emxArray_boolean_T *bestInliers;
  emxArray_int32_T *r;
  emxArray_real_T b_samplePoints_data;
  emxArray_real_T *b_allPoints;
  emxArray_real_T *dis;
  double modelParams[9];
  double samplePoints_data[8];
  const double *allPoints_data;
  double bestDis;
  double j;
  double *dis_data;
  int samplePoints_size[3];
  int ib;
  int idxTrial;
  int jlast;
  int k;
  int lastBlockLength;
  int nblocks;
  int numPts;
  int numTrials;
  int nz;
  int skipTrials;
  int *r1;
  boolean_T b_data[9];
  boolean_T exitg1;
  boolean_T isFound;
  boolean_T isValidModel;
  boolean_T *bestInliers_data;
  allPoints_data = allPoints->data;
  numPts = allPoints->size[0];
  idxTrial = 1;
  numTrials = 1000;
  bestDis = 1.5 * (double)allPoints->size[0];
  bestModelParams_size[0] = 0;
  bestModelParams_size[1] = 0;
  skipTrials = 0;
  emxInit_boolean_T(&bestInliers, 1);
  lastBlockLength = bestInliers->size[0];
  bestInliers->size[0] = allPoints->size[0];
  emxEnsureCapacity_boolean_T(bestInliers, lastBlockLength);
  bestInliers_data = bestInliers->data;
  nz = allPoints->size[0];
  for (lastBlockLength = 0; lastBlockLength < nz; lastBlockLength++) {
    bestInliers_data[lastBlockLength] = false;
  }
  emxInit_real_T(&dis, 1);
  while ((idxTrial <= numTrials) && (skipTrials < 10000)) {
    double indices_data[2];
    double selectedLoc;
    boolean_T b[9];
    indices_data[1] = 0.0;
    if (numPts <= 2) {
      indices_data[0] = 1.0;
      j = b_rand() * 2.0;
      j = floor(j);
      indices_data[1] = indices_data[(int)(j + 1.0) - 1];
      indices_data[(int)(j + 1.0) - 1] = 2.0;
    } else if ((double)numPts / 4.0 <= 2.0) {
      double loc_data_idx_0;
      double t;
      t = 0.0;
      selectedLoc = numPts;
      loc_data_idx_0 = 2.0 / (double)numPts;
      j = b_rand();
      while (j > loc_data_idx_0) {
        t++;
        selectedLoc--;
        loc_data_idx_0 += (1.0 - loc_data_idx_0) * (2.0 / selectedLoc);
      }
      t++;
      j = b_rand();
      j = floor(j);
      indices_data[0] = 0.0;
      indices_data[(int)(j + 1.0) - 1] = t;
      selectedLoc = (double)numPts - t;
      loc_data_idx_0 = 1.0 / selectedLoc;
      j = b_rand();
      while (j > loc_data_idx_0) {
        t++;
        selectedLoc--;
        loc_data_idx_0 += (1.0 - loc_data_idx_0) * (1.0 / selectedLoc);
      }
      t++;
      j = b_rand() * 2.0;
      j = floor(j);
      indices_data[1] = indices_data[(int)(j + 1.0) - 1];
      indices_data[(int)(j + 1.0) - 1] = t;
    } else {
      double loc_data_idx_0;
      signed char hashTbl_data[2];
      hashTbl_data[0] = 0;
      hashTbl_data[1] = 0;
      selectedLoc = b_rand() * (((double)numPts - 1.0) + 1.0);
      selectedLoc = floor(selectedLoc);
      if (selectedLoc == 0.0) {
        j = 0.0;
      } else {
        j = fmod(selectedLoc, 2.0);
        if (j == 0.0) {
          j = 0.0;
        }
      }
      indices_data[0] = selectedLoc + 1.0;
      loc_data_idx_0 = selectedLoc;
      hashTbl_data[(int)(j + 1.0) - 1] = 1;
      jlast = hashTbl_data[(int)fmod((double)numPts - 1.0, 2.0)];
      while ((jlast > 0) && (selectedLoc != (double)numPts - 1.0)) {
        jlast = 0;
      }
      if (jlast > 0) {
        jlast = 0;
      } else {
        jlast = numPts - 1;
      }
      selectedLoc = b_rand() * (((double)numPts - 2.0) + 1.0);
      selectedLoc = floor(selectedLoc);
      if (selectedLoc == 0.0) {
        j = 0.0;
      } else {
        j = fmod(selectedLoc, 2.0);
        if (j == 0.0) {
          j = 0.0;
        }
      }
      j = hashTbl_data[(int)(j + 1.0) - 1];
      while ((j > 0.0) && (loc_data_idx_0 != selectedLoc)) {
        j = 0.0;
      }
      if (j > 0.0) {
        indices_data[1] = (double)jlast + 1.0;
      } else {
        indices_data[1] = selectedLoc + 1.0;
      }
    }
    samplePoints_size[0] = 2;
    samplePoints_size[1] = 2;
    samplePoints_size[2] = 2;
    j = indices_data[0];
    selectedLoc = indices_data[1];
    for (lastBlockLength = 0; lastBlockLength < 2; lastBlockLength++) {
      samplePoints_data[4 * lastBlockLength] =
          allPoints_data[((int)j + allPoints->size[0] * 2 * lastBlockLength) -
                         1];
      samplePoints_data[4 * lastBlockLength + 1] = allPoints_data
          [((int)selectedLoc + allPoints->size[0] * 2 * lastBlockLength) - 1];
      samplePoints_data[4 * lastBlockLength + 2] =
          allPoints_data[(((int)j + allPoints->size[0]) +
                          allPoints->size[0] * 2 * lastBlockLength) -
                         1];
      samplePoints_data[4 * lastBlockLength + 3] =
          allPoints_data[(((int)selectedLoc + allPoints->size[0]) +
                          allPoints->size[0] * 2 * lastBlockLength) -
                         1];
    }
    b_samplePoints_data.data = &samplePoints_data[0];
    b_samplePoints_data.size = &samplePoints_size[0];
    b_samplePoints_data.allocatedSize = 8;
    b_samplePoints_data.numDimensions = 3;
    b_samplePoints_data.canFreeData = false;
    computeRigid2d(&b_samplePoints_data, modelParams);
    for (lastBlockLength = 0; lastBlockLength < 9; lastBlockLength++) {
      j = modelParams[lastBlockLength];
      b_data[lastBlockLength] = rtIsInf(j);
      b[lastBlockLength] = rtIsNaN(j);
    }
    isValidModel = true;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 9)) {
      if (b_data[k] || b[k]) {
        isValidModel = false;
        exitg1 = true;
      } else {
        k++;
      }
    }
    if (isValidModel) {
      evaluateTform2d(modelParams, allPoints, dis);
      dis_data = dis->data;
      nz = dis->size[0] - 1;
      for (nblocks = 0; nblocks <= nz; nblocks++) {
        if (dis_data[nblocks] > 1.5) {
          dis_data[nblocks] = 1.5;
        }
      }
      if (dis->size[0] == 0) {
        j = 0.0;
      } else {
        if (dis->size[0] <= 1024) {
          jlast = dis->size[0];
          lastBlockLength = 0;
          nblocks = 1;
        } else {
          jlast = 1024;
          nblocks = (int)((unsigned int)dis->size[0] >> 10);
          lastBlockLength = dis->size[0] - (nblocks << 10);
          if (lastBlockLength > 0) {
            nblocks++;
          } else {
            lastBlockLength = 1024;
          }
        }
        j = dis_data[0];
        for (k = 2; k <= jlast; k++) {
          j += dis_data[k - 1];
        }
        for (ib = 2; ib <= nblocks; ib++) {
          jlast = (ib - 1) << 10;
          selectedLoc = dis_data[jlast];
          if (ib == nblocks) {
            nz = lastBlockLength;
          } else {
            nz = 1024;
          }
          for (k = 2; k <= nz; k++) {
            selectedLoc += dis_data[(jlast + k) - 1];
          }
          j += selectedLoc;
        }
      }
      if (j < bestDis) {
        bestDis = j;
        lastBlockLength = bestInliers->size[0];
        bestInliers->size[0] = dis->size[0];
        emxEnsureCapacity_boolean_T(bestInliers, lastBlockLength);
        bestInliers_data = bestInliers->data;
        nz = dis->size[0];
        for (lastBlockLength = 0; lastBlockLength < nz; lastBlockLength++) {
          bestInliers_data[lastBlockLength] = (dis_data[lastBlockLength] < 1.5);
        }
        bestModelParams_size[0] = 3;
        bestModelParams_size[1] = 3;
        memcpy(&bestModelParams_data[0], &modelParams[0], 9U * sizeof(double));
        jlast = bestInliers->size[0];
        if (bestInliers->size[0] == 0) {
          nz = 0;
        } else {
          nz = bestInliers_data[0];
          for (k = 2; k <= jlast; k++) {
            nz += bestInliers_data[k - 1];
          }
        }
        j = rt_powd_snf((double)nz / (double)numPts, 2.0);
        if (j < 2.2204460492503131E-16) {
          jlast = MAX_int32_T;
        } else {
          j = ceil(-1.9999999999999996 / log10(1.0 - j));
          if (j < 2.147483648E+9) {
            jlast = (int)j;
          } else if (j >= 2.147483648E+9) {
            jlast = MAX_int32_T;
          } else {
            jlast = 0;
          }
        }
        if (numTrials > jlast) {
          numTrials = jlast;
        }
      }
      idxTrial++;
    } else {
      skipTrials++;
    }
  }
  nz = bestModelParams_size[0] * bestModelParams_size[1];
  for (lastBlockLength = 0; lastBlockLength < nz; lastBlockLength++) {
    j = bestModelParams_data[lastBlockLength];
    b_data[lastBlockLength] = ((!rtIsInf(j)) && (!rtIsNaN(j)));
  }
  isValidModel = true;
  jlast = 1;
  exitg1 = false;
  while ((!exitg1) && (jlast <= nz)) {
    if (!b_data[jlast - 1]) {
      isValidModel = false;
      exitg1 = true;
    } else {
      jlast++;
    }
  }
  if (isValidModel && (bestInliers->size[0] != 0)) {
    jlast = bestInliers->size[0];
    nz = bestInliers_data[0];
    for (k = 2; k <= jlast; k++) {
      nz += bestInliers_data[k - 1];
    }
    if (nz >= 2) {
      isFound = true;
    } else {
      isFound = false;
    }
  } else {
    isFound = false;
  }
  emxInit_int32_T(&r);
  emxInit_real_T(&b_allPoints, 3);
  if (isFound) {
    boolean_T guard1 = false;
    nz = bestInliers->size[0] - 1;
    jlast = 0;
    for (nblocks = 0; nblocks <= nz; nblocks++) {
      if (bestInliers_data[nblocks]) {
        jlast++;
      }
    }
    lastBlockLength = r->size[0];
    r->size[0] = jlast;
    emxEnsureCapacity_int32_T(r, lastBlockLength);
    r1 = r->data;
    jlast = 0;
    for (nblocks = 0; nblocks <= nz; nblocks++) {
      if (bestInliers_data[nblocks]) {
        r1[jlast] = nblocks;
        jlast++;
      }
    }
    lastBlockLength =
        b_allPoints->size[0] * b_allPoints->size[1] * b_allPoints->size[2];
    b_allPoints->size[0] = r->size[0];
    b_allPoints->size[1] = 2;
    b_allPoints->size[2] = 2;
    emxEnsureCapacity_real_T(b_allPoints, lastBlockLength);
    dis_data = b_allPoints->data;
    nz = r->size[0];
    for (lastBlockLength = 0; lastBlockLength < 2; lastBlockLength++) {
      for (jlast = 0; jlast < 2; jlast++) {
        for (nblocks = 0; nblocks < nz; nblocks++) {
          dis_data[(nblocks + b_allPoints->size[0] * jlast) +
                   b_allPoints->size[0] * 2 * lastBlockLength] =
              allPoints_data[(r1[nblocks] + allPoints->size[0] * jlast) +
                             allPoints->size[0] * 2 * lastBlockLength];
        }
      }
    }
    computeRigid2d(b_allPoints, modelParams);
    evaluateTform2d(modelParams, allPoints, dis);
    dis_data = dis->data;
    nz = dis->size[0] - 1;
    for (nblocks = 0; nblocks <= nz; nblocks++) {
      if (dis_data[nblocks] > 1.5) {
        dis_data[nblocks] = 1.5;
      }
    }
    bestModelParams_size[0] = 3;
    bestModelParams_size[1] = 3;
    memcpy(&bestModelParams_data[0], &modelParams[0], 9U * sizeof(double));
    lastBlockLength = inliers->size[0];
    inliers->size[0] = dis->size[0];
    emxEnsureCapacity_boolean_T(inliers, lastBlockLength);
    bestInliers_data = inliers->data;
    nz = dis->size[0];
    for (lastBlockLength = 0; lastBlockLength < nz; lastBlockLength++) {
      bestInliers_data[lastBlockLength] = (dis_data[lastBlockLength] < 1.5);
    }
    for (lastBlockLength = 0; lastBlockLength < 9; lastBlockLength++) {
      j = modelParams[lastBlockLength];
      b_data[lastBlockLength] = ((!rtIsInf(j)) && (!rtIsNaN(j)));
    }
    isValidModel = true;
    jlast = 1;
    exitg1 = false;
    while ((!exitg1) && (jlast <= 9)) {
      if (!b_data[jlast - 1]) {
        isValidModel = false;
        exitg1 = true;
      } else {
        jlast++;
      }
    }
    guard1 = false;
    if (!isValidModel) {
      guard1 = true;
    } else {
      isValidModel = false;
      jlast = 1;
      exitg1 = false;
      while ((!exitg1) && (jlast <= inliers->size[0])) {
        if (bestInliers_data[jlast - 1]) {
          isValidModel = true;
          exitg1 = true;
        } else {
          jlast++;
        }
      }
      if (!isValidModel) {
        guard1 = true;
      }
    }
    if (guard1) {
      isFound = false;
      lastBlockLength = inliers->size[0];
      inliers->size[0] = allPoints->size[0];
      emxEnsureCapacity_boolean_T(inliers, lastBlockLength);
      bestInliers_data = inliers->data;
      nz = allPoints->size[0];
      for (lastBlockLength = 0; lastBlockLength < nz; lastBlockLength++) {
        bestInliers_data[lastBlockLength] = false;
      }
    }
  } else {
    lastBlockLength = inliers->size[0];
    inliers->size[0] = allPoints->size[0];
    emxEnsureCapacity_boolean_T(inliers, lastBlockLength);
    bestInliers_data = inliers->data;
    nz = allPoints->size[0];
    for (lastBlockLength = 0; lastBlockLength < nz; lastBlockLength++) {
      bestInliers_data[lastBlockLength] = false;
    }
  }
  emxFree_real_T(&b_allPoints);
  emxFree_int32_T(&r);
  emxFree_real_T(&dis);
  emxFree_boolean_T(&bestInliers);
  return isFound;
}

static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
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
    y = atan2(i, i1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }
  return y;
}

static double rt_hypotd_snf(double u0, double u1)
{
  double a;
  double b;
  double y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = rtNaN;
  } else {
    y = a * 1.4142135623730951;
  }
  return y;
}

static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    double d;
    double d1;
    d = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
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
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }
  return y;
}

static double rt_remd_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = rtNaN;
  } else if (rtIsInf(u1)) {
    y = u0;
  } else if ((u1 != 0.0) && (u1 != trunc(u1))) {
    double q;
    q = fabs(u0 / u1);
    if (!(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q)) {
      y = 0.0 * u0;
    } else {
      y = fmod(u0, u1);
    }
  } else {
    y = fmod(u0, u1);
  }
  return y;
}

static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }
  return y;
}

static void svd(const double A[4], double U[4], double s[2], double V[4])
{
  double b_s[2];
  double e[2];
  double A_idx_3;
  double f;
  double nrm;
  double rt;
  double sm;
  double snorm;
  double sqds;
  double temp;
  int iter;
  int k;
  int kase;
  int m;
  int q;
  int qs;
  temp = A[0];
  sm = A[1];
  sqds = A[2];
  A_idx_3 = A[3];
  nrm = xnrm2(A);
  if (nrm > 0.0) {
    if (A[0] < 0.0) {
      b_s[0] = -nrm;
    } else {
      b_s[0] = nrm;
    }
    if (fabs(b_s[0]) >= 1.0020841800044864E-292) {
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
      rt = fabs(nrm);
      nrm /= rt;
      b_s[q] = rt;
      if (q + 1 < 2) {
        e[0] /= nrm;
      }
      kase = q << 1;
      qs = kase + 2;
      for (k = kase + 1; k <= qs; k++) {
        U[k - 1] *= nrm;
      }
    }
    if ((q + 1 < 2) && (e[0] != 0.0)) {
      rt = fabs(e[0]);
      nrm = rt / e[0];
      e[0] = rt;
      b_s[1] *= nrm;
      V[2] *= nrm;
      V[3] *= nrm;
    }
  }
  iter = 0;
  snorm =
      fmax(fmax(0.0, fmax(fabs(b_s[0]), fabs(e[0]))), fmax(fabs(b_s[1]), 0.0));
  while ((m > 0) && (iter < 75)) {
    int ii_tmp_tmp;
    boolean_T exitg1;
    ii_tmp_tmp = m - 1;
    q = m - 1;
    exitg1 = false;
    while (!(exitg1 || (q == 0))) {
      nrm = fabs(e[0]);
      if ((nrm <= 2.2204460492503131E-16 * (fabs(b_s[0]) + fabs(b_s[1]))) ||
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
            nrm = fabs(e[0]);
          }
          if (kase > q + 1) {
            nrm += fabs(e[0]);
          }
          rt = fabs(b_s[kase - 1]);
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
    case 1:
      f = e[0];
      e[0] = 0.0;
      for (k = ii_tmp_tmp; k >= q + 1; k--) {
        A_idx_3 = xrotg(&b_s[0], &f, &rt);
        kase = (m - 1) << 1;
        temp = A_idx_3 * V[0] + rt * V[kase];
        V[kase] = A_idx_3 * V[kase] - rt * V[0];
        V[0] = temp;
        nrm = V[kase + 1];
        temp = A_idx_3 * V[1] + rt * nrm;
        V[kase + 1] = A_idx_3 * nrm - rt * V[1];
        V[1] = temp;
      }
      break;
    case 2:
      f = e[q - 1];
      e[q - 1] = 0.0;
      for (k = q + 1; k <= m; k++) {
        A_idx_3 = xrotg(&b_s[k - 1], &f, &sm);
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
      break;
    case 3: {
      double scale;
      nrm = b_s[m - 1];
      scale = fmax(
          fmax(fmax(fmax(fabs(nrm), fabs(b_s[0])), fabs(e[0])), fabs(b_s[q])),
          fabs(e[q]));
      sm = nrm / scale;
      rt = b_s[0] / scale;
      nrm = e[0] / scale;
      sqds = b_s[q] / scale;
      temp = ((rt + sm) * (rt - sm) + nrm * nrm) / 2.0;
      A_idx_3 = sm * nrm;
      A_idx_3 *= A_idx_3;
      if ((temp != 0.0) || (A_idx_3 != 0.0)) {
        rt = sqrt(temp * temp + A_idx_3);
        if (temp < 0.0) {
          rt = -rt;
        }
        rt = A_idx_3 / (temp + rt);
      } else {
        rt = 0.0;
      }
      f = (sqds + sm) * (sqds - sm) + rt;
      rt = sqds * (e[q] / scale);
      for (k = q + 1; k < 2; k++) {
        A_idx_3 = xrotg(&f, &rt, &sm);
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
        A_idx_3 = xrotg(&b_s[0], &rt, &sm);
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
        for (k = kase + 1; k <= qs; k++) {
          V[k - 1] = -V[k - 1];
        }
      }
      while ((q + 1 < 2) && (b_s[0] < b_s[1])) {
        rt = b_s[0];
        b_s[0] = b_s[1];
        b_s[1] = rt;
        xswap(V);
        xswap(U);
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

static double xnrm2(const double x[4])
{
  double absxk;
  double scale;
  double t;
  double y;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }
  absxk = fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  return scale * sqrt(y);
}

static double xrotg(double *a, double *b, double *s)
{
  double absa;
  double absb;
  double c;
  double roe;
  double scale;
  roe = *b;
  absa = fabs(*a);
  absb = fabs(*b);
  if (absa > absb) {
    roe = *a;
  }
  scale = absa + absb;
  if (scale == 0.0) {
    *s = 0.0;
    c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    double ads;
    double bds;
    ads = absa / scale;
    bds = absb / scale;
    scale *= sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      scale = -scale;
    }
    c = *a / scale;
    *s = *b / scale;
    if (absa > absb) {
      *b = *s;
    } else if (c != 0.0) {
      *b = 1.0 / c;
    } else {
      *b = 1.0;
    }
    *a = scale;
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

void estimateAffineRigid2D(const emxArray_real_T *pts1,
                           const emxArray_real_T *pts2, double tform2x3[6],
                           emxArray_boolean_T *inlierIndex, int *status)
{
  static const signed char iv[3] = {0, 0, 1};
  emxArray_boolean_T *inlierIdx;
  emxArray_real_T *points;
  double failedMatrix[9];
  double tmatrix_data[9];
  double x_data[9];
  const double *pts1_data;
  const double *pts2_data;
  double r;
  double s;
  double smax;
  double wpr;
  double *points_data;
  int tmatrix_size[2];
  int i;
  int i1;
  int ijA;
  int j;
  int k;
  int yk;
  boolean_T isodd;
  boolean_T *inlierIdx_data;
  boolean_T *inlierIndex_data;
  if (!isInitialized_estimateAffineRigid2D) {
    estimateAffineRigid2D_initialize();
  }
  pts2_data = pts2->data;
  pts1_data = pts1->data;
  *status = (pts1->size[0] < 2);
  memset(&failedMatrix[0], 0, 9U * sizeof(double));
  failedMatrix[0] = 1.0;
  failedMatrix[4] = 1.0;
  failedMatrix[8] = 1.0;
  if (*status == 0) {
    boolean_T guard1 = false;
    emxInit_real_T(&points, 3);
    i = points->size[0] * points->size[1] * points->size[2];
    points->size[0] = pts1->size[0];
    points->size[1] = 2;
    points->size[2] = 2;
    emxEnsureCapacity_real_T(points, i);
    points_data = points->data;
    i = pts1->size[0] << 1;
    for (j = 0; j < i; j++) {
      points_data[j] = pts1_data[j];
    }
    i1 = pts2->size[0] << 1;
    for (j = 0; j < i1; j++) {
      points_data[i + j] = pts2_data[j];
    }
    emxInit_boolean_T(&inlierIdx, 1);
    isodd = msac(points, tmatrix_data, tmatrix_size, inlierIdx);
    inlierIdx_data = inlierIdx->data;
    emxFree_real_T(&points);
    i = inlierIndex->size[0] * inlierIndex->size[1];
    inlierIndex->size[0] = inlierIdx->size[0];
    inlierIndex->size[1] = 1;
    emxEnsureCapacity_boolean_T(inlierIndex, i);
    inlierIndex_data = inlierIndex->data;
    yk = inlierIdx->size[0];
    for (i = 0; i < yk; i++) {
      inlierIndex_data[i] = inlierIdx_data[i];
    }
    if (!isodd) {
      *status = 2;
    }
    if ((tmatrix_size[0] == 0) || (tmatrix_size[1] == 0)) {
      smax = 1.0;
    } else {
      int m;
      int n;
      int u1;
      int x_size_idx_0;
      signed char ipiv_data[3];
      m = tmatrix_size[0];
      n = tmatrix_size[1] - 2;
      x_size_idx_0 = tmatrix_size[0];
      yk = tmatrix_size[0] * tmatrix_size[1];
      memcpy(&x_data[0], &tmatrix_data[0], (unsigned int)yk * sizeof(double));
      yk = tmatrix_size[0];
      u1 = tmatrix_size[1];
      if (yk <= u1) {
        u1 = yk;
      }
      ipiv_data[0] = 1;
      yk = 1;
      for (k = 2; k <= u1; k++) {
        yk++;
        ipiv_data[k - 1] = (signed char)yk;
      }
      if (tmatrix_size[0] - 1 <= tmatrix_size[1]) {
        i = tmatrix_size[0];
      } else {
        i = 2;
      }
      for (j = 0; j <= i - 2; j++) {
        int b_tmp;
        int ipiv_tmp;
        int jA;
        int jp1j;
        int mmj;
        mmj = m - j;
        b_tmp = j * (m + 1);
        jp1j = b_tmp + 2;
        if (mmj < 1) {
          yk = -1;
        } else {
          yk = 0;
          if (mmj > 1) {
            smax = fabs(x_data[b_tmp]);
            for (k = 2; k <= mmj; k++) {
              s = fabs(x_data[(b_tmp + k) - 1]);
              if (s > smax) {
                yk = k - 1;
                smax = s;
              }
            }
          }
        }
        if (x_data[b_tmp + yk] != 0.0) {
          if (yk != 0) {
            ipiv_tmp = j + yk;
            ipiv_data[j] = (signed char)(ipiv_tmp + 1);
            for (k = 0; k <= n + 1; k++) {
              yk = k * m;
              jA = j + yk;
              smax = x_data[jA];
              i1 = ipiv_tmp + yk;
              x_data[jA] = x_data[i1];
              x_data[i1] = smax;
            }
          }
          i1 = b_tmp + mmj;
          for (yk = jp1j; yk <= i1; yk++) {
            x_data[yk - 1] /= x_data[b_tmp];
          }
        }
        yk = n - j;
        ipiv_tmp = b_tmp + m;
        jA = ipiv_tmp;
        for (k = 0; k <= yk; k++) {
          smax = x_data[ipiv_tmp + k * m];
          if (smax != 0.0) {
            i1 = jA + 2;
            jp1j = mmj + jA;
            for (ijA = i1; ijA <= jp1j; ijA++) {
              x_data[ijA - 1] += x_data[((b_tmp + ijA) - jA) - 1] * -smax;
            }
          }
          jA += m;
        }
      }
      smax = x_data[0];
      for (k = 0; k <= x_size_idx_0 - 2; k++) {
        smax *= x_data[(k + x_size_idx_0 * (k + 1)) + 1];
      }
      isodd = false;
      for (k = 0; k <= u1 - 2; k++) {
        if (ipiv_data[k] > k + 1) {
          isodd = !isodd;
        }
      }
      if (isodd) {
        smax = -smax;
      }
    }
    guard1 = false;
    if (smax == 0.0) {
      guard1 = true;
    } else {
      boolean_T exitg1;
      yk = tmatrix_size[0] * tmatrix_size[1];
      i = inlierIdx->size[0];
      inlierIdx->size[0] = yk;
      emxEnsureCapacity_boolean_T(inlierIdx, i);
      inlierIdx_data = inlierIdx->data;
      for (i = 0; i < yk; i++) {
        smax = tmatrix_data[i];
        inlierIdx_data[i] = (rtIsInf(smax) || rtIsNaN(smax));
      }
      isodd = false;
      yk = 1;
      exitg1 = false;
      while ((!exitg1) && (yk <= inlierIdx->size[0])) {
        if (inlierIdx_data[yk - 1]) {
          isodd = true;
          exitg1 = true;
        } else {
          yk++;
        }
      }
      if (isodd) {
        guard1 = true;
      }
    }
    if (guard1) {
      *status = 2;
      tmatrix_size[0] = 3;
      memcpy(&tmatrix_data[0], &failedMatrix[0], 9U * sizeof(double));
    }
    emxFree_boolean_T(&inlierIdx);
  } else {
    i = inlierIndex->size[0] * inlierIndex->size[1];
    inlierIndex->size[0] = pts1->size[0];
    inlierIndex->size[1] = pts1->size[0];
    emxEnsureCapacity_boolean_T(inlierIndex, i);
    inlierIndex_data = inlierIndex->data;
    yk = pts1->size[0] * pts1->size[0];
    for (i = 0; i < yk; i++) {
      inlierIndex_data[i] = false;
    }
    tmatrix_size[0] = 3;
    memcpy(&tmatrix_data[0], &failedMatrix[0], 9U * sizeof(double));
  }
  if (*status != 0) {
    tmatrix_size[0] = 3;
    memcpy(&tmatrix_data[0], &failedMatrix[0], 9U * sizeof(double));
  }
  i = tmatrix_size[0];
  for (i1 = 0; i1 < 3; i1++) {
    failedMatrix[3 * i1] = tmatrix_data[i1];
    failedMatrix[3 * i1 + 1] = tmatrix_data[i1 + i];
    failedMatrix[3 * i1 + 2] = tmatrix_data[i1 + i * 2];
  }
  double b_failedMatrix[4];
  double dv[4];
  b_failedMatrix[0] = failedMatrix[0];
  b_failedMatrix[1] = failedMatrix[1];
  b_failedMatrix[2] = failedMatrix[3];
  b_failedMatrix[3] = failedMatrix[4];
  smax = constrainToRotationMatrix2D(b_failedMatrix, dv);
  if (rtIsNaN(smax + 180.0) || rtIsInf(smax + 180.0)) {
    r = rtNaN;
  } else if (smax + 180.0 == 0.0) {
    r = 0.0;
  } else {
    r = fmod(smax + 180.0, 360.0);
    if (r == 0.0) {
      r = 0.0;
    } else if (smax + 180.0 < 0.0) {
      r += 360.0;
    }
  }
  wpr = rt_roundd_snf(r - 180.0);
  if (r - 180.0 == wpr) {
    isodd = true;
  } else {
    s = fabs((r - 180.0) - wpr);
    if ((r - 180.0 == 0.0) || (wpr == 0.0)) {
      isodd = (s < 4.94065645841247E-324);
    } else {
      smax = fabs(r - 180.0) + fabs(wpr);
      if (smax < 2.2250738585072014E-308) {
        isodd = (s < 4.94065645841247E-324);
      } else {
        isodd =
            (s / fmin(smax, 1.7976931348623157E+308) < 2.2204460492503131E-16);
      }
    }
  }
  s = r - 180.0;
  if (isodd) {
    s = wpr;
  }
  smax = s;
  b_cosd(&smax);
  b_sind(&s);
  tmatrix_data[0] = smax;
  tmatrix_data[3] = -s;
  tmatrix_data[6] = failedMatrix[6];
  tmatrix_data[1] = s;
  tmatrix_data[4] = smax;
  tmatrix_data[7] = failedMatrix[7];
  for (i = 0; i < 3; i++) {
    tmatrix_data[3 * i + 2] = iv[i];
    yk = i << 1;
    tform2x3[yk] = tmatrix_data[3 * i];
    tform2x3[yk + 1] = tmatrix_data[3 * i + 1];
  }
}

void estimateAffineRigid2D_initialize(void)
{
  c_eml_rand_mt19937ar_stateful_i();
  isInitialized_estimateAffineRigid2D = true;
}

void estimateAffineRigid2D_terminate(void)
{
  isInitialized_estimateAffineRigid2D = false;
}
