#include "selectUniform2.h"
#include "rt_nonfinite.h"
#include "selectUniform2_emxutil.h"
#include "selectUniform2_types.h"
#include "rt_nonfinite.h"
#include <math.h>

static void selectPoints(const emxArray_real_T *points,
                         const double imageSize[2],
                         const emxArray_real_T *metric, double numPoints,
                         emxArray_boolean_T *pointsIdx);

static void selectPoints(const emxArray_real_T *points,
                         const double imageSize[2],
                         const emxArray_real_T *metric, double numPoints,
                         emxArray_boolean_T *pointsIdx)
{
  emxArray_int32_T *r;
  emxArray_uint32_T *binIdx;
  const double *metric_data;
  const double *points_data;
  int b_i;
  int i;
  int loop_ub_tmp;
  unsigned int *binIdx_data;
  int *r1;
  boolean_T *pointsIdx_data;
  metric_data = metric->data;
  points_data = points->data;
  if (numPoints == 1.0) {
    int idx;
    int last;
    last = metric->size[0];
    if (metric->size[0] <= 2) {
      if (metric->size[0] == 1) {
        idx = 1;
      } else {
        double d;
        d = metric_data[metric->size[0] - 1];
        if ((metric_data[0] < d) ||
            (rtIsNaN(metric_data[0]) && (!rtIsNaN(d)))) {
          idx = metric->size[0];
        } else {
          idx = 1;
        }
      }
    } else {
      if (!rtIsNaN(metric_data[0])) {
        idx = 1;
      } else {
        boolean_T exitg1;
        idx = 0;
        loop_ub_tmp = 2;
        exitg1 = false;
        while ((!exitg1) && (loop_ub_tmp <= last)) {
          if (!rtIsNaN(metric_data[loop_ub_tmp - 1])) {
            idx = loop_ub_tmp;
            exitg1 = true;
          } else {
            loop_ub_tmp++;
          }
        }
      }
      if (idx == 0) {
        idx = 1;
      } else {
        double aspectRatio;
        aspectRatio = metric_data[idx - 1];
        i = idx + 1;
        for (loop_ub_tmp = i; loop_ub_tmp <= last; loop_ub_tmp++) {
          double d;
          d = metric_data[loop_ub_tmp - 1];
          if (aspectRatio < d) {
            aspectRatio = d;
            idx = loop_ub_tmp;
          }
        }
      }
    }
    i = pointsIdx->size[0];
    pointsIdx->size[0] = points->size[0];
    emxEnsureCapacity_boolean_T(pointsIdx, i);
    pointsIdx_data = pointsIdx->data;
    last = points->size[0];
    for (i = 0; i < last; i++) {
      pointsIdx_data[i] = false;
    }
    pointsIdx_data[idx - 1] = true;
  } else {
    double aspectRatio;
    double d;
    double gridStep_idx_1;
    double h;
    int last;
    aspectRatio = imageSize[0] / imageSize[1];
    h = fmax(floor(sqrt(numPoints / aspectRatio)), 1.0);
    d = fmax(floor(h * aspectRatio), 1.0);
    aspectRatio = imageSize[0] / (d + 1.0);
    gridStep_idx_1 = imageSize[1] / (h + 1.0);
    emxInit_uint32_T(&binIdx, 2);
    i = binIdx->size[0] * binIdx->size[1];
    binIdx->size[0] = (int)d;
    binIdx->size[1] = (int)h;
    emxEnsureCapacity_uint32_T(binIdx, i);
    binIdx_data = binIdx->data;
    loop_ub_tmp = (int)d * (int)h;
    for (i = 0; i < loop_ub_tmp; i++) {
      binIdx_data[i] = 0U;
    }
    i = points->size[0];
    last = pointsIdx->size[0];
    pointsIdx->size[0] = points->size[0];
    emxEnsureCapacity_boolean_T(pointsIdx, last);
    pointsIdx_data = pointsIdx->data;
    for (b_i = 0; b_i < i; b_i++) {
      double d1;
      double whichBin_idx_0;
      double whichBin_idx_1;
      int idx;
      boolean_T p;
      d1 = floor(points_data[b_i] / aspectRatio);
      whichBin_idx_0 = d1 + 1.0;
      if (rtIsNaN(d)) {
        p = false;
      } else if (rtIsNaN(d1 + 1.0)) {
        p = true;
      } else {
        p = (d1 + 1.0 > d);
      }
      if (p) {
        whichBin_idx_0 = d;
      }
      d1 = floor(points_data[b_i + points->size[0]] / gridStep_idx_1);
      whichBin_idx_1 = d1 + 1.0;
      if (rtIsNaN(h)) {
        p = false;
      } else if (rtIsNaN(d1 + 1.0)) {
        p = true;
      } else {
        p = (d1 + 1.0 > h);
      }
      if (p) {
        whichBin_idx_1 = h;
      }
      idx = (int)binIdx_data[((int)whichBin_idx_0 +
                              binIdx->size[0] * ((int)whichBin_idx_1 - 1)) -
                             1];
      if ((idx < 1) || (metric_data[idx - 1] < metric_data[b_i])) {
        binIdx_data[((int)whichBin_idx_0 +
                     binIdx->size[0] * ((int)whichBin_idx_1 - 1)) -
                    1] = (unsigned int)(b_i + 1);
      }
      pointsIdx_data[b_i] = false;
    }
    loop_ub_tmp--;
    last = 0;
    for (b_i = 0; b_i <= loop_ub_tmp; b_i++) {
      if ((int)binIdx_data[b_i] > 0) {
        last++;
      }
    }
    emxInit_int32_T(&r);
    i = r->size[0];
    r->size[0] = last;
    emxEnsureCapacity_int32_T(r, i);
    r1 = r->data;
    last = 0;
    for (b_i = 0; b_i <= loop_ub_tmp; b_i++) {
      if ((int)binIdx_data[b_i] > 0) {
        r1[last] = (int)binIdx_data[b_i];
        last++;
      }
    }
    emxFree_uint32_T(&binIdx);
    last = r->size[0];
    for (i = 0; i < last; i++) {
      pointsIdx_data[r1[i] - 1] = true;
    }
    emxFree_int32_T(&r);
  }
}

void selectUniform2(emxArray_real_T *points, emxArray_real_T *responses,
                    double N, const double imageSize[2],
                    emxArray_real_T *pointsOut, emxArray_real_T *b_index)
{
  emxArray_boolean_T *idx;
  emxArray_int32_T *r;
  emxArray_int32_T *r2;
  emxArray_int32_T *r3;
  emxArray_real_T *b_points;
  emxArray_real_T *b_responses;
  emxArray_uint32_T *b_origIdx;
  emxArray_uint32_T *idxOut;
  emxArray_uint32_T *origIdx;
  double first;
  double *pointsOut_data;
  double *points_data;
  double *responses_data;
  int end;
  int i;
  int u0;
  unsigned int *b_origIdx_data;
  unsigned int *idxOut_data;
  unsigned int *origIdx_data;
  int *r1;
  boolean_T *idx_data;
  responses_data = responses->data;
  points_data = points->data;
  emxInit_uint32_T(&idxOut, 1);
  i = idxOut->size[0];
  idxOut->size[0] = responses->size[0];
  emxEnsureCapacity_uint32_T(idxOut, i);
  idxOut_data = idxOut->data;
  u0 = points->size[0];
  if (u0 < 2) {
    u0 = 2;
  }
  if (points->size[0] == 0) {
    u0 = 0;
  }
  emxInit_uint32_T(&origIdx, 2);
  origIdx_data = origIdx->data;
  if (u0 < 1) {
    origIdx->size[0] = 1;
    origIdx->size[1] = 0;
  } else {
    i = origIdx->size[0] * origIdx->size[1];
    origIdx->size[0] = 1;
    origIdx->size[1] = u0;
    emxEnsureCapacity_uint32_T(origIdx, i);
    origIdx_data = origIdx->data;
    u0--;
    for (i = 0; i <= u0; i++) {
      origIdx_data[i] = (unsigned int)i + 1U;
    }
  }
  emxInit_boolean_T(&idx);
  selectPoints(points, imageSize, responses, N, idx);
  idx_data = idx->data;
  end = idx->size[0] - 1;
  u0 = 0;
  for (i = 0; i <= end; i++) {
    if (idx_data[i]) {
      u0++;
    }
  }
  emxInit_int32_T(&r);
  i = r->size[0];
  r->size[0] = u0;
  emxEnsureCapacity_int32_T(r, i);
  r1 = r->data;
  u0 = 0;
  for (i = 0; i <= end; i++) {
    if (idx_data[i]) {
      r1[u0] = i;
      u0++;
    }
  }
  if (r->size[0] < 1) {
    u0 = 0;
  } else {
    u0 = r->size[0];
  }
  for (i = 0; i < u0; i++) {
    idxOut_data[i] = origIdx_data[r1[i]];
  }
  first = (double)r->size[0] + 1.0;
  emxFree_int32_T(&r);
  emxInit_int32_T(&r2);
  emxInit_int32_T(&r3);
  emxInit_uint32_T(&b_origIdx, 2);
  emxInit_real_T(&b_points, 2);
  emxInit_real_T(&b_responses, 1);
  while (first <= N) {
    double d;
    end = idx->size[0] - 1;
    u0 = 0;
    for (i = 0; i <= end; i++) {
      if (!idx_data[i]) {
        u0++;
      }
    }
    i = r3->size[0];
    r3->size[0] = u0;
    emxEnsureCapacity_int32_T(r3, i);
    r1 = r3->data;
    u0 = 0;
    for (i = 0; i <= end; i++) {
      if (!idx_data[i]) {
        r1[u0] = i;
        u0++;
      }
    }
    i = b_origIdx->size[0] * b_origIdx->size[1];
    b_origIdx->size[0] = 1;
    b_origIdx->size[1] = r3->size[0];
    emxEnsureCapacity_uint32_T(b_origIdx, i);
    b_origIdx_data = b_origIdx->data;
    u0 = r3->size[0];
    for (i = 0; i < u0; i++) {
      b_origIdx_data[i] = origIdx_data[r1[i]];
    }
    i = origIdx->size[0] * origIdx->size[1];
    origIdx->size[0] = 1;
    origIdx->size[1] = b_origIdx->size[1];
    emxEnsureCapacity_uint32_T(origIdx, i);
    origIdx_data = origIdx->data;
    u0 = b_origIdx->size[1];
    for (i = 0; i < u0; i++) {
      origIdx_data[i] = b_origIdx_data[i];
    }
    i = b_points->size[0] * b_points->size[1];
    b_points->size[0] = r3->size[0];
    b_points->size[1] = 2;
    emxEnsureCapacity_real_T(b_points, i);
    pointsOut_data = b_points->data;
    u0 = r3->size[0];
    for (i = 0; i < 2; i++) {
      for (end = 0; end < u0; end++) {
        pointsOut_data[end + b_points->size[0] * i] =
            points_data[r1[end] + points->size[0] * i];
      }
    }
    i = points->size[0] * points->size[1];
    points->size[0] = b_points->size[0];
    points->size[1] = 2;
    emxEnsureCapacity_real_T(points, i);
    points_data = points->data;
    u0 = b_points->size[0] << 1;
    for (i = 0; i < u0; i++) {
      points_data[i] = pointsOut_data[i];
    }
    i = b_responses->size[0];
    b_responses->size[0] = r3->size[0];
    emxEnsureCapacity_real_T(b_responses, i);
    pointsOut_data = b_responses->data;
    u0 = r3->size[0];
    for (i = 0; i < u0; i++) {
      pointsOut_data[i] = responses_data[r1[i]];
    }
    i = responses->size[0];
    responses->size[0] = b_responses->size[0];
    emxEnsureCapacity_real_T(responses, i);
    responses_data = responses->data;
    u0 = b_responses->size[0];
    for (i = 0; i < u0; i++) {
      responses_data[i] = pointsOut_data[i];
    }
    selectPoints(points, imageSize, responses, N - (first - 1.0), idx);
    idx_data = idx->data;
    end = idx->size[0] - 1;
    u0 = 0;
    for (i = 0; i <= end; i++) {
      if (idx_data[i]) {
        u0++;
      }
    }
    i = r2->size[0];
    r2->size[0] = u0;
    emxEnsureCapacity_int32_T(r2, i);
    r1 = r2->data;
    u0 = 0;
    for (i = 0; i <= end; i++) {
      if (idx_data[i]) {
        r1[u0] = i;
        u0++;
      }
    }
    d = (first + (double)r2->size[0]) - 1.0;
    if (first > d) {
      i = 0;
      end = 0;
    } else {
      i = (int)first - 1;
      end = (int)d;
    }
    u0 = end - i;
    for (end = 0; end < u0; end++) {
      idxOut_data[i + end] = origIdx_data[r1[end]];
    }
    first += (double)r2->size[0];
  }
  emxFree_real_T(&b_responses);
  emxFree_real_T(&b_points);
  emxFree_uint32_T(&b_origIdx);
  emxFree_int32_T(&r3);
  emxFree_int32_T(&r2);
  emxFree_boolean_T(&idx);
  emxFree_uint32_T(&origIdx);
  if (N < 1.0) {
    u0 = 0;
  } else {
    u0 = (int)N;
  }
  i = pointsOut->size[0] * pointsOut->size[1];
  pointsOut->size[0] = u0;
  pointsOut->size[1] = 2;
  emxEnsureCapacity_real_T(pointsOut, i);
  pointsOut_data = pointsOut->data;
  for (i = 0; i < 2; i++) {
    for (end = 0; end < u0; end++) {
      pointsOut_data[end + pointsOut->size[0] * i] =
          points_data[((int)idxOut_data[end] + points->size[0] * i) - 1];
    }
  }
  i = b_index->size[0];
  b_index->size[0] = u0;
  emxEnsureCapacity_real_T(b_index, i);
  pointsOut_data = b_index->data;
  for (i = 0; i < u0; i++) {
    pointsOut_data[i] = idxOut_data[i];
  }
  emxFree_uint32_T(&idxOut);
}

void selectUniform2_initialize(void)
{
}

void selectUniform2_terminate(void)
{
}
