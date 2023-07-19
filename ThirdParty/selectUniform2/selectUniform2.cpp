#include "selectUniform2.h"
#include "coder_array.h"
#include <cmath>

namespace selectUniform2 {
static void selectPoints(const ::coder::array<double, 2U> &points,
                         const double imageSize[2],
                         const ::coder::array<double, 1U> &metric,
                         double numPoints,
                         ::coder::array<boolean_T, 1U> &pointsIdx);

}

namespace selectUniform2 {
static void selectPoints(const ::coder::array<double, 2U> &points,
                         const double imageSize[2],
                         const ::coder::array<double, 1U> &metric,
                         double numPoints,
                         ::coder::array<boolean_T, 1U> &pointsIdx)
{
  ::coder::array<unsigned int, 2U> binIdx;
  ::coder::array<int, 1U> r;
  if (numPoints == 1.0) {
    int b_i;
    int idx;
    int last;
    last = metric.size(0);
    if (metric.size(0) <= 2) {
      if (metric.size(0) == 1) {
        idx = 1;
      } else {
        double d;
        d = metric[metric.size(0) - 1];
        if ((metric[0] < d) || (std::isnan(metric[0]) && (!std::isnan(d)))) {
          idx = metric.size(0);
        } else {
          idx = 1;
        }
      }
    } else {
      int i;
      if (!std::isnan(metric[0])) {
        idx = 1;
      } else {
        boolean_T exitg1;
        idx = 0;
        i = 2;
        exitg1 = false;
        while ((!exitg1) && (i <= last)) {
          if (!std::isnan(metric[i - 1])) {
            idx = i;
            exitg1 = true;
          } else {
            i++;
          }
        }
      }
      if (idx == 0) {
        idx = 1;
      } else {
        double aspectRatio;
        aspectRatio = metric[idx - 1];
        b_i = idx + 1;
        for (i = b_i; i <= last; i++) {
          double d;
          d = metric[i - 1];
          if (aspectRatio < d) {
            aspectRatio = d;
            idx = i;
          }
        }
      }
    }
    pointsIdx.set_size(points.size(0));
    last = points.size(0);
    for (b_i = 0; b_i < last; b_i++) {
      pointsIdx[b_i] = false;
    }
    pointsIdx[idx - 1] = true;
  } else {
    double aspectRatio;
    double d;
    double gridStep_idx_1;
    double h;
    int b_i;
    int idx;
    int last;
    aspectRatio = imageSize[0] / imageSize[1];
    h = std::fmax(std::floor(std::sqrt(numPoints / aspectRatio)), 1.0);
    d = std::fmax(std::floor(h * aspectRatio), 1.0);
    aspectRatio = imageSize[0] / (d + 1.0);
    gridStep_idx_1 = imageSize[1] / (h + 1.0);
    binIdx.set_size(static_cast<int>(d), static_cast<int>(h));
    last = static_cast<int>(d) * static_cast<int>(h);
    for (b_i = 0; b_i < last; b_i++) {
      binIdx[b_i] = 0U;
    }
    b_i = points.size(0);
    pointsIdx.set_size(points.size(0));
    for (int i{0}; i < b_i; i++) {
      double d1;
      double whichBin_idx_0;
      double whichBin_idx_1;
      boolean_T p;
      d1 = std::floor(points[i] / aspectRatio);
      whichBin_idx_0 = d1 + 1.0;
      if (std::isnan(d)) {
        p = false;
      } else if (std::isnan(d1 + 1.0)) {
        p = true;
      } else {
        p = (d1 + 1.0 > d);
      }
      if (p) {
        whichBin_idx_0 = d;
      }
      d1 = std::floor(points[i + points.size(0)] / gridStep_idx_1);
      whichBin_idx_1 = d1 + 1.0;
      if (std::isnan(h)) {
        p = false;
      } else if (std::isnan(d1 + 1.0)) {
        p = true;
      } else {
        p = (d1 + 1.0 > h);
      }
      if (p) {
        whichBin_idx_1 = h;
      }
      idx = static_cast<int>(
          binIdx[(static_cast<int>(whichBin_idx_0) +
                  binIdx.size(0) * (static_cast<int>(whichBin_idx_1) - 1)) -
                 1]);
      if ((idx < 1) || (metric[idx - 1] < metric[i])) {
        binIdx[(static_cast<int>(whichBin_idx_0) +
                binIdx.size(0) * (static_cast<int>(whichBin_idx_1) - 1)) -
               1] = static_cast<unsigned int>(i + 1);
      }
      pointsIdx[i] = false;
    }
    idx = last - 1;
    last = 0;
    for (int i{0}; i <= idx; i++) {
      if (static_cast<int>(binIdx[i]) > 0) {
        last++;
      }
    }
    r.set_size(last);
    last = 0;
    for (int i{0}; i <= idx; i++) {
      if (static_cast<int>(binIdx[i]) > 0) {
        r[last] = static_cast<int>(binIdx[i]);
        last++;
      }
    }
    last = r.size(0);
    for (b_i = 0; b_i < last; b_i++) {
      pointsIdx[r[b_i] - 1] = true;
    }
  }
}

void selectUniform2(::coder::array<double, 2U> &points,
                    ::coder::array<double, 1U> &responses, double N,
                    const double imageSize[2],
                    ::coder::array<double, 2U> &pointsOut,
                    ::coder::array<double, 1U> &b_index)
{
  ::coder::array<double, 2U> b_points;
  ::coder::array<double, 1U> b_responses;
  ::coder::array<unsigned int, 2U> b_origIdx;
  ::coder::array<unsigned int, 2U> origIdx;
  ::coder::array<unsigned int, 1U> idxOut;
  ::coder::array<int, 1U> r;
  ::coder::array<int, 1U> r1;
  ::coder::array<int, 1U> r2;
  ::coder::array<boolean_T, 1U> idx;
  int end;
  int i;
  int u0;
  idxOut.set_size(responses.size(0));
  u0 = points.size(0);
  if (u0 < 2) {
    u0 = 2;
  }
  if (points.size(0) == 0) {
    u0 = 0;
  }
  if (u0 < 1) {
    origIdx.set_size(1, 0);
  } else {
    origIdx.set_size(1, u0);
    u0--;
    for (end = 0; end <= u0; end++) {
      origIdx[end] = static_cast<unsigned int>(end) + 1U;
    }
  }
  selectPoints(points, imageSize, responses, N, idx);
  end = idx.size(0) - 1;
  u0 = 0;
  for (i = 0; i <= end; i++) {
    if (idx[i]) {
      u0++;
    }
  }
  r.set_size(u0);
  u0 = 0;
  for (i = 0; i <= end; i++) {
    if (idx[i]) {
      r[u0] = i;
      u0++;
    }
  }
  if (r.size(0) < 1) {
    u0 = 0;
  } else {
    u0 = r.size(0);
  }
  for (end = 0; end < u0; end++) {
    idxOut[end] = origIdx[r[end]];
  }
  for (double b_first{static_cast<double>(r.size(0)) + 1.0}; b_first <= N;
       b_first += static_cast<double>(r2.size(0))) {
    double d;
    end = idx.size(0) - 1;
    u0 = 0;
    for (i = 0; i <= end; i++) {
      if (!idx[i]) {
        u0++;
      }
    }
    r1.set_size(u0);
    u0 = 0;
    for (i = 0; i <= end; i++) {
      if (!idx[i]) {
        r1[u0] = i;
        u0++;
      }
    }
    b_origIdx.set_size(1, r1.size(0));
    u0 = r1.size(0);
    for (end = 0; end < u0; end++) {
      b_origIdx[end] = origIdx[r1[end]];
    }
    origIdx.set_size(1, b_origIdx.size(1));
    u0 = b_origIdx.size(1);
    for (end = 0; end < u0; end++) {
      origIdx[end] = b_origIdx[end];
    }
    b_points.set_size(r1.size(0), 2);
    u0 = r1.size(0);
    for (end = 0; end < 2; end++) {
      for (i = 0; i < u0; i++) {
        b_points[i + b_points.size(0) * end] =
            points[r1[i] + points.size(0) * end];
      }
    }
    points.set_size(b_points.size(0), 2);
    u0 = b_points.size(0) << 1;
    for (end = 0; end < u0; end++) {
      points[end] = b_points[end];
    }
    b_responses.set_size(r1.size(0));
    u0 = r1.size(0);
    for (end = 0; end < u0; end++) {
      b_responses[end] = responses[r1[end]];
    }
    responses.set_size(b_responses.size(0));
    u0 = b_responses.size(0);
    for (end = 0; end < u0; end++) {
      responses[end] = b_responses[end];
    }
    selectPoints(points, imageSize, responses, N - (b_first - 1.0), idx);
    end = idx.size(0) - 1;
    u0 = 0;
    for (i = 0; i <= end; i++) {
      if (idx[i]) {
        u0++;
      }
    }
    r2.set_size(u0);
    u0 = 0;
    for (i = 0; i <= end; i++) {
      if (idx[i]) {
        r2[u0] = i;
        u0++;
      }
    }
    d = (b_first + static_cast<double>(r2.size(0))) - 1.0;
    if (b_first > d) {
      end = 0;
      i = 0;
    } else {
      end = static_cast<int>(b_first) - 1;
      i = static_cast<int>(d);
    }
    u0 = i - end;
    for (i = 0; i < u0; i++) {
      idxOut[end + i] = origIdx[r2[i]];
    }
  }
  if (N < 1.0) {
    u0 = 0;
  } else {
    u0 = static_cast<int>(N);
  }
  pointsOut.set_size(u0, 2);
  for (end = 0; end < 2; end++) {
    for (i = 0; i < u0; i++) {
      pointsOut[i + pointsOut.size(0) * end] =
          points[(static_cast<int>(idxOut[i]) + points.size(0) * end) - 1];
    }
  }
  b_index.set_size(u0);
  for (end = 0; end < u0; end++) {
    b_index[end] = idxOut[end];
  }
}

void selectUniform2_initialize()
{
}

void selectUniform2_terminate()
{
}

} // namespace selectUniform2
