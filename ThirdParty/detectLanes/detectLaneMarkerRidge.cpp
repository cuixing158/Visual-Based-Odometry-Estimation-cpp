#include "detectLaneMarkerRidge.h"
#include "detectLaneMarkerRidge_types.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <cstring>

namespace detectLaneMarkerRidge {
omp_nest_lock_t detectLaneMarkerRidge_nestLockGlobal;

static boolean_T isInitialized_detectLaneMarkerRidge{false};

} // namespace detectLaneMarkerRidge

namespace detectLaneMarkerRidge {
static void b_binary_expand_op(::coder::array<short, 2U> &in1,
                               const ::coder::array<float, 2U> &in2,
                               const ::coder::array<unsigned char, 2U> &in3);

static void b_binary_expand_op(::coder::array<float, 2U> &in1,
                               const ::coder::array<unsigned char, 2U> &in2,
                               const ::coder::array<unsigned char, 2U> &in3,
                               const ::coder::array<unsigned char, 2U> &in4);

static void binary_expand_op(::coder::array<boolean_T, 2U> &in1,
                             const ::coder::array<float, 2U> &in2,
                             const ::coder::array<double, 2U> &in3,
                             const ::coder::array<unsigned char, 2U> &in4);

namespace coder {
static void bwareaopen(::coder::array<boolean_T, 2U> &varargin_1,
                       double varargin_2);

static void findNonZero(const ::coder::array<boolean_T, 2U> &BW, double numRow,
                        double numCol,
                        ::coder::array<int, 2U> &nonZeroPixelMatrix,
                        ::coder::array<int, 1U> &numNonZeros);

static void findNonZeroOmp(const ::coder::array<boolean_T, 2U> &BW,
                           ::coder::array<int, 2U> &nonZeroPixelMatrix,
                           ::coder::array<int, 1U> &tempNumsVector);

static int getHoughPixelsOmp(const ::coder::array<int, 2U> &nonZero,
                             double firstRho, double slope, int peak1,
                             int peak2,
                             const ::coder::array<int, 1U> &tempNumsVector,
                             ::coder::array<int, 2U> &houghPix);

static int getLocationOfMax(const ::coder::array<double, 2U> &A, int &jMax);

static void houghlines(const ::coder::array<boolean_T, 2U> &varargin_1,
                       const ::coder::array<double, 2U> &varargin_3,
                       const double varargin_4_data[],
                       const int varargin_4_size[2],
                       ::coder::array<struct0_T, 2U> &lines);

static void houghpeaks(const ::coder::array<double, 2U> &varargin_1,
                       double varargin_4, double peaks_data[],
                       int peaks_size[2]);

static double labelingWu_parallel(const ::coder::array<boolean_T, 2U> &im,
                                  double M, double N,
                                  ::coder::array<double, 2U> &L);

static void sortrows(::coder::array<int, 2U> &y, const double varargin_1[2]);

static void
standardHoughTransformOptimized(const ::coder::array<boolean_T, 2U> &BW,
                                const ::coder::array<double, 2U> &rho,
                                double numRow, double numCol,
                                ::coder::array<double, 2U> &H);

} // namespace coder
} // namespace detectLaneMarkerRidge

namespace detectLaneMarkerRidge {
static void b_binary_expand_op(::coder::array<short, 2U> &in1,
                               const ::coder::array<float, 2U> &in2,
                               const ::coder::array<unsigned char, 2U> &in3)
{
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }
  in1.set_size(loop_ub, in1.size(1));
  if (in3.size(1) == 1) {
    b_loop_ub = in2.size(1);
  } else {
    b_loop_ub = in3.size(1);
  }
  in1.set_size(in1.size(0), b_loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_0 = (in3.size(0) != 1);
  stride_1_1 = (in3.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (int i{0}; i < b_loop_ub; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      in1[i1 + in1.size(0) * i] = static_cast<short>(
          static_cast<short>(in2[i1 * stride_0_0 + in2.size(0) * aux_0_1]) -
          in3[i1 * stride_1_0 + in3.size(0) * aux_1_1]);
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

static void b_binary_expand_op(::coder::array<float, 2U> &in1,
                               const ::coder::array<unsigned char, 2U> &in2,
                               const ::coder::array<unsigned char, 2U> &in3,
                               const ::coder::array<unsigned char, 2U> &in4)
{
  ::coder::array<float, 2U> r;
  int aux_0_1;
  int aux_1_1;
  int aux_2_1;
  int aux_3_1;
  int b_loop_ub;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  int stride_2_0;
  int stride_2_1;
  int stride_3_0;
  int stride_3_1;
  if (in3.size(0) == 1) {
    i = in1.size(0);
  } else {
    i = in3.size(0);
  }
  if (in4.size(0) == 1) {
    if (i == 1) {
      loop_ub = in2.size(0);
    } else {
      loop_ub = i;
    }
  } else {
    loop_ub = in4.size(0);
  }
  if (in3.size(1) == 1) {
    i = in1.size(1);
  } else {
    i = in3.size(1);
  }
  if (in4.size(1) == 1) {
    if (i == 1) {
      b_loop_ub = in2.size(1);
    } else {
      b_loop_ub = i;
    }
  } else {
    b_loop_ub = in4.size(1);
  }
  r.set_size(loop_ub, b_loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_0 = (in1.size(0) != 1);
  stride_1_1 = (in1.size(1) != 1);
  stride_2_0 = (in3.size(0) != 1);
  stride_2_1 = (in3.size(1) != 1);
  stride_3_0 = (in4.size(0) != 1);
  stride_3_1 = (in4.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  aux_2_1 = 0;
  aux_3_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      r[i1 + r.size(0) * i] = static_cast<short>(
          static_cast<short>(
              static_cast<short>(
                  (in2[i1 * stride_0_0 + in2.size(0) * aux_0_1] << 1) -
                  static_cast<short>(
                      in1[i1 * stride_1_0 + in1.size(0) * aux_1_1])) -
              in3[i1 * stride_2_0 + in3.size(0) * aux_2_1]) -
          in4[i1 * stride_3_0 + in4.size(0) * aux_3_1]);
    }
    aux_3_1 += stride_3_1;
    aux_2_1 += stride_2_1;
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(r.size(0), r.size(1));
  loop_ub = r.size(1);
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = r.size(0);
    for (int i1{0}; i1 < b_loop_ub; i1++) {
      in1[i1 + in1.size(0) * i] = r[i1 + r.size(0) * i];
    }
  }
}

static void binary_expand_op(::coder::array<boolean_T, 2U> &in1,
                             const ::coder::array<float, 2U> &in2,
                             const ::coder::array<double, 2U> &in3,
                             const ::coder::array<unsigned char, 2U> &in4)
{
  int aux_0_1;
  int aux_1_1;
  int aux_2_1;
  int b_loop_ub;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  int stride_2_0;
  int stride_2_1;
  if (in4.size(0) == 1) {
    if (in3.size(0) == 1) {
      loop_ub = in2.size(0);
    } else {
      loop_ub = in3.size(0);
    }
  } else {
    loop_ub = in4.size(0);
  }
  in1.set_size(loop_ub, in1.size(1));
  if (in4.size(1) == 1) {
    if (in3.size(1) == 1) {
      b_loop_ub = in2.size(1);
    } else {
      b_loop_ub = in3.size(1);
    }
  } else {
    b_loop_ub = in4.size(1);
  }
  in1.set_size(in1.size(0), b_loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_0 = (in3.size(0) != 1);
  stride_1_1 = (in3.size(1) != 1);
  stride_2_0 = (in4.size(0) != 1);
  stride_2_1 = (in4.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  aux_2_1 = 0;
  for (int i{0}; i < b_loop_ub; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      in1[i1 + in1.size(0) * i] =
          ((in2[i1 * stride_0_0 + in2.size(0) * aux_0_1] >
            in3[i1 * stride_1_0 + in3.size(0) * aux_1_1]) &&
           (in4[i1 * stride_2_0 + in4.size(0) * aux_2_1] != 0));
    }
    aux_2_1 += stride_2_1;
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

namespace coder {
static void bwareaopen(::coder::array<boolean_T, 2U> &varargin_1,
                       double varargin_2)
{
  ::coder::array<double, 2U> L;
  ::coder::array<int, 1U> regionLengths;
  double numComponents;
  int i;
  int nRows;
  int numPixels;
  int r;
  numComponents = 0.0;
  if ((varargin_1.size(0) == 0) || (varargin_1.size(1) == 0)) {
    unsigned int unnamed_idx_0;
    unsigned int unnamed_idx_1;
    unnamed_idx_0 = static_cast<unsigned int>(varargin_1.size(0));
    unnamed_idx_1 = static_cast<unsigned int>(varargin_1.size(1));
    L.set_size(static_cast<int>(unnamed_idx_0),
               static_cast<int>(unnamed_idx_1));
    nRows = static_cast<int>(unnamed_idx_0) * static_cast<int>(unnamed_idx_1);
    for (i = 0; i < nRows; i++) {
      L[i] = 0.0;
    }
  } else {
    numComponents =
        labelingWu_parallel(varargin_1, static_cast<double>(varargin_1.size(0)),
                            static_cast<double>(varargin_1.size(1)), L);
  }
  numPixels = L.size(0) * L.size(1);
  nRows = static_cast<int>(numComponents);
  regionLengths.set_size(nRows);
  for (i = 0; i < nRows; i++) {
    regionLengths[i] = 0;
  }
  for (nRows = 0; nRows < numPixels; nRows++) {
    i = static_cast<int>(L[nRows]);
    if (i > 0) {
      regionLengths[i - 1] = regionLengths[i - 1] + 1;
    }
  }
  nRows = varargin_1.size(0);
  numPixels = varargin_1.size(1) - 1;
#pragma omp parallel for num_threads(omp_get_max_threads()) private(r)

  for (int c = 0; c <= numPixels; c++) {
    for (r = 0; r < nRows; r++) {
      if (varargin_1[r + varargin_1.size(0) * c] &&
          (regionLengths[static_cast<int>(L[r + L.size(0) * c]) - 1] <
           varargin_2)) {
        varargin_1[r + varargin_1.size(0) * c] = false;
      }
    }
  }
}

static void findNonZero(const ::coder::array<boolean_T, 2U> &BW, double numRow,
                        double numCol,
                        ::coder::array<int, 2U> &nonZeroPixelMatrix,
                        ::coder::array<int, 1U> &numNonZeros)
{
  ::coder::array<unsigned int, 1U> tempBin;
  int i;
  int k;
  unsigned int tempNum;
  int ub_loop;
  nonZeroPixelMatrix.set_size(BW.size(0), BW.size(1));
  numNonZeros.set_size(static_cast<int>(numCol));
  ub_loop = static_cast<int>(numCol) - 1;
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
    tempBin, tempNum, k, i)

  for (int j = 0; j <= ub_loop; j++) {
    k = static_cast<int>(numRow);
    tempBin.set_size(k);
    tempNum = 0U;
    for (i = 0; i < k; i++) {
      if (BW[i + BW.size(0) * j]) {
        tempNum++;
        tempBin[static_cast<int>(tempNum) - 1] =
            static_cast<unsigned int>(i) + 1U;
      }
    }
    numNonZeros[j] = static_cast<int>(tempNum);
    k = 0;
    while ((k <= static_cast<int>(numRow) - 1) &&
           (!(static_cast<double>(k) + 1.0 > tempNum))) {
      nonZeroPixelMatrix[k + nonZeroPixelMatrix.size(0) * j] =
          static_cast<int>(tempBin[k]);
      k++;
    }
  }
}

static void findNonZeroOmp(const ::coder::array<boolean_T, 2U> &BW,
                           ::coder::array<int, 2U> &nonZeroPixelMatrix,
                           ::coder::array<int, 1U> &tempNumsVector)
{
  ::coder::array<int, 1U> tempBin;
  int i;
  int numRow;
  unsigned int tempNum;
  int ub_loop;
  numRow = BW.size(0);
  nonZeroPixelMatrix.set_size(BW.size(0), BW.size(1));
  tempNumsVector.set_size(BW.size(1));
  ub_loop = BW.size(1) - 1;
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
    tempBin, tempNum, i)

  for (int j = 0; j <= ub_loop; j++) {
    tempBin.set_size(numRow);
    tempNum = 0U;
    for (i = 0; i < numRow; i++) {
      if (BW[i + BW.size(0) * j]) {
        tempNum++;
        tempBin[static_cast<int>(tempNum) - 1] = i;
      }
    }
    tempNumsVector[j] = static_cast<int>(tempNum);
    i = 0;
    while ((i <= numRow - 1) && (i + 1 <= static_cast<int>(tempNum))) {
      nonZeroPixelMatrix[i + nonZeroPixelMatrix.size(0) * j] = tempBin[i];
      i++;
    }
  }
}

static int getHoughPixelsOmp(const ::coder::array<int, 2U> &nonZero,
                             double firstRho, double slope, int peak1,
                             int peak2,
                             const ::coder::array<int, 1U> &tempNumsVector,
                             ::coder::array<int, 2U> &houghPix)
{
  ::coder::array<int, 2U> houghPixTemp;
  ::coder::array<int, 1U> tempBin;
  ::coder::array<int, 1U> tempHoughPixNumsVector;
  double sortingOrder[2];
  double colMaxPrime;
  double colMin;
  double colMinPrime;
  double cosTheta;
  double rhoVal;
  double rowMaxPrime;
  double rowMin;
  double rowMinPrime;
  double thetaVal;
  int b_i;
  int b_j;
  int colMax;
  int i;
  int j;
  int k;
  unsigned int n;
  int numCol;
  int numHoughPix;
  int numHoughPixPrime;
  int numRow;
  int rhoVal_tmp;
  int rowMax;
  int tempNum;
  int ub_loop;
  numRow = nonZero.size(0);
  numCol = nonZero.size(1) - 1;
  numHoughPix = 0;
  thetaVal =
      ((static_cast<double>(peak2) - 1.0) - 90.0) * 3.1415926535897931 / 180.0;
  cosTheta = std::cos(thetaVal);
  thetaVal = std::sin(thetaVal);
  rowMax = 0;
  rowMin = rtInf;
  colMax = 0;
  colMin = rtInf;
  houghPixTemp.set_size(nonZero.size(0), nonZero.size(1));
  tempHoughPixNumsVector.set_size(nonZero.size(1));
  ub_loop = nonZero.size(1) - 1;
#pragma omp parallel num_threads(omp_get_max_threads()) private(               \
    tempBin, colMinPrime, colMaxPrime, rowMinPrime, rowMaxPrime,               \
    numHoughPixPrime, tempNum, rhoVal, b_i, b_j, rhoVal_tmp)
  {
    numHoughPixPrime = 0;
    rowMaxPrime = rtMinusInf;
    rowMinPrime = rtInf;
    colMaxPrime = rtMinusInf;
    colMinPrime = rtInf;
#pragma omp for nowait
    for (k = 0; k <= ub_loop; k++) {
      tempBin.set_size(numRow);
      tempNum = 0;
      b_i = tempNumsVector[k];
      for (b_j = 0; b_j < b_i; b_j++) {
        rhoVal_tmp = nonZero[b_j + nonZero.size(0) * k];
        rhoVal = ((static_cast<double>(k) + 1.0) - 1.0) * cosTheta +
                 static_cast<double>(rhoVal_tmp) * thetaVal;
        if (static_cast<int>((slope * (rhoVal - firstRho) + 1.0) + 0.5) ==
            peak1) {
          tempNum++;
          tempBin[tempNum - 1] = rhoVal_tmp + 1;
        }
      }
      tempHoughPixNumsVector[k] = tempNum;
      numHoughPixPrime += tempNum;
      if (tempNum != 0) {
        rowMaxPrime =
            std::fmax(rowMaxPrime, static_cast<double>(tempBin[tempNum - 1]));
        rowMinPrime = std::fmin(rowMinPrime, static_cast<double>(tempBin[0]));
        colMaxPrime = std::fmax(colMaxPrime, static_cast<double>(k) + 1.0);
        colMinPrime = std::fmin(colMinPrime, static_cast<double>(k) + 1.0);
      }
      b_i = 0;
      while ((b_i <= numRow - 1) && (b_i + 1 <= tempNum)) {
        houghPixTemp[b_i + houghPixTemp.size(0) * k] = tempBin[b_i];
        b_i++;
      }
    }
    omp_set_nest_lock(&detectLaneMarkerRidge_nestLockGlobal);
    {

      numHoughPix += numHoughPixPrime;
      rowMax =
          static_cast<int>(std::fmax(static_cast<double>(rowMax), rowMaxPrime));
      rowMin = std::fmin(rowMin, rowMinPrime);
      colMax =
          static_cast<int>(std::fmax(static_cast<double>(colMax), colMaxPrime));
      colMin = std::fmin(colMin, colMinPrime);
    }
    omp_unset_nest_lock(&detectLaneMarkerRidge_nestLockGlobal);
  }
  if (numHoughPix < 1) {
    houghPix.set_size(0, 0);
  } else {
    houghPix.set_size(nonZero.size(0) * nonZero.size(1), 2);
    n = 0U;
    for (numRow = 0; numRow <= numCol; numRow++) {
      ub_loop = tempHoughPixNumsVector[numRow];
      for (j = 0; j < ub_loop; j++) {
        i = static_cast<int>(n + static_cast<unsigned int>(j));
        houghPix[i] = houghPixTemp[j + houghPixTemp.size(0) * numRow];
        houghPix[i + houghPix.size(0)] = numRow + 1;
      }
      n += static_cast<unsigned int>(tempHoughPixNumsVector[numRow]);
    }
    if (static_cast<double>(rowMax) - rowMin >
        static_cast<double>(colMax) - colMin) {
      sortingOrder[0] = 1.0;
      sortingOrder[1] = 2.0;
    } else {
      sortingOrder[0] = 2.0;
      sortingOrder[1] = 1.0;
    }
    for (ub_loop = 0; ub_loop < 2; ub_loop++) {
      for (i = 0; i < numHoughPix; i++) {
        houghPix[i + numHoughPix * ub_loop] =
            houghPix[i + houghPix.size(0) * ub_loop];
      }
    }
    houghPix.set_size(numHoughPix, 2);
    sortrows(houghPix, sortingOrder);
  }
  return numHoughPix;
}

static int getLocationOfMax(const ::coder::array<double, 2U> &A, int &jMax)
{
  double maxarray[180];
  double maxVal;
  double temp;
  double tempmax;
  unsigned int indarray[180];
  int N;
  int iMax;
  int j;
  int tempidx;
  iMax = 0;
  jMax = 0;
  N = A.size(0);
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
    temp, tempidx, tempmax, j)

  for (int i = 0; i < 180; i++) {
    tempmax = A[A.size(0) * i];
    tempidx = 1;
    for (j = 2; j <= N; j++) {
      temp = A[(j + A.size(0) * i) - 1];
      if (temp > tempmax) {
        tempmax = temp;
        tempidx = j;
      }
    }
    maxarray[i] = tempmax;
    indarray[i] = static_cast<unsigned int>(tempidx);
  }
  maxVal = -1.0;
  for (N = 0; N < 180; N++) {
    double d;
    d = maxarray[N];
    if (d > maxVal) {
      maxVal = d;
      jMax = N + 1;
      iMax = static_cast<int>(indarray[N]);
    }
  }
  return iMax;
}

static void houghlines(const ::coder::array<boolean_T, 2U> &varargin_1,
                       const ::coder::array<double, 2U> &varargin_3,
                       const double varargin_4_data[],
                       const int varargin_4_size[2],
                       ::coder::array<struct0_T, 2U> &lines)
{
  ::coder::array<double, 1U> distances2;
  ::coder::array<float, 1U> rhoArray;
  ::coder::array<int, 2U> b_point1Array;
  ::coder::array<int, 2U> houghPix;
  ::coder::array<int, 2U> nonZeroPixelMatrix;
  ::coder::array<int, 2U> point1Array;
  ::coder::array<int, 2U> point2Array;
  ::coder::array<unsigned int, 1U> indices;
  ::coder::array<int, 1U> tempNumsVector;
  ::coder::array<signed char, 1U> thetaArray;
  struct0_T s;
  double firstRho;
  double slope;
  int b_i;
  int i;
  int i1;
  int j;
  int numLines;
  unsigned int tempNum;
  boolean_T useParfor;
  useParfor = (varargin_1.size(0) * varargin_1.size(1) > 2500);
  if (useParfor) {
    findNonZeroOmp(varargin_1, nonZeroPixelMatrix, tempNumsVector);
  } else {
    nonZeroPixelMatrix.set_size(varargin_1.size(0), varargin_1.size(1));
    tempNumsVector.set_size(varargin_1.size(1));
    i = varargin_1.size(1);
    for (j = 0; j < i; j++) {
      tempNum = 0U;
      i1 = varargin_1.size(0);
      for (b_i = 0; b_i < i1; b_i++) {
        if (varargin_1[b_i + varargin_1.size(0) * j]) {
          tempNum++;
          nonZeroPixelMatrix[(static_cast<int>(tempNum) +
                              nonZeroPixelMatrix.size(0) * j) -
                             1] = b_i;
        }
      }
      tempNumsVector[j] = static_cast<int>(tempNum);
    }
  }
  numLines = 0;
  point1Array.set_size(0, 2);
  point2Array.set_size(0, 2);
  thetaArray.set_size(0);
  rhoArray.set_size(0);
  firstRho = varargin_3[0];
  slope = (static_cast<double>(varargin_3.size(1)) - 1.0) /
          (varargin_3[varargin_3.size(1) - 1] - varargin_3[0]);
  i = varargin_4_size[0];
  for (int peakIdx{0}; peakIdx < i; peakIdx++) {
    double thetaVal;
    int colMax;
    int i2;
    int i3;
    int numHoughPix;
    int peak1;
    int peak2;
    int rowMax;
    int rowMax_tmp;
    peak1 = static_cast<int>(varargin_4_data[peakIdx]);
    peak2 = static_cast<int>(varargin_4_data[peakIdx + varargin_4_size[0]]);
    if (useParfor) {
      numHoughPix = getHoughPixelsOmp(nonZeroPixelMatrix, firstRho, slope,
                                      peak1, peak2, tempNumsVector, houghPix);
    } else {
      double colMin;
      double cosTheta;
      double rowMin;
      numHoughPix = -1;
      thetaVal = ((static_cast<double>(peak2) - 1.0) - 90.0) *
                 3.1415926535897931 / 180.0;
      cosTheta = std::cos(thetaVal);
      thetaVal = std::sin(thetaVal);
      rowMax = 0;
      rowMin = rtInf;
      colMax = 0;
      colMin = rtInf;
      i1 = nonZeroPixelMatrix.size(1);
      houghPix.set_size(nonZeroPixelMatrix.size(0) * nonZeroPixelMatrix.size(1),
                        2);
      for (j = 0; j < i1; j++) {
        i2 = tempNumsVector[j];
        for (b_i = 0; b_i < i2; b_i++) {
          i3 = nonZeroPixelMatrix[b_i + nonZeroPixelMatrix.size(0) * j];
          if (static_cast<int>(
                  (slope * ((((static_cast<double>(j) + 1.0) - 1.0) * cosTheta +
                             static_cast<double>(i3) * thetaVal) -
                            firstRho) +
                   1.0) +
                  0.5) == peak1) {
            numHoughPix++;
            houghPix[numHoughPix] = i3 + 1;
            houghPix[numHoughPix + houghPix.size(0)] = j + 1;
            rowMax_tmp = houghPix[numHoughPix];
            rowMax = static_cast<int>(std::fmax(
                static_cast<double>(rowMax), static_cast<double>(rowMax_tmp)));
            rowMin = std::fmin(rowMin, static_cast<double>(rowMax_tmp));
            colMax = static_cast<int>(std::fmax(static_cast<double>(colMax),
                                                static_cast<double>(j) + 1.0));
            colMin = std::fmin(colMin, static_cast<double>(j) + 1.0);
          }
        }
      }
      if (numHoughPix + 1 < 1) {
        houghPix.set_size(0, 0);
      } else {
        double sortingOrder[2];
        if (static_cast<double>(rowMax) - rowMin >
            static_cast<double>(colMax) - colMin) {
          sortingOrder[0] = 1.0;
          sortingOrder[1] = 2.0;
        } else {
          sortingOrder[0] = 2.0;
          sortingOrder[1] = 1.0;
        }
        for (i1 = 0; i1 < 2; i1++) {
          for (i2 = 0; i2 <= numHoughPix; i2++) {
            houghPix[i2 + (numHoughPix + 1) * i1] =
                houghPix[i2 + houghPix.size(0) * i1];
          }
        }
        houghPix.set_size(numHoughPix + 1, 2);
        sortrows(houghPix, sortingOrder);
      }
      numHoughPix++;
    }
    if (numHoughPix >= 1) {
      int exitg1;
      distances2.set_size(houghPix.size(0) - 1);
      thetaVal = 0.0;
      i1 = houghPix.size(0);
      for (int k{0}; k <= i1 - 2; k++) {
        rowMax_tmp = houghPix[k + 1] - houghPix[k];
        numHoughPix = 1;
        rowMax = 2;
        do {
          exitg1 = 0;
          if ((rowMax & 1) != 0) {
            numHoughPix *= rowMax_tmp;
          }
          rowMax >>= 1;
          if (rowMax == 0) {
            exitg1 = 1;
          } else {
            rowMax_tmp *= rowMax_tmp;
          }
        } while (exitg1 == 0);
        rowMax_tmp = houghPix[(k + houghPix.size(0)) + 1] -
                     houghPix[k + houghPix.size(0)];
        colMax = 1;
        rowMax = 2;
        do {
          exitg1 = 0;
          if ((rowMax & 1) != 0) {
            colMax *= rowMax_tmp;
          }
          rowMax >>= 1;
          if (rowMax == 0) {
            exitg1 = 1;
          } else {
            rowMax_tmp *= rowMax_tmp;
          }
        } while (exitg1 == 0);
        i2 = numHoughPix + colMax;
        distances2[k] = i2;
        if (i2 > 3600) {
          thetaVal++;
        }
      }
      indices.set_size(static_cast<int>(thetaVal + 2.0));
      indices[0] = 0U;
      indices[static_cast<int>(thetaVal + 2.0) - 1] =
          static_cast<unsigned int>(houghPix.size(0));
      tempNum = 1U;
      i1 = houghPix.size(0);
      for (int k{0}; k <= i1 - 2; k++) {
        if (distances2[k] > 3600.0) {
          tempNum++;
          indices[static_cast<int>(tempNum) - 1] =
              static_cast<unsigned int>(k + 1);
        }
      }
      i1 = indices.size(0);
      for (int k{0}; k <= i1 - 2; k++) {
        int a_tmp;
        int b_a_tmp;
        a_tmp = static_cast<int>(indices[k + 1]) - 1;
        j = houghPix[static_cast<int>(indices[k])];
        b_a_tmp = houghPix[a_tmp];
        rowMax_tmp = b_a_tmp - j;
        numHoughPix = 1;
        rowMax = 2;
        do {
          exitg1 = 0;
          if ((rowMax & 1) != 0) {
            numHoughPix *= rowMax_tmp;
          }
          rowMax >>= 1;
          if (rowMax == 0) {
            exitg1 = 1;
          } else {
            rowMax_tmp *= rowMax_tmp;
          }
        } while (exitg1 == 0);
        b_i = houghPix[static_cast<int>(indices[k]) + houghPix.size(0)];
        a_tmp = houghPix[a_tmp + houghPix.size(0)];
        rowMax_tmp = a_tmp - b_i;
        colMax = 1;
        rowMax = 2;
        do {
          exitg1 = 0;
          if ((rowMax & 1) != 0) {
            colMax *= rowMax_tmp;
          }
          rowMax >>= 1;
          if (rowMax == 0) {
            exitg1 = 1;
          } else {
            rowMax_tmp *= rowMax_tmp;
          }
        } while (exitg1 == 0);
        if (numHoughPix + colMax >= 100) {
          numLines++;
          numHoughPix = point1Array.size(0);
          b_point1Array.set_size(point1Array.size(0) + 1, 2);
          for (i2 = 0; i2 < 2; i2++) {
            for (i3 = 0; i3 < numHoughPix; i3++) {
              b_point1Array[i3 + b_point1Array.size(0) * i2] =
                  point1Array[i3 + point1Array.size(0) * i2];
            }
          }
          b_point1Array[point1Array.size(0)] = b_i;
          b_point1Array[point1Array.size(0) + b_point1Array.size(0)] = j;
          point1Array.set_size(b_point1Array.size(0), 2);
          rowMax_tmp = b_point1Array.size(0) << 1;
          for (i2 = 0; i2 < rowMax_tmp; i2++) {
            point1Array[i2] = b_point1Array[i2];
          }
          numHoughPix = point2Array.size(0);
          b_point1Array.set_size(point2Array.size(0) + 1, 2);
          for (i2 = 0; i2 < 2; i2++) {
            for (i3 = 0; i3 < numHoughPix; i3++) {
              b_point1Array[i3 + b_point1Array.size(0) * i2] =
                  point2Array[i3 + point2Array.size(0) * i2];
            }
          }
          b_point1Array[point2Array.size(0)] = a_tmp;
          b_point1Array[point2Array.size(0) + b_point1Array.size(0)] = b_a_tmp;
          point2Array.set_size(b_point1Array.size(0), 2);
          for (i2 = 0; i2 < rowMax_tmp; i2++) {
            point2Array[i2] = b_point1Array[i2];
          }
          i2 = thetaArray.size(0);
          thetaArray.set_size(thetaArray.size(0) + 1);
          thetaArray[i2] = static_cast<signed char>(peak2 - 91);
          i2 = rhoArray.size(0);
          rhoArray.set_size(rhoArray.size(0) + 1);
          rhoArray[i2] = static_cast<float>(varargin_3[peak1 - 1]);
        }
      }
    }
  }
  s.point1[0] = 0.0;
  s.point2[0] = 0.0;
  s.point1[1] = 0.0;
  s.point2[1] = 0.0;
  s.theta = 0.0;
  s.rho = 0.0;
  lines.set_size(1, numLines);
  for (int k{0}; k < numLines; k++) {
    lines[k] = s;
    lines[k].point1[0] = point1Array[k];
    lines[k].point2[0] = point2Array[k];
    lines[k].point1[1] = point1Array[k + point1Array.size(0)];
    lines[k].point2[1] = point2Array[k + point2Array.size(0)];
    lines[k].theta = thetaArray[k];
    lines[k].rho = rhoArray[k];
  }
}

static void houghpeaks(const ::coder::array<double, 2U> &varargin_1,
                       double varargin_4, double peaks_data[],
                       int peaks_size[2])
{
  ::coder::array<double, 2U> hnew;
  int nhoodSizeDefault_idx_0;
  nhoodSizeDefault_idx_0 =
      2 * static_cast<int>(
              std::ceil(static_cast<double>(varargin_1.size(0)) / 50.0 / 2.0)) +
      1;
  if (varargin_1.size(0) < 3) {
    nhoodSizeDefault_idx_0 = 1;
  }
  if (varargin_1.size(0) * 180 < 1) {
    peaks_size[0] = 0;
    peaks_size[1] = 0;
  } else {
    int peakCoordinates[10];
    int jPeak;
    int loop_ub_tmp;
    int nhoodCenter_i;
    int numRowH;
    int peakIdx;
    int thetaToRemove;
    boolean_T isDone;
    boolean_T isThetaAntisymmetric;
    numRowH = varargin_1.size(0);
    isDone = false;
    hnew.set_size(varargin_1.size(0), 180);
    loop_ub_tmp = varargin_1.size(0) * 180;
    for (jPeak = 0; jPeak < loop_ub_tmp; jPeak++) {
      hnew[jPeak] = varargin_1[jPeak];
    }
    nhoodCenter_i = static_cast<int>(
        (static_cast<double>(nhoodSizeDefault_idx_0) - 1.0) / 2.0);
    peakIdx = 0;
    nhoodSizeDefault_idx_0 = -90;
    jPeak = -90;
    for (thetaToRemove = 0; thetaToRemove < 179; thetaToRemove++) {
      if (nhoodSizeDefault_idx_0 >
          static_cast<signed char>(static_cast<signed char>(thetaToRemove + 1) -
                                   90)) {
        nhoodSizeDefault_idx_0 = static_cast<signed char>(
            static_cast<signed char>(thetaToRemove + 1) - 90);
      }
      if (jPeak < static_cast<signed char>(
                      static_cast<signed char>(thetaToRemove + 1) - 90)) {
        jPeak = static_cast<signed char>(
            static_cast<signed char>(thetaToRemove + 1) - 90);
      }
    }
    isThetaAntisymmetric =
        (std::abs(
             static_cast<double>(nhoodSizeDefault_idx_0) +
             std::abs(static_cast<double>(jPeak - nhoodSizeDefault_idx_0)) /
                 179.0 * 5.0) <= jPeak);
    while (!isDone) {
      if (numRowH < 1200) {
        if (hnew.size(0) * 180 <= 2) {
          if (hnew.size(0) * 180 == 1) {
            nhoodSizeDefault_idx_0 = 1;
          } else {
            double d;
            d = hnew[hnew.size(0) * 180 - 1];
            if ((hnew[0] < d) || (std::isnan(hnew[0]) && (!std::isnan(d)))) {
              nhoodSizeDefault_idx_0 = loop_ub_tmp;
            } else {
              nhoodSizeDefault_idx_0 = 1;
            }
          }
        } else {
          if (!std::isnan(hnew[0])) {
            nhoodSizeDefault_idx_0 = 1;
          } else {
            boolean_T exitg1;
            nhoodSizeDefault_idx_0 = 0;
            thetaToRemove = 2;
            exitg1 = false;
            while ((!exitg1) && (thetaToRemove <= loop_ub_tmp)) {
              if (!std::isnan(hnew[thetaToRemove - 1])) {
                nhoodSizeDefault_idx_0 = thetaToRemove;
                exitg1 = true;
              } else {
                thetaToRemove++;
              }
            }
          }
          if (nhoodSizeDefault_idx_0 == 0) {
            nhoodSizeDefault_idx_0 = 1;
          } else {
            double ex;
            ex = hnew[nhoodSizeDefault_idx_0 - 1];
            jPeak = nhoodSizeDefault_idx_0 + 1;
            for (thetaToRemove = jPeak; thetaToRemove <= loop_ub_tmp;
                 thetaToRemove++) {
              double d;
              d = hnew[thetaToRemove - 1];
              if (ex < d) {
                ex = d;
                nhoodSizeDefault_idx_0 = thetaToRemove;
              }
            }
          }
        }
        if (static_cast<unsigned int>(hnew.size(0)) == 0U) {
          jPeak = MAX_int32_T;
        } else {
          jPeak = static_cast<int>(
              static_cast<unsigned int>(nhoodSizeDefault_idx_0 - 1) /
              static_cast<unsigned int>(hnew.size(0)));
        }
        nhoodSizeDefault_idx_0 -= jPeak * hnew.size(0);
        jPeak++;
      } else {
        nhoodSizeDefault_idx_0 = getLocationOfMax(hnew, jPeak);
      }
      if (hnew[(nhoodSizeDefault_idx_0 + hnew.size(0) * (jPeak - 1)) - 1] >=
          varargin_4) {
        int rhoMax;
        int rhoMin;
        int thetaMin;
        peakIdx++;
        peakCoordinates[peakIdx - 1] = nhoodSizeDefault_idx_0;
        peakCoordinates[peakIdx + 4] = jPeak;
        rhoMin = nhoodSizeDefault_idx_0 - nhoodCenter_i;
        rhoMax = nhoodSizeDefault_idx_0 + nhoodCenter_i;
        thetaMin = jPeak - 2;
        nhoodSizeDefault_idx_0 = jPeak + 2;
        if (rhoMin < 1) {
          rhoMin = 1;
        }
        if (rhoMax > numRowH) {
          rhoMax = numRowH;
        }
        for (int theta{thetaMin}; theta <= nhoodSizeDefault_idx_0; theta++) {
          for (int rho{rhoMin}; rho <= rhoMax; rho++) {
            jPeak = rho - 1;
            thetaToRemove = theta;
            if (isThetaAntisymmetric) {
              if (theta > 180) {
                jPeak = numRowH - rho;
                thetaToRemove = theta - 180;
              } else if (theta < 1) {
                jPeak = numRowH - rho;
                thetaToRemove = theta + 180;
              }
            }
            if ((thetaToRemove <= 180) && (thetaToRemove >= 1)) {
              hnew[jPeak + hnew.size(0) * (thetaToRemove - 1)] = 0.0;
            }
          }
        }
        isDone = (peakIdx == 5);
      } else {
        isDone = true;
      }
    }
    if (peakIdx == 0) {
      peaks_size[0] = 0;
      peaks_size[1] = 0;
    } else {
      peaks_size[0] = peakIdx;
      peaks_size[1] = 2;
      for (jPeak = 0; jPeak < 2; jPeak++) {
        for (nhoodSizeDefault_idx_0 = 0; nhoodSizeDefault_idx_0 < peakIdx;
             nhoodSizeDefault_idx_0++) {
          peaks_data[nhoodSizeDefault_idx_0 + peakIdx * jPeak] =
              peakCoordinates[nhoodSizeDefault_idx_0 + 5 * jPeak];
        }
      }
    }
  }
}

static double labelingWu_parallel(const ::coder::array<boolean_T, 2U> &im,
                                  double M, double N,
                                  ::coder::array<double, 2U> &L)
{
  ::coder::array<double, 1U> P;
  ::coder::array<int, 1U> chunksSizeAndLabels;
  double b_c;
  double c_tmp;
  double d;
  double firstLabel;
  double k;
  double label;
  double rootj;
  double stripeWidth;
  int c;
  int c_c;
  int exitg1;
  int i;
  int i1;
  int nParallelStripes;
  int r;
  nParallelStripes = static_cast<int>(M);
  i = static_cast<int>(M);
  L.set_size(nParallelStripes, static_cast<int>(N));
  chunksSizeAndLabels.set_size(static_cast<int>(N + 8.0));
  nParallelStripes =
      static_cast<int>(std::fmax(1.0, std::fmin(std::floor(N / 4.0), 8.0)));
  stripeWidth = std::ceil(N / static_cast<double>(nParallelStripes));
  P.set_size(static_cast<int>(std::ceil(M * N / 2.0) + 1.0));
  P[0] = 0.0;
  nParallelStripes--;
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
    rootj, b_c, firstLabel, label, c_tmp, d, i1, c_c, r, exitg1)

  for (int thread = 0; thread <= nParallelStripes; thread++) {
    c_tmp = static_cast<double>(thread) * stripeWidth + 1.0;
    d = (static_cast<double>(thread) + 1.0) * stripeWidth;
    rootj = std::round(std::fmin(d + 1.0, N + 1.0));
    if (rootj < 2.147483648E+9) {
      if (rootj >= -2.147483648E+9) {
        i1 = static_cast<int>(rootj);
      } else {
        i1 = MIN_int32_T;
      }
    } else if (rootj >= 2.147483648E+9) {
      i1 = MAX_int32_T;
    } else {
      i1 = 0;
    }
    chunksSizeAndLabels[static_cast<int>(c_tmp) - 1] = i1;
    label = std::floor(c_tmp / 2.0) * std::floor((M + 1.0) / 2.0) + 1.0;
    firstLabel = label;
    i1 = static_cast<int>(std::fmin(d, N) + (1.0 - c_tmp));
    for (c_c = 0; c_c < i1; c_c++) {
      b_c = c_tmp + static_cast<double>(c_c);
      for (r = 0; r < i; r++) {
        if (im[r + im.size(0) * (static_cast<int>(b_c) - 1)]) {
          if ((b_c > c_tmp) &&
              im[r + im.size(0) * (static_cast<int>(b_c) - 2)]) {
            L[r + L.size(0) * (static_cast<int>(b_c) - 1)] =
                L[r + L.size(0) * (static_cast<int>(b_c) - 2)];
          } else if ((static_cast<double>(r) + 1.0 < M) && (b_c > c_tmp) &&
                     im[(r + im.size(0) * (static_cast<int>(b_c) - 2)) + 1]) {
            if ((b_c > c_tmp) && (static_cast<unsigned int>(r) + 1U > 1U) &&
                im[(r + im.size(0) * (static_cast<int>(b_c) - 2)) - 1]) {
              L[r + L.size(0) * (static_cast<int>(b_c) - 1)] =
                  L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) - 1];
              do {
                exitg1 = 0;
                d = L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
                if (P[static_cast<int>(d + 1.0) - 1] < d) {
                  L[r + L.size(0) * (static_cast<int>(b_c) - 1)] =
                      P[static_cast<int>(d + 1.0) - 1];
                } else {
                  exitg1 = 1;
                }
              } while (exitg1 == 0);
              rootj = L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) + 1];
              if (L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) - 1] !=
                  rootj) {
                while (P[static_cast<int>(rootj + 1.0) - 1] < rootj) {
                  rootj = P[static_cast<int>(rootj + 1.0) - 1];
                }
                if (d > rootj) {
                  L[r + L.size(0) * (static_cast<int>(b_c) - 1)] = rootj;
                }
                do {
                  exitg1 = 0;
                  d = L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) + 1];
                  rootj = P[static_cast<int>(d + 1.0) - 1];
                  if (rootj < d) {
                    P[static_cast<int>(d + 1.0) - 1] =
                        L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
                    L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) + 1] =
                        rootj;
                  } else {
                    exitg1 = 1;
                  }
                } while (exitg1 == 0);
                P[static_cast<int>(d + 1.0) - 1] =
                    L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
              }
              do {
                exitg1 = 0;
                d = L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) - 1];
                rootj = P[static_cast<int>(d + 1.0) - 1];
                if (rootj < d) {
                  P[static_cast<int>(d + 1.0) - 1] =
                      L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
                  L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) - 1] = rootj;
                } else {
                  exitg1 = 1;
                }
              } while (exitg1 == 0);
              P[static_cast<int>(d + 1.0) - 1] =
                  L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
            } else if ((static_cast<unsigned int>(r) + 1U > 1U) &&
                       im[(r + im.size(0) * (static_cast<int>(b_c) - 1)) - 1]) {
              L[r + L.size(0) * (static_cast<int>(b_c) - 1)] =
                  L[(r + L.size(0) * (static_cast<int>(b_c) - 1)) - 1];
              do {
                exitg1 = 0;
                d = L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
                if (P[static_cast<int>(d + 1.0) - 1] < d) {
                  L[r + L.size(0) * (static_cast<int>(b_c) - 1)] =
                      P[static_cast<int>(d + 1.0) - 1];
                } else {
                  exitg1 = 1;
                }
              } while (exitg1 == 0);
              rootj = L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) + 1];
              if (L[(r + L.size(0) * (static_cast<int>(b_c) - 1)) - 1] !=
                  rootj) {
                while (P[static_cast<int>(rootj + 1.0) - 1] < rootj) {
                  rootj = P[static_cast<int>(rootj + 1.0) - 1];
                }
                if (d > rootj) {
                  L[r + L.size(0) * (static_cast<int>(b_c) - 1)] = rootj;
                }
                do {
                  exitg1 = 0;
                  d = L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) + 1];
                  rootj = P[static_cast<int>(d + 1.0) - 1];
                  if (rootj < d) {
                    P[static_cast<int>(d + 1.0) - 1] =
                        L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
                    L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) + 1] =
                        rootj;
                  } else {
                    exitg1 = 1;
                  }
                } while (exitg1 == 0);
                P[static_cast<int>(d + 1.0) - 1] =
                    L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
              }
              do {
                exitg1 = 0;
                d = L[(r + L.size(0) * (static_cast<int>(b_c) - 1)) - 1];
                rootj = P[static_cast<int>(d + 1.0) - 1];
                if (rootj < d) {
                  P[static_cast<int>(d + 1.0) - 1] =
                      L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
                  L[(r + L.size(0) * (static_cast<int>(b_c) - 1)) - 1] = rootj;
                } else {
                  exitg1 = 1;
                }
              } while (exitg1 == 0);
              P[static_cast<int>(d + 1.0) - 1] =
                  L[r + L.size(0) * (static_cast<int>(b_c) - 1)];
            } else {
              L[r + L.size(0) * (static_cast<int>(b_c) - 1)] =
                  L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) + 1];
            }
          } else if ((b_c > c_tmp) &&
                     (static_cast<unsigned int>(r) + 1U > 1U) &&
                     im[(r + im.size(0) * (static_cast<int>(b_c) - 2)) - 1]) {
            L[r + L.size(0) * (static_cast<int>(b_c) - 1)] =
                L[(r + L.size(0) * (static_cast<int>(b_c) - 2)) - 1];
          } else if ((static_cast<unsigned int>(r) + 1U > 1U) &&
                     im[(r + im.size(0) * (static_cast<int>(b_c) - 1)) - 1]) {
            L[r + L.size(0) * (static_cast<int>(b_c) - 1)] =
                L[(r + L.size(0) * (static_cast<int>(b_c) - 1)) - 1];
          } else {
            L[r + L.size(0) * (static_cast<int>(b_c) - 1)] = label;
            P[static_cast<int>(label + 1.0) - 1] = label;
            label++;
          }
        } else {
          L[r + L.size(0) * (static_cast<int>(b_c) - 1)] = 0.0;
        }
      }
    }
    d = label - firstLabel;
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        i1 = static_cast<int>(d);
      } else {
        i1 = MIN_int32_T;
      }
    } else if (d >= 2.147483648E+9) {
      i1 = MAX_int32_T;
    } else {
      i1 = 0;
    }
    chunksSizeAndLabels[static_cast<int>(c_tmp + 1.0) - 1] = i1;
  }
  for (c = chunksSizeAndLabels[0] - 1; c + 1 <= N;
       c = chunksSizeAndLabels[c] - 1) {
    for (nParallelStripes = 0; nParallelStripes < i; nParallelStripes++) {
      if (L[nParallelStripes + L.size(0) * c] != 0.0) {
        double b_i;
        double j;
        double root;
        if (static_cast<unsigned int>(nParallelStripes) + 1U > 1U) {
          b_i = L[(nParallelStripes + L.size(0) * (c - 1)) - 1];
          if (b_i != 0.0) {
            j = L[nParallelStripes + L.size(0) * c];
            root = b_i;
            while (P[static_cast<int>(root + 1.0) - 1] < root) {
              root = P[static_cast<int>(root + 1.0) - 1];
            }
            if (b_i != j) {
              stripeWidth = j;
              while (P[static_cast<int>(stripeWidth + 1.0) - 1] < stripeWidth) {
                stripeWidth = P[static_cast<int>(stripeWidth + 1.0) - 1];
              }
              if (root > stripeWidth) {
                root = stripeWidth;
              }
              do {
                exitg1 = 0;
                stripeWidth = P[static_cast<int>(j + 1.0) - 1];
                if (stripeWidth < j) {
                  P[static_cast<int>(j + 1.0) - 1] = root;
                  j = stripeWidth;
                } else {
                  exitg1 = 1;
                }
              } while (exitg1 == 0);
              P[static_cast<int>(j + 1.0) - 1] = root;
            }
            do {
              exitg1 = 0;
              stripeWidth = P[static_cast<int>(b_i + 1.0) - 1];
              if (stripeWidth < b_i) {
                P[static_cast<int>(b_i + 1.0) - 1] = root;
                b_i = stripeWidth;
              } else {
                exitg1 = 1;
              }
            } while (exitg1 == 0);
            P[static_cast<int>(b_i + 1.0) - 1] = root;
            L[nParallelStripes + L.size(0) * c] = root;
          }
        }
        if (static_cast<double>(nParallelStripes) + 1.0 < M) {
          b_i = L[(nParallelStripes + L.size(0) * (c - 1)) + 1];
          if (b_i != 0.0) {
            j = L[nParallelStripes + L.size(0) * c];
            root = b_i;
            while (P[static_cast<int>(root + 1.0) - 1] < root) {
              root = P[static_cast<int>(root + 1.0) - 1];
            }
            if (b_i != j) {
              stripeWidth = j;
              while (P[static_cast<int>(stripeWidth + 1.0) - 1] < stripeWidth) {
                stripeWidth = P[static_cast<int>(stripeWidth + 1.0) - 1];
              }
              if (root > stripeWidth) {
                root = stripeWidth;
              }
              do {
                exitg1 = 0;
                stripeWidth = P[static_cast<int>(j + 1.0) - 1];
                if (stripeWidth < j) {
                  P[static_cast<int>(j + 1.0) - 1] = root;
                  j = stripeWidth;
                } else {
                  exitg1 = 1;
                }
              } while (exitg1 == 0);
              P[static_cast<int>(j + 1.0) - 1] = root;
            }
            do {
              exitg1 = 0;
              stripeWidth = P[static_cast<int>(b_i + 1.0) - 1];
              if (stripeWidth < b_i) {
                P[static_cast<int>(b_i + 1.0) - 1] = root;
                b_i = stripeWidth;
              } else {
                exitg1 = 1;
              }
            } while (exitg1 == 0);
            P[static_cast<int>(b_i + 1.0) - 1] = root;
            L[nParallelStripes + L.size(0) * c] = root;
          }
        }
        b_i = L[nParallelStripes + L.size(0) * (c - 1)];
        if (b_i != 0.0) {
          j = L[nParallelStripes + L.size(0) * c];
          root = b_i;
          while (P[static_cast<int>(root + 1.0) - 1] < root) {
            root = P[static_cast<int>(root + 1.0) - 1];
          }
          if (b_i != j) {
            stripeWidth = j;
            while (P[static_cast<int>(stripeWidth + 1.0) - 1] < stripeWidth) {
              stripeWidth = P[static_cast<int>(stripeWidth + 1.0) - 1];
            }
            if (root > stripeWidth) {
              root = stripeWidth;
            }
            do {
              exitg1 = 0;
              stripeWidth = P[static_cast<int>(j + 1.0) - 1];
              if (stripeWidth < j) {
                P[static_cast<int>(j + 1.0) - 1] = root;
                j = stripeWidth;
              } else {
                exitg1 = 1;
              }
            } while (exitg1 == 0);
            P[static_cast<int>(j + 1.0) - 1] = root;
          }
          do {
            exitg1 = 0;
            stripeWidth = P[static_cast<int>(b_i + 1.0) - 1];
            if (stripeWidth < b_i) {
              P[static_cast<int>(b_i + 1.0) - 1] = root;
              b_i = stripeWidth;
            } else {
              exitg1 = 1;
            }
          } while (exitg1 == 0);
          P[static_cast<int>(b_i + 1.0) - 1] = root;
          L[nParallelStripes + L.size(0) * c] = root;
        }
      }
    }
  }
  k = 1.0;
  c = 1;
  while (c <= N) {
    int b_qY;
    int qY;
    if (c < -2147483647) {
      qY = MIN_int32_T;
    } else {
      qY = c - 1;
    }
    stripeWidth =
        std::round(static_cast<double>(qY) / 2.0) * std::floor((M + 1.0) / 2.0);
    if (stripeWidth < 2.147483648E+9) {
      if (stripeWidth >= -2.147483648E+9) {
        nParallelStripes = static_cast<int>(stripeWidth);
      } else {
        nParallelStripes = MIN_int32_T;
      }
    } else if (stripeWidth >= 2.147483648E+9) {
      nParallelStripes = MAX_int32_T;
    } else {
      nParallelStripes = 0;
    }
    if (nParallelStripes > 2147483646) {
      qY = MAX_int32_T;
    } else {
      qY = nParallelStripes + 1;
    }
    if (qY > 2147483646) {
      b_qY = MAX_int32_T;
    } else {
      b_qY = qY + 1;
    }
    if (c > 2147483646) {
      nParallelStripes = MAX_int32_T;
    } else {
      nParallelStripes = c + 1;
    }
    nParallelStripes = chunksSizeAndLabels[nParallelStripes - 1];
    if ((qY < 0) && (nParallelStripes < MIN_int32_T - qY)) {
      qY = MIN_int32_T;
    } else if ((qY > 0) && (nParallelStripes > MAX_int32_T - qY)) {
      qY = MAX_int32_T;
    } else {
      qY += nParallelStripes;
    }
    for (nParallelStripes = b_qY; nParallelStripes <= qY; nParallelStripes++) {
      stripeWidth = P[nParallelStripes - 1];
      if (stripeWidth < static_cast<double>(nParallelStripes) - 1.0) {
        P[nParallelStripes - 1] = P[static_cast<int>(stripeWidth + 1.0) - 1];
      } else {
        P[nParallelStripes - 1] = k;
        k++;
      }
    }
    c = chunksSizeAndLabels[c - 1];
  }
  k--;
  nParallelStripes = static_cast<int>(N) - 1;
#pragma omp parallel for num_threads(omp_get_max_threads()) private(r)

  for (c_c = 0; c_c <= nParallelStripes; c_c++) {
    for (r = 0; r < i; r++) {
      L[r + L.size(0) * c_c] =
          P[static_cast<int>(L[r + L.size(0) * c_c] + 1.0) - 1];
    }
  }
  return k;
}

static void sortrows(::coder::array<int, 2U> &y, const double varargin_1[2])
{
  ::coder::array<int, 1U> idx;
  ::coder::array<int, 1U> iwork;
  int col[2];
  int b_i;
  int i;
  int j;
  int k;
  int loop_ub;
  int n;
  int qEnd;
  int v1;
  int v2;
  boolean_T exitg1;
  boolean_T p;
  col[0] = static_cast<int>(varargin_1[0]);
  col[1] = static_cast<int>(varargin_1[1]);
  n = y.size(0) + 1;
  idx.set_size(y.size(0));
  loop_ub = y.size(0);
  for (i = 0; i < loop_ub; i++) {
    idx[i] = 0;
  }
  iwork.set_size(y.size(0));
  i = y.size(0) - 1;
  for (k = 1; k <= i; k += 2) {
    p = true;
    loop_ub = 0;
    exitg1 = false;
    while ((!exitg1) && (loop_ub < 2)) {
      v1 = y[(k + y.size(0) * (col[loop_ub] - 1)) - 1];
      v2 = y[k + y.size(0) * (col[loop_ub] - 1)];
      if (v1 == v2) {
        loop_ub++;
      } else {
        p = (v1 <= v2);
        exitg1 = true;
      }
    }
    if (p) {
      idx[k - 1] = k;
      idx[k] = k + 1;
    } else {
      idx[k - 1] = k + 1;
      idx[k] = k;
    }
  }
  if ((y.size(0) & 1) != 0) {
    idx[y.size(0) - 1] = y.size(0);
  }
  b_i = 2;
  while (b_i < n - 1) {
    int i2;
    i2 = b_i << 1;
    j = 1;
    for (int pEnd{b_i + 1}; pEnd < n; pEnd = qEnd + b_i) {
      int b_p;
      int kEnd;
      int q;
      b_p = j;
      q = pEnd;
      qEnd = j + i2;
      if (qEnd > n) {
        qEnd = n;
      }
      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        p = true;
        loop_ub = 0;
        exitg1 = false;
        while ((!exitg1) && (loop_ub < 2)) {
          v1 = y[(idx[b_p - 1] + y.size(0) * (col[loop_ub] - 1)) - 1];
          v2 = y[(idx[q - 1] + y.size(0) * (col[loop_ub] - 1)) - 1];
          if (v1 == v2) {
            loop_ub++;
          } else {
            p = (v1 <= v2);
            exitg1 = true;
          }
        }
        if (p) {
          iwork[k] = idx[b_p - 1];
          b_p++;
          if (b_p == pEnd) {
            while (q < qEnd) {
              k++;
              iwork[k] = idx[q - 1];
              q++;
            }
          }
        } else {
          iwork[k] = idx[q - 1];
          q++;
          if (q == qEnd) {
            while (b_p < pEnd) {
              k++;
              iwork[k] = idx[b_p - 1];
              b_p++;
            }
          }
        }
        k++;
      }
      for (k = 0; k < kEnd; k++) {
        idx[(j + k) - 1] = iwork[k];
      }
      j = qEnd;
    }
    b_i = i2;
  }
  iwork.set_size(y.size(0));
  for (j = 0; j < 2; j++) {
    for (b_i = 0; b_i <= i; b_i++) {
      iwork[b_i] = y[(idx[b_i] + y.size(0) * j) - 1];
    }
    for (b_i = 0; b_i <= i; b_i++) {
      y[b_i + y.size(0) * j] = iwork[b_i];
    }
  }
}

static void
standardHoughTransformOptimized(const ::coder::array<boolean_T, 2U> &BW,
                                const ::coder::array<double, 2U> &rho,
                                double numRow, double numCol,
                                ::coder::array<double, 2U> &H)
{
  ::coder::array<int, 2U> nonZeroPixelMatrix;
  ::coder::array<int, 1U> numNonZeros;
  ::coder::array<int, 1U> rhoIdxVector;
  double cost[180];
  double sint[180];
  double firstRho;
  double myRho;
  double slope;
  int b_i;
  int c_i;
  int i1;
  int j;
  int n;
  int rhoLength;
  rhoLength = rho.size(1);
  firstRho = rho[0];
  H.set_size(rho.size(1), 180);
  for (int i{0}; i < 180; i++) {
    slope = (((static_cast<double>(i) + 1.0) - 1.0) - 90.0) *
            3.1415926535897931 / 180.0;
    cost[i] = std::cos(slope);
    sint[i] = std::sin(slope);
  }
  slope = (static_cast<double>(rho.size(1)) - 1.0) /
          (rho[rho.size(1) - 1] - rho[0]);
  findNonZero(BW, numRow, numCol, nonZeroPixelMatrix, numNonZeros);
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
    rhoIdxVector, myRho, n, b_i, j, i1, c_i)

  for (int thetaIdx = 0; thetaIdx < 180; thetaIdx++) {
    rhoIdxVector.set_size(rhoLength);
    for (b_i = 0; b_i < rhoLength; b_i++) {
      rhoIdxVector[b_i] = 0;
    }
    b_i = static_cast<int>(numCol);
    for (j = 0; j < b_i; j++) {
      i1 = numNonZeros[j];
      for (c_i = 0; c_i < i1; c_i++) {
        n = nonZeroPixelMatrix[c_i + nonZeroPixelMatrix.size(0) * j];
        myRho = ((static_cast<double>(j) + 1.0) - 1.0) * cost[thetaIdx] +
                (static_cast<double>(n) - 1.0) * sint[thetaIdx];
        n = static_cast<int>(slope * (myRho - firstRho) + 0.5);
        rhoIdxVector[n] = rhoIdxVector[n] + 1;
      }
    }
    n = rhoIdxVector.size(0);
    for (b_i = 0; b_i < n; b_i++) {
      H[b_i + H.size(0) * thetaIdx] = rhoIdxVector[b_i];
    }
  }
}

} // namespace coder
void detectLaneMarkerRidge(const ::coder::array<unsigned char, 2U> &bevImage,
                           double approxLaneWidthPixels,
                           const ::coder::array<unsigned char, 2U> &mask,
                           double sensitivity,
                           ::coder::array<struct0_T, 2U> &lines)
{
  ::coder::array<double, 2U> H;
  ::coder::array<double, 2U> R;
  ::coder::array<double, 2U> r1;
  ::coder::array<float, 2U> Ileft;
  ::coder::array<int, 2U> idxA;
  ::coder::array<unsigned int, 2U> idxDir;
  ::coder::array<unsigned int, 2U> y;
  ::coder::array<short, 2U> x;
  ::coder::array<unsigned char, 2U> Ipad;
  ::coder::array<unsigned char, 2U> r;
  ::coder::array<boolean_T, 2U> BW;
  double counts[256];
  double tmp_data[10];
  double q;
  double sizeB_idx_1;
  int ni[256];
  int tmp_size[2];
  int high_i;
  int i;
  int k;
  int mid_i;
  int nx;
  int nx_tmp;
  unsigned int sizeA_idx_1;
  int small;
  short ii_data;
  boolean_T exitg1;
  if (!isInitialized_detectLaneMarkerRidge) {
    detectLaneMarkerRidge_initialize();
  }
  if ((bevImage.size(0) == 0) || (bevImage.size(1) == 0)) {
    sizeB_idx_1 = static_cast<double>(bevImage.size(1)) +
                  2.0 * (approxLaneWidthPixels + 1.0);
    Ipad.set_size(bevImage.size(0), static_cast<int>(sizeB_idx_1));
    small = bevImage.size(0) * static_cast<int>(sizeB_idx_1);
    for (i = 0; i < small; i++) {
      Ipad[i] = 0U;
    }
  } else {
    sizeA_idx_1 = static_cast<unsigned int>(bevImage.size(1));
    sizeB_idx_1 = 2.0 * (approxLaneWidthPixels + 1.0) +
                  static_cast<double>(bevImage.size(1));
    if (!(bevImage.size(0) < sizeB_idx_1)) {
      sizeB_idx_1 = bevImage.size(0);
    }
    idxA.set_size(static_cast<int>(sizeB_idx_1), 2);
    y.set_size(1, bevImage.size(0));
    small = bevImage.size(0) - 1;
    for (i = 0; i <= small; i++) {
      y[i] = static_cast<unsigned int>(i) + 1U;
    }
    idxDir.set_size(1, y.size(1));
    small = y.size(1);
    for (i = 0; i < small; i++) {
      idxDir[i] = y[i];
    }
    small = idxDir.size(1);
    for (i = 0; i < small; i++) {
      idxA[i] = static_cast<int>(idxDir[i]);
    }
    y.set_size(1, bevImage.size(1));
    small = bevImage.size(1) - 1;
    for (i = 0; i <= small; i++) {
      y[i] = static_cast<unsigned int>(i) + 1U;
    }
    idxDir.set_size(
        1, (static_cast<int>(approxLaneWidthPixels + 1.0) + y.size(1)) +
               static_cast<int>(approxLaneWidthPixels + 1.0));
    high_i = static_cast<int>(approxLaneWidthPixels + 1.0);
    for (i = 0; i < high_i; i++) {
      idxDir[i] = 1U;
    }
    small = y.size(1);
    for (i = 0; i < small; i++) {
      idxDir[i + static_cast<int>(approxLaneWidthPixels + 1.0)] = y[i];
    }
    for (i = 0; i < high_i; i++) {
      idxDir[(i + static_cast<int>(approxLaneWidthPixels + 1.0)) + y.size(1)] =
          sizeA_idx_1;
    }
    small = idxDir.size(1);
    for (i = 0; i < small; i++) {
      idxA[i + idxA.size(0)] = static_cast<int>(idxDir[i]);
    }
    sizeB_idx_1 = static_cast<double>(bevImage.size(1)) +
                  2.0 * (approxLaneWidthPixels + 1.0);
    i = static_cast<int>(sizeB_idx_1);
    Ipad.set_size(bevImage.size(0), i);
    for (nx = 0; nx < i; nx++) {
      k = Ipad.size(0);
      for (small = 0; small < k; small++) {
        Ipad[small + Ipad.size(0) * nx] =
            bevImage[(idxA[small] +
                      bevImage.size(0) * (idxA[nx + idxA.size(0)] - 1)) -
                     1];
      }
    }
  }
  q = 2.0 * (approxLaneWidthPixels + 1.0);
  sizeB_idx_1 = static_cast<double>(Ipad.size(1)) - q;
  if (sizeB_idx_1 < 1.0) {
    small = 0;
  } else {
    small = static_cast<int>(sizeB_idx_1);
  }
  Ileft.set_size(Ipad.size(0), small);
  for (i = 0; i < small; i++) {
    nx = Ipad.size(0);
    for (k = 0; k < nx; k++) {
      Ileft[k + Ileft.size(0) * i] = Ipad[k + Ipad.size(0) * i];
    }
  }
  if (q + 1.0 > Ipad.size(1)) {
    i = 0;
    k = 0;
  } else {
    i = static_cast<int>(q + 1.0) - 1;
    k = Ipad.size(1);
  }
  nx = Ipad.size(0);
  high_i = k - i;
  for (k = 0; k < high_i; k++) {
    for (mid_i = 0; mid_i < nx; mid_i++) {
      Ipad[mid_i + nx * k] = Ipad[mid_i + Ipad.size(0) * (i + k)];
    }
  }
  Ipad.set_size(Ipad.size(0), high_i);
  if ((Ileft.size(0) == Ipad.size(0)) && (Ileft.size(1) == high_i)) {
    x.set_size(Ileft.size(0), Ileft.size(1));
    small = Ileft.size(0) * Ileft.size(1);
    for (i = 0; i < small; i++) {
      x[i] = static_cast<short>(static_cast<int>(Ileft[i]) - Ipad[i]);
    }
  } else {
    b_binary_expand_op(x, Ileft, Ipad);
  }
  nx = x.size(0) * x.size(1);
  r.set_size(x.size(0), x.size(1));
  for (k = 0; k < nx; k++) {
    r[k] = static_cast<unsigned char>(std::abs(static_cast<float>(x[k])));
  }
  if (Ileft.size(0) == 1) {
    i = Ipad.size(0);
  } else {
    i = Ileft.size(0);
  }
  if (Ileft.size(1) == 1) {
    k = Ipad.size(1);
  } else {
    k = Ileft.size(1);
  }
  if (bevImage.size(0) == 1) {
    small = i;
  } else {
    small = bevImage.size(0);
  }
  if (bevImage.size(1) == 1) {
    high_i = k;
  } else {
    high_i = bevImage.size(1);
  }
  if ((Ileft.size(0) == Ipad.size(0)) && (Ileft.size(1) == Ipad.size(1)) &&
      (bevImage.size(0) == i) && (bevImage.size(1) == k) &&
      (small == r.size(0)) && (high_i == r.size(1))) {
    small = bevImage.size(0) * bevImage.size(1);
    Ileft.set_size(bevImage.size(0), bevImage.size(1));
    for (i = 0; i < small; i++) {
      Ileft[i] = static_cast<short>(
          static_cast<short>(static_cast<short>((bevImage[i] << 1) -
                                                static_cast<short>(Ileft[i])) -
                             Ipad[i]) -
          r[i]);
    }
  } else {
    b_binary_expand_op(Ileft, bevImage, Ipad, r);
  }
  high_i = Ileft.size(0) * Ileft.size(1);
  if (high_i <= 2) {
    if (high_i == 1) {
      small = static_cast<short>(Ileft[0]);
      nx = static_cast<short>(Ileft[0]);
    } else {
      nx = static_cast<short>(Ileft[high_i - 1]);
      if (static_cast<short>(Ileft[0]) > nx) {
        small = nx;
      } else {
        small = static_cast<short>(Ileft[0]);
      }
      nx = static_cast<short>(Ileft[high_i - 1]);
      if (static_cast<short>(Ileft[0]) >= nx) {
        nx = static_cast<short>(Ileft[0]);
      }
    }
  } else {
    small = static_cast<short>(Ileft[0]);
    for (k = 2; k <= high_i; k++) {
      i = static_cast<short>(Ileft[k - 1]);
      if (small > i) {
        small = i;
      }
    }
    nx = static_cast<short>(Ileft[0]);
    for (k = 2; k <= high_i; k++) {
      i = static_cast<short>(Ileft[k - 1]);
      if (nx < i) {
        nx = i;
      }
    }
  }
  nx -= small;
  if (std::abs(static_cast<float>(nx)) < 2.22044605E-15F) {
    unsigned int sizeA_idx_0;
    sizeA_idx_0 = static_cast<unsigned int>(Ileft.size(0));
    sizeA_idx_1 = static_cast<unsigned int>(Ileft.size(1));
    Ileft.set_size(static_cast<int>(sizeA_idx_0),
                   static_cast<int>(sizeA_idx_1));
    small = static_cast<int>(sizeA_idx_0) * static_cast<int>(sizeA_idx_1);
    for (i = 0; i < small; i++) {
      Ileft[i] = 0.0F;
    }
  } else {
    for (i = 0; i < high_i; i++) {
      Ileft[i] =
          (Ileft[i] - static_cast<float>(small)) / static_cast<float>(nx);
    }
  }
  std::memset(&ni[0], 0, 256U * sizeof(int));
  nx_tmp = Ileft.size(0) * Ileft.size(1);
  for (k = 0; k < nx_tmp; k++) {
    if ((Ileft[k] >= 0.0F) && (Ileft[k] <= 1.0F)) {
      float bGuess;
      bGuess = std::ceil(Ileft[k] / 0.00390625F);
      if ((bGuess >= 1.0F) && (bGuess < 257.0F) &&
          (Ileft[k] >= 0.00390625 * (bGuess - 1.0)) &&
          (Ileft[k] < 0.00390625 * bGuess)) {
        nx = static_cast<int>(bGuess) - 1;
        ni[nx]++;
      } else {
        nx = 1;
        small = 2;
        high_i = 257;
        while (high_i > small) {
          mid_i = (nx + high_i) >> 1;
          if (Ileft[k] >= 0.00390625 * (static_cast<double>(mid_i) - 1.0)) {
            nx = mid_i;
            small = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }
        ni[nx - 1]++;
      }
    }
  }
  for (i = 0; i < 256; i++) {
    counts[i] = ni[i];
  }
  for (k = 0; k < 255; k++) {
    counts[k + 1] += counts[k];
  }
  sizeB_idx_1 =
      ((1.0 - sensitivity) / 10.0 + 0.9) * static_cast<double>(nx_tmp);
  small = 0;
  k = 1;
  nx = 0;
  exitg1 = false;
  while ((!exitg1) && (nx < 256)) {
    if (counts[nx] >= sizeB_idx_1) {
      small = 1;
      ii_data = static_cast<short>(nx + 1);
      exitg1 = true;
    } else {
      nx++;
    }
  }
  if (small == 0) {
    k = 0;
  }
  for (i = 0; i < k; i++) {
    sizeB_idx_1 = static_cast<double>(ii_data) / 256.0;
  }
  nx = Ileft.size(0);
  i = Ileft.size(0);
  r1.set_size(Ileft.size(0), k * Ileft.size(1));
  high_i = Ileft.size(1);
  for (small = 0; small < high_i; small++) {
    mid_i = small * (nx * k);
    for (int jcol{0}; jcol < k; jcol++) {
      for (int itilerow{0}; itilerow < i; itilerow++) {
        r1[mid_i + itilerow] = sizeB_idx_1;
      }
    }
  }
  if (Ileft.size(0) == 1) {
    i = r1.size(0);
  } else {
    i = Ileft.size(0);
  }
  if (Ileft.size(1) == 1) {
    small = r1.size(1);
  } else {
    small = Ileft.size(1);
  }
  if ((Ileft.size(0) == r1.size(0)) && (Ileft.size(1) == r1.size(1)) &&
      (i == mask.size(0)) && (small == mask.size(1))) {
    BW.set_size(Ileft.size(0), Ileft.size(1));
    for (i = 0; i < nx_tmp; i++) {
      BW[i] = ((Ileft[i] > r1[i]) && (mask[i] != 0));
    }
  } else {
    binary_expand_op(BW, Ileft, r1, mask);
  }
  coder::bwareaopen(
      BW, std::fmax(
              std::round(0.0001 * static_cast<double>(BW.size(0) * BW.size(1))),
              30.0));
  q = std::ceil(std::sqrt((static_cast<double>(BW.size(0)) - 1.0) *
                              (static_cast<double>(BW.size(0)) - 1.0) +
                          (static_cast<double>(BW.size(1)) - 1.0) *
                              (static_cast<double>(BW.size(1)) - 1.0)));
  sizeB_idx_1 = 2.0 * q + 1.0;
  if (sizeB_idx_1 >= 0.0) {
    R.set_size(1, static_cast<int>(sizeB_idx_1));
    nx = static_cast<int>(sizeB_idx_1) - 1;
    R[static_cast<int>(sizeB_idx_1) - 1] = q;
    if (R.size(1) >= 2) {
      R[0] = -q;
      if (R.size(1) >= 3) {
        if (-q == -q) {
          sizeB_idx_1 = q / (static_cast<double>(R.size(1)) - 1.0);
          for (k = 2; k <= nx; k++) {
            R[k - 1] =
                static_cast<double>(((k << 1) - R.size(1)) - 1) * sizeB_idx_1;
          }
          if ((R.size(1) & 1) == 1) {
            R[R.size(1) >> 1] = 0.0;
          }
        } else {
          sizeB_idx_1 = (q - (-q)) / (static_cast<double>(R.size(1)) - 1.0);
          i = R.size(1);
          for (k = 0; k <= i - 3; k++) {
            R[k + 1] = -q + (static_cast<double>(k) + 1.0) * sizeB_idx_1;
          }
        }
      }
    }
  }
  coder::standardHoughTransformOptimized(BW, R, static_cast<double>(BW.size(0)),
                                         static_cast<double>(BW.size(1)), H);
  nx = H.size(0) * 180;
  if (H.size(0) * 180 <= 2) {
    if (H.size(0) * 180 == 1) {
      sizeB_idx_1 = H[0];
    } else {
      sizeB_idx_1 = H[H.size(0) * 180 - 1];
      if ((!(H[0] < sizeB_idx_1)) &&
          ((!std::isnan(H[0])) || std::isnan(sizeB_idx_1))) {
        sizeB_idx_1 = H[0];
      }
    }
  } else {
    if (!std::isnan(H[0])) {
      small = 1;
    } else {
      small = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= nx)) {
        if (!std::isnan(H[k - 1])) {
          small = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (small == 0) {
      sizeB_idx_1 = H[0];
    } else {
      sizeB_idx_1 = H[small - 1];
      i = small + 1;
      for (k = i; k <= nx; k++) {
        q = H[k - 1];
        if (sizeB_idx_1 < q) {
          sizeB_idx_1 = q;
        }
      }
    }
  }
  coder::houghpeaks(H, std::ceil(0.3 * sizeB_idx_1), tmp_data, tmp_size);
  coder::houghlines(BW, R, tmp_data, tmp_size, lines);
}

void detectLaneMarkerRidge_initialize()
{
  omp_init_nest_lock(&detectLaneMarkerRidge_nestLockGlobal);
  isInitialized_detectLaneMarkerRidge = true;
}

void detectLaneMarkerRidge_terminate()
{
  omp_destroy_nest_lock(&detectLaneMarkerRidge_nestLockGlobal);
  isInitialized_detectLaneMarkerRidge = false;
}

} // namespace detectLaneMarkerRidge
