///
/// @file           : mtimes.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

/// @include file    : Include Files
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cstring>

/// Function Definitions
///
/// @fn             : mtimes
/// @brief          :
/// @param          : const double A_data[]
///                   const int A_size[2]
///                   const ::coder::array<double, 2U> &B
///                   double C_data[]
///                   int C_size[2]
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace internal {
namespace blas {
void mtimes(const double A_data[], const int A_size[2],
            const ::coder::array<double, 2U> &B, double C_data[], int C_size[2])
{
  int inner;
  int mc;
  int nc;
  mc = A_size[1];
  inner = A_size[0];
  nc = B.size(1);
  C_size[0] = A_size[1];
  C_size[1] = B.size(1);
  for (int j{0}; j < nc; j++) {
    int boffset;
    int coffset;
    coffset = j * mc;
    boffset = j * B.size(0);
    std::memset(&C_data[coffset], 0,
                static_cast<unsigned int>((mc + coffset) - coffset) *
                    sizeof(double));
    for (int k{0}; k < inner; k++) {
      double bkj;
      bkj = B[boffset + k];
      for (int i{0}; i < mc; i++) {
        int b_i;
        b_i = coffset + i;
        C_data[b_i] += A_data[i * A_size[0] + k] * bkj;
      }
    }
  }
}

///
/// @fn             : mtimes
/// @brief          :
/// @param          : const double A_data[]
///                   const int A_size[2]
///                   const double B_data[]
///                   const int B_size[2]
///                   double C_data[]
///                   int C_size[2]
/// @return         : void
///
void mtimes(const double A_data[], const int A_size[2], const double B_data[],
            const int B_size[2], double C_data[], int C_size[2])
{
  int inner;
  int mc;
  int nc;
  mc = A_size[0];
  inner = A_size[1];
  nc = B_size[1];
  C_size[0] = A_size[0];
  C_size[1] = B_size[1];
  for (int j{0}; j < nc; j++) {
    int boffset;
    int coffset;
    coffset = j * mc;
    boffset = j * B_size[0];
    std::memset(&C_data[coffset], 0,
                static_cast<unsigned int>((mc + coffset) - coffset) *
                    sizeof(double));
    for (int k{0}; k < inner; k++) {
      double bkj;
      int aoffset;
      aoffset = k * A_size[0];
      bkj = B_data[boffset + k];
      for (int i{0}; i < mc; i++) {
        int b_i;
        b_i = coffset + i;
        C_data[b_i] += A_data[aoffset + i] * bkj;
      }
    }
  }
}

} // namespace blas
} // namespace internal
} // namespace coder
} // namespace SlamGraph2D

///
/// File trailer for mtimes.cpp
///
/// [EOF]
///
