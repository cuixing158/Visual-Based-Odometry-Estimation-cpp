#ifndef SELECTUNIFORM2_H
#define SELECTUNIFORM2_H

#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

namespace selectUniform2 {
extern void selectUniform2(::coder::array<double, 2U> &points,
                           ::coder::array<double, 1U> &responses, double N,
                           const double imageSize[2],
                           ::coder::array<double, 2U> &pointsOut,
                           ::coder::array<double, 1U> &b_index);

extern void selectUniform2_initialize();

extern void selectUniform2_terminate();

} // namespace selectUniform2

#endif
