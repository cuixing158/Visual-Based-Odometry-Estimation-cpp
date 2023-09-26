#ifndef ESTGEOTFORM2DFORPTSANDLINES_H
#define ESTGEOTFORM2DFORPTSANDLINES_H

#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

namespace estgeotform2dForPtsAndLines {
extern void estgeotform2dForPtsAndLines(
    const ::coder::array<double, 2U> &matchedPoints1,
    const ::coder::array<double, 2U> &matchedPoints2,
    const ::coder::array<double, 2U> &matchedLines1,
    const ::coder::array<double, 2U> &matchedLines2, double tform[6],
    ::coder::array<bool, 1U> &inlierIndex, double *status);

extern void estgeotform2dForPtsAndLines_initialize();

extern void estgeotform2dForPtsAndLines_terminate();

} // namespace estgeotform2dForPtsAndLines

#endif
