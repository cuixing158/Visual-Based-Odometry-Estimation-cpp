#ifndef ESTIMATEAFFINERIGID2D_H
#define ESTIMATEAFFINERIGID2D_H

#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

namespace estimateAffineRigid2D {
extern void estimateAffineRigid2D(const ::coder::array<double, 2U> &pts1,
                                  const ::coder::array<double, 2U> &pts2,
                                  double tform2x3[6],
                                  ::coder::array<boolean_T, 2U> &inlierIndex,
                                  int *status);

extern void estimateAffineRigid2D_initialize();

extern void estimateAffineRigid2D_terminate();

} // namespace estimateAffineRigid2D

#endif
