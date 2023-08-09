#ifndef POSEGRAPHOPTIMIZE_H
#define POSEGRAPHOPTIMIZE_H

#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

namespace poseGraphOptimize {
extern void poseGraphOptimize(const ::coder::array<double, 2U> &absposes,
                              const ::coder::array<double, 2U> &loopNodePairs,
                              const ::coder::array<double, 2U> &loopPoses,
                              ::coder::array<double, 2U> &updatedPoses);

extern void poseGraphOptimize_initialize();

extern void poseGraphOptimize_terminate();

} // namespace poseGraphOptimize

#endif
