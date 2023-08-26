#ifndef DETECTLANEMARKERRIDGE_H
#define DETECTLANEMARKERRIDGE_H

#include "detectLaneMarkerRidge_types.h"
#include "rtwtypes.h"
#include "coder_array.h"
#include "omp.h"
#include <cstddef>
#include <cstdlib>

namespace detectLaneMarkerRidge {
extern omp_nest_lock_t detectLaneMarkerRidge_nestLockGlobal;
}

namespace detectLaneMarkerRidge {
extern void
detectLaneMarkerRidge(const ::coder::array<unsigned char, 2U> &bevImage,
                      double approxLaneWidthPixels,
                      const ::coder::array<unsigned char, 2U> &mask,
                      double sensitivity, ::coder::array<struct0_T, 2U> &lines);

extern void detectLaneMarkerRidge_initialize();

extern void detectLaneMarkerRidge_terminate();

} // namespace detectLaneMarkerRidge

#endif
