#ifndef DETECTLANEMARKERRIDGE_TYPES_H
#define DETECTLANEMARKERRIDGE_TYPES_H

#include "rtwtypes.h"
#define MAX_THREADS omp_get_max_threads()

namespace detectLaneMarkerRidge {
struct struct0_T {
  double point1[2];
  double point2[2];
  double theta;
  double rho;
};

} // namespace detectLaneMarkerRidge

#endif
