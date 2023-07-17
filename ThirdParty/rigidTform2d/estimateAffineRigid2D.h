#ifndef ESTIMATEAFFINERIGID2D_H
#define ESTIMATEAFFINERIGID2D_H

#include "estimateAffineRigid2D_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void estimateAffineRigid2D(const emxArray_real_T *pts1,
                                  const emxArray_real_T *pts2,
                                  double tform2x3[6],
                                  emxArray_boolean_T *inlierIndex, int *status);

extern void estimateAffineRigid2D_initialize(void);

extern void estimateAffineRigid2D_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
