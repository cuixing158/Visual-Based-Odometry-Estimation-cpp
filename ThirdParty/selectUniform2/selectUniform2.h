#ifndef SELECTUNIFORM2_H
#define SELECTUNIFORM2_H

#include "rtwtypes.h"
#include "selectUniform2_types.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void selectUniform2(emxArray_real_T *points, emxArray_real_T *responses,
                           double N, const double imageSize[2],
                           emxArray_real_T *pointsOut,
                           emxArray_real_T *b_index);

extern void selectUniform2_initialize(void);

extern void selectUniform2_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
