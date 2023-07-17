
#include "rt_nonfinite.h"
#include <math.h>

#if defined(__ICL) && __ICL == 1700
#pragma warning(disable : 264)
#endif

real_T rtNaN = (real_T)NAN;
real_T rtInf = (real_T)INFINITY;
real_T rtMinusInf = -(real_T)INFINITY;
real32_T rtNaNF = (real32_T)NAN;
real32_T rtInfF = (real32_T)INFINITY;
real32_T rtMinusInfF = -(real32_T)INFINITY;

#if defined(__ICL) && __ICL == 1700
#pragma warning(default : 264)
#endif

boolean_T rtIsInf(real_T value)
{
  return (isinf(value) != 0U);
}

boolean_T rtIsInfF(real32_T value)
{
  return (isinf((real_T)value) != 0U);
}

boolean_T rtIsNaN(real_T value)
{
  return (isnan(value) != 0U);
}

boolean_T rtIsNaNF(real32_T value)
{
  return (isnan((real_T)value) != 0U);
}
