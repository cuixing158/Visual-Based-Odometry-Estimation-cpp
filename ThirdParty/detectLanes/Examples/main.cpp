/*************************************************************************/
/* This automatically generated example C++ main file shows how to call  */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

#include "main.h"
#include "detectLaneMarkerRidge.h"
#include "detectLaneMarkerRidge_types.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

static coder::array<unsigned char, 2U> argInit_UnboundedxUnbounded_uint8_T();

static double argInit_real_T();

static unsigned char argInit_uint8_T();

static coder::array<unsigned char, 2U> argInit_UnboundedxUnbounded_uint8_T()
{
  coder::array<unsigned char, 2U> result;

  result.set_size(2, 2);

  for (int idx0{0}; idx0 < result.size(0); idx0++) {
    for (int idx1{0}; idx1 < result.size(1); idx1++) {

      result[idx0 + result.size(0) * idx1] = argInit_uint8_T();
    }
  }
  return result;
}

static double argInit_real_T()
{
  return 0.0;
}

static unsigned char argInit_uint8_T()
{
  return 0U;
}

int main(int, char **)
{

  main_detectLaneMarkerRidge();

  detectLaneMarkerRidge::detectLaneMarkerRidge_terminate();
  return 0;
}

void main_detectLaneMarkerRidge()
{
  coder::array<detectLaneMarkerRidge::struct0_T, 2U> lines;
  coder::array<unsigned char, 2U> bevImage_tmp;
  double approxLaneWidthPixels_tmp;

  bevImage_tmp = argInit_UnboundedxUnbounded_uint8_T();
  approxLaneWidthPixels_tmp = argInit_real_T();

  detectLaneMarkerRidge::detectLaneMarkerRidge(
      bevImage_tmp, approxLaneWidthPixels_tmp, bevImage_tmp,
      approxLaneWidthPixels_tmp, lines);
}
