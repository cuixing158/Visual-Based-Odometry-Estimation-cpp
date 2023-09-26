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
#include "estgeotform2dForPtsAndLines.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

static coder::array<double, 2U> argInit_Unboundedx2_real_T();

static coder::array<double, 2U> argInit_Unboundedx4_real_T();

static double argInit_real_T();

static coder::array<double, 2U> argInit_Unboundedx2_real_T()
{
  coder::array<double, 2U> result;

  result.set_size(2, 2);

  for (int idx0{0}; idx0 < result.size(0); idx0++) {
    for (int idx1{0}; idx1 < 2; idx1++) {

      result[idx0 + result.size(0) * idx1] = argInit_real_T();
    }
  }
  return result;
}

static coder::array<double, 2U> argInit_Unboundedx4_real_T()
{
  coder::array<double, 2U> result;

  result.set_size(2, 4);

  for (int idx0{0}; idx0 < result.size(0); idx0++) {
    for (int idx1{0}; idx1 < 4; idx1++) {

      result[idx0 + result.size(0) * idx1] = argInit_real_T();
    }
  }
  return result;
}

static double argInit_real_T()
{
  return 0.0;
}

int main(int, char **)
{

  main_estgeotform2dForPtsAndLines();

  estgeotform2dForPtsAndLines::estgeotform2dForPtsAndLines_terminate();
  return 0;
}

void main_estgeotform2dForPtsAndLines()
{
  coder::array<double, 2U> matchedLines1_tmp;
  coder::array<double, 2U> matchedPoints1_tmp;
  coder::array<bool, 1U> inlierIndex;
  double tform[6];
  double status;

  matchedPoints1_tmp = argInit_Unboundedx2_real_T();

  matchedLines1_tmp = argInit_Unboundedx4_real_T();

  estgeotform2dForPtsAndLines::estgeotform2dForPtsAndLines(
      matchedPoints1_tmp, matchedPoints1_tmp, matchedLines1_tmp,
      matchedLines1_tmp, tform, inlierIndex, &status);
}
