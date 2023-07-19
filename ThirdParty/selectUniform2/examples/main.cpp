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
#include "selectUniform2.h"
#include "coder_array.h"

static void argInit_1x2_real_T(double result[2]);

static coder::array<double, 1U> argInit_Unboundedx1_real_T();

static coder::array<double, 2U> argInit_Unboundedx2_real_T();

static double argInit_real_T();

static void argInit_1x2_real_T(double result[2])
{

  for (int idx1{0}; idx1 < 2; idx1++) {

    result[idx1] = argInit_real_T();
  }
}

static coder::array<double, 1U> argInit_Unboundedx1_real_T()
{
  coder::array<double, 1U> result;

  result.set_size(2);

  for (int idx0{0}; idx0 < result.size(0); idx0++) {

    result[idx0] = argInit_real_T();
  }
  return result;
}

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

static double argInit_real_T()
{
  return 0.0;
}

int main(int, char **)
{

  main_selectUniform2();

  selectUniform2::selectUniform2_terminate();
  return 0;
}

void main_selectUniform2()
{
  coder::array<double, 2U> points;
  coder::array<double, 2U> pointsOut;
  coder::array<double, 1U> b_index;
  coder::array<double, 1U> responses;
  double dv[2];

  points = argInit_Unboundedx2_real_T();

  responses = argInit_Unboundedx1_real_T();

  argInit_1x2_real_T(dv);
  selectUniform2::selectUniform2(points, responses, argInit_real_T(), dv,
                                 pointsOut, b_index);
}
