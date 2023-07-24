///
/// @file           : main.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

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

/// @include file    : Include Files
#include "main.h"
#include "myGraph.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_types.h"
#include "coder_array.h"

/// Function Declarations
static void argInit_1x3_real_T(double result[3]);

static double argInit_real_T();

static SlamGraph2D::struct0_T argInit_struct0_T();

/// Function Definitions
///
/// @fn             : argInit_1x3_real_T
/// @brief          :
/// @param          : double result[3]
/// @return         : void
///
static void argInit_1x3_real_T(double result[3])
{
  // Loop over the array to initialize each element.
  for (int idx1{0}; idx1 < 3; idx1++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx1] = argInit_real_T();
  }
}

///
/// @fn             : argInit_real_T
/// @brief          :
/// @param          : void
/// @return         : double
///
static double argInit_real_T()
{
  return 0.0;
}

///
/// @fn             : argInit_struct0_T
/// @brief          :
/// @param          : void
/// @return         : SlamGraph2D::struct0_T
///
static SlamGraph2D::struct0_T argInit_struct0_T()
{
  SlamGraph2D::struct0_T result;
  double result_tmp;
  // Set the value of each structure field.
  // Change this value to the value that the application requires.
  result_tmp = argInit_real_T();
  result.MaxNumNodes = result_tmp;
  result.MaxNumEdges = result_tmp;
  return result;
}

///
/// @fn             : main
/// @brief          :
/// @param          : int argc
///                   char **argv
/// @return         : int
///
int main(int, char **)
{
  SlamGraph2D::myGraph *classInstance;
  classInstance = new SlamGraph2D::myGraph;
  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_slamMapGraph2D(classInstance);
  delete classInstance;
  return 0;
}

///
/// @fn             : main_slamMapGraph2D
/// @brief          :
/// @param          : SlamGraph2D::myGraph *instancePtr
/// @return         : void
///
void main_slamMapGraph2D(SlamGraph2D::myGraph *instancePtr)
{
  coder::array<double, 2U> poses1;
  SlamGraph2D::struct0_T r;
  double measurement_tmp[3];
  double fromID_tmp;
  // Initialize function 'slamMapGraph2D' input arguments.
  // Initialize function input argument 'poseParams'.
  // Initialize function input argument 'measurement'.
  argInit_1x3_real_T(measurement_tmp);
  fromID_tmp = argInit_real_T();
  // Initialize function input argument 'guesspose'.
  // Call the entry-point 'slamMapGraph2D'.
  r = argInit_struct0_T();
  instancePtr->slamMapGraph2D(&r, measurement_tmp, fromID_tmp, fromID_tmp,
                              measurement_tmp, poses1);
}

///
/// File trailer for main.cpp
///
/// [EOF]
///
