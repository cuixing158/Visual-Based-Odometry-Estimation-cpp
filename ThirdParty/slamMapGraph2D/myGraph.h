///
/// @file           : myGraph.h
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

#ifndef MYGRAPH_H
#define MYGRAPH_H

/// @include file    : Include Files
#include "rtwtypes.h"
#include "slamMapGraph2D_types.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

/// Type Definitions
namespace SlamGraph2D {
class myGraph {
   public:
    myGraph();
    ~myGraph();
    void slamMapGraph2D(const struct0_T *poseParams, const double measurement[3],
                        double fromID, double toID, const double guesspose[3],
                        ::coder::array<double, 2U> &poses1);
    slamMapGraph2DStackData *getStackData();

   private:
    slamMapGraph2DPersistentData pd_;
    slamMapGraph2DStackData SD_;
};

}  // namespace SlamGraph2D

#endif
///
/// File trailer for myGraph.h
///
/// [EOF]
///
