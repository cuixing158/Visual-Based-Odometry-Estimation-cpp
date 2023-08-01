///
/// @file           : main.cpp
/// @target         : Texas Instruments->C6000
/// @details        : for path build map algorithms(loop+pose)
/// @author         : cuixingxing
/// @email          : xingxing.cui@long-horn.com
/// @date           : 02-Mar-2023 10:06:28
/// @version        : V0.1.2
/// @copyright      : Copyright (C) 2023 Long-Horn Inc.All rights reserved.
///

#include <algorithm>
#include <iostream>
#include <fstream>
#include "opencv2/opencv.hpp"

#include "c_cpp_utils/path.h"
#include "src/HDMapping.h"

int main(int, char **) {
    std::string imagePath = "opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/116_new_undistort/116";  // opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/116_new_undistort/116";  //"/opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/map_R_new_undistort/map_R";
    std::vector<std::string> imagePaths;
    size_t numImgs = getFullNames(filesystem::path(imagePath), imagePaths, ".jpg");
    // std::sort(imagePaths.begin(), imagePaths.end(),
    //           [](std::string p1, std::string p2) { return atoi(filesystem::path(p1).filenameNoExt().substr(6).c_str()) < atoi(filesystem::path(p2).filenameNoExt().substr(6).c_str()); });
    std::sort(imagePaths.begin(), imagePaths.end(),
              [](std::string p1, std::string p2) { return atoi(filesystem::path(p1).filenameNoExt().c_str()) < atoi(filesystem::path(p2).filenameNoExt().c_str()); });

    buildMapping::HDMapping obj;
    double cumTime = 0.0;
    for (size_t i = 1; i < numImgs; i++) {  // from 351, to 1160 for 116_new_undistort/116
        std::cout << imagePaths[i] << std::endl;
        cv::Mat srcImage = cv::imread(imagePaths[i]);
        double t1 = cv::getTickCount();
        obj.constructWorldMap(srcImage);
        double elapseTime = (double)(cv::getTickCount() - t1) / cv::getTickFrequency();
        cumTime += elapseTime;
        printf("Elapsed second Time:%.5f,avg time:%.6f\n", elapseTime, cumTime / (i + 1));
    }

    return 0;
}