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
#include "estimateAffineRigid2D.h"

int main(int, char **) {
    std::string imagePath = "opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/116_new_undistort/116";  // opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/116_new_undistort/116";  //"/opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/map_R_new_undistort/map_R";
    std::vector<std::string> imagePaths;
    size_t numImgs = getFullNames(filesystem::path(imagePath), imagePaths, ".jpg");
    // std::sort(imagePaths.begin(), imagePaths.end(),
    //           [](std::string p1, std::string p2) { return atoi(filesystem::path(p1).filenameNoExt().substr(6).c_str()) < atoi(filesystem::path(p2).filenameNoExt().substr(6).c_str()); });
    std::sort(imagePaths.begin(), imagePaths.end(),
              [](std::string p1, std::string p2) { return atoi(filesystem::path(p1).filenameNoExt().c_str()) < atoi(filesystem::path(p2).filenameNoExt().c_str()); });

    int status = 0;
    coder::array<double, 2U> pts1_tmp, pts2_tmp;
    coder::array<boolean_T, 2U> inlierIndex;
    double tform2x3_[6];

    pts1_tmp.set_size(5, 2);
    pts2_tmp.set_size(5, 2);
    for (size_t i = 0; i < 5; i++) {
        pts1_tmp[i] = 1;
        pts1_tmp[i + 5] = 2;

        pts2_tmp[i] = 1;
        pts2_tmp[i + 5] = 2;
    }
    estimateAffineRigid2D::estimateAffineRigid2D(pts1_tmp, pts2_tmp, tform2x3_, inlierIndex, &status);
    cv::Mat inliers = cv::Mat(5, 1, CV_8U, inlierIndex.data()).clone();
    cv::Mat tform2x3 = (cv::Mat_<double>(2, 3) << tform2x3_[0], tform2x3_[2], tform2x3_[4], tform2x3_[1], tform2x3_[3], tform2x3_[5]);
    std::cout << "tform2x3:" << tform2x3 << std::endl;

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