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

void matConvertToYUV() {}

void YUVConvertToMat(const char* yuvImageFile, unsigned char* image_data) {
    FILE* fp = fopen(yuvImageFile, "rb");
    unsigned char image_data[480 * 640];
    fread(image_data, sizeof(unsigned char), 480 * 640, fp);
    fclose(fp);
}

int main(int, char**) {
    std::string imagePath = "/opt_disk2/rd22946/vscode_work/python_works/mytest/test";  // /opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/116_new_undistort/116";  //"/opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/map_R_new_undistort/map_R";
    std::vector<std::string> imagePaths;
    size_t numImgs = getFullNames(filesystem::path(imagePath), imagePaths, ".yuv");
    // std::sort(imagePaths.begin(), imagePaths.end(),
    //           [](std::string p1, std::string p2) { return atoi(filesystem::path(p1).filenameNoExt().substr(6).c_str()) < atoi(filesystem::path(p2).filenameNoExt().substr(6).c_str()); });
    std::sort(imagePaths.begin(), imagePaths.end(),
              [](std::string p1, std::string p2) { return atoi(filesystem::path(p1).filenameNoExt().c_str()) < atoi(filesystem::path(p2).filenameNoExt().c_str()); });
    std::ofstream fid("imageFilesList.txt");

    buildMapping::HDMapping obj;
    double cumTime = 0.0;
    bool flag = false;  // 建图终止标志，暴露给用户控制
    buildMapping::HDMapping::buildMapStatus statusB = buildMapping::HDMapping::buildMapStatus::BUILD_MAP_PROCESSING;
    buildMapping::HDMapping::localizeMapStatus statusL = buildMapping::HDMapping::localizeMapStatus::LOCALIZE_MAP_PROCESSING;
    for (size_t i = 351; i < 1160; i++) {  // from 351, to 1160 for 116_new_undistort/116
        fid << imagePaths[i] << std::endl;
        // cv::Mat srcImage = cv::imread(imagePaths[i]);
        cv::Mat srcImage = cv::Mat::zeros(480, 640, CV_8UC1);
        YUVConvertToMat(char(imagePaths[i]), srcImage.data);
        if (i == 1159) {
            flag = true;
        }

        // main loop
        double t1 = cv::getTickCount();
        statusB = obj.constructWorldMap(srcImage, flag);
        // statusL = obj.localizeWorldMap(srcImage);
        double elapseTime = (double)(cv::getTickCount() - t1) / cv::getTickFrequency();
        cumTime += elapseTime;
        printf("%zd,Elapsed second Time:%.5f,avg time:%.6f,status:%d\n", i, elapseTime, cumTime / (i + 1), statusB);

        // result
        if (statusB == buildMapping::HDMapping::buildMapStatus::BUILD_MAP_SUCCESSFUL) {
            break;
        }
        if (statusL == buildMapping::HDMapping::localizeMapStatus::LOCALIZE_MAP_SUCCESSFUL) {
            for (int j = 0; j < 3; j++) {
                printf("m_locVehiclePose[%d]=%.3f\n", j, obj.m_locVehiclePose[j]);
            }
            break;
        }
    }
    fid.close();
    return 0;
}
