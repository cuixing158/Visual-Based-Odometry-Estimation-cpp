#include "external.h"
#include <iostream>
#include <algorithm>
#include <fstream>

#include "opencv2/opencv.hpp"

#include "c_cpp_utils/path.h"
#include "src/HDMapping.h"
#ifdef __cplusplus
extern "C" {
#endif
static buildMapping::HDMapping obj;

void maping_world(int width,int height,unsigned char*image_data,char flag )
{

    cv::Mat srcImage(width,height,CV_8UC1,image_data);
    cv::imshow("show",srcImage);
    cv::waitKey(1);
    buildMapping::HDMapping::buildMapStatus statusB = buildMapping::HDMapping::buildMapStatus::BUILD_MAP_PROCESSING;
    statusB = obj.constructWorldMap(srcImage, flag);
}

void init_map_src()
{
    std::string imagePath = "/home/linzhiqiang/D/buildMapping_CPP/map_R_new_undistort/map_R";  // /opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/116_new_undistort/116";  //"/opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/map_R_new_undistort/map_R";
    std::vector<std::string> imagePaths;
    size_t numImgs = getFullNames(filesystem::path(imagePath), imagePaths, ".jpg");
    std::sort(imagePaths.begin(), imagePaths.end(),
               [](std::string p1, std::string p2) { return atoi(filesystem::path(p1).filenameNoExt().c_str()) < atoi(filesystem::path(p2).filenameNoExt().c_str()); });
    std::ofstream fid("imageFilesList.txt");
    for (size_t i = 350; i < 1160; i++) {  // from 351, to 1160 for 116_new_undistort/116
        fid << imagePaths[i] << std::endl;
    }
    fid.close();
}

#ifdef __cplusplus
};
#endif
