#ifndef CONSTRUCTWORLDMAP_TYPES_H
#define CONSTRUCTWORLDMAP_TYPES_H

#include "opencv2/opencv.hpp"

#define MAX_THREADS omp_get_max_threads()

namespace buildMapping {
struct constructWorldMapPersistentData;

}

namespace buildMapping {
struct struct0_T {
    unsigned char undistortImage[307200];
    double currFrontBasePose[3];
    bool isuseGT;
};

struct uint80m_T {
    unsigned long chunks[2];
};

struct uint160m_T {
    unsigned long chunks[4];
};

struct imref2d {
    double XWorldLimits[2];
    double YWorldLimits[2];
    double ImageSize[2];
};

struct HDmap {
    cv::Mat bigImg;
    imref2d ref;
};

struct cerealPoses {
    double x;
    double y;
    double theta;
};
struct struct1_T {
    HDmap HDmapData;
    std::vector<cv::Vec3f> vehiclePoses;
    double cumDist;
    double pixelExtentInWorldXY;
    bool isBuildMap;
    double buildMapStopFrame;
    bool isBuildMapOver;
    bool isLocSuccess;
    double locVehiclePose[3];
};

struct imresize {
    unsigned char out[76800];
};

struct helperDetectAndExtractFeatures {
    unsigned char Iu8[307200];
};

struct fuseOptimizeHDMap {
    unsigned char currImg[307200];
};

struct constructWorldMapStackData {
    imresize f0;
    helperDetectAndExtractFeatures f1;
    fuseOptimizeHDMap f2;
    constructWorldMapPersistentData *pd;
};

}  // namespace buildMapping

#endif
