#ifndef HDMAPPING_H
#define HDMAPPING_H

#include <vector>
#include <random>
#include "opencv2/opencv.hpp"
#include "constructWorldMap_types.h"

// #include "spdlog/spdlog.h"

#define DEBUG_SHOW_IMAGE 0
namespace buildMapping {
class HDMapping {
   public:
    HDMapping();
    ~HDMapping();
    void constructWorldMap(const cv::Mat& srcImage);

   public:
    HDmap HDmapOutput;
    std::vector<cv::Vec3f> vehiclePoses;
    double cumDist;
    double pixelExtentInWorldXY;
    bool isBuildMap;
    double buildMapStopFrame;
    bool isBuildMapOver;
    bool isLocSuccess;
    double locVehiclePose[3];

   private:
    // 预定义变量
    cv::Mat orbDetectMask, BW;  // 分别对应MATLAB中的vehiclePolygon所围成的mask，BW变量，只用于更新的ROI
    std::vector<cv::KeyPoint> preKeypts, currKeypts;
    cv::Mat preDescriptors, currDescriptors;
    cv::Mat prevImg;
    cv::Mat preRelTform, relTform, previousImgPose;  // 2*3转换矩阵

    //
    cv::Ptr<cv::Feature2D> orbDetector;

    void selectUniform(std::vector<cv::KeyPoint>& keypts, cv::Mat& Descriptions, size_t numPoints, std::vector<cv::KeyPoint>& outKeypts, cv::Mat& outDescriptions);
};

}  // namespace buildMapping

#endif
