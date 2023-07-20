#ifndef HDMAPPING_H
#define HDMAPPING_H

#include <vector>
#include <random>
#include "opencv2/opencv.hpp"
#include "constructWorldMap_types.h"

// #include "spdlog/spdlog.h"

namespace buildMapping {
class HDMapping {
   public:
    HDMapping();
    ~HDMapping();
    void constructWorldMap(const cv::Mat& srcImage);

   public:
    HDmap HDmapOutput;
    std::vector<cv::Vec3d> vehiclePoses;  // [x,y,theta]位姿,theta为弧度
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
    cv::Vec3d initViclePtPose;                       // [x,y,theta]位姿

    //
    cv::Ptr<cv::Feature2D> orbDetector;

    void selectUniform(std::vector<cv::KeyPoint>& keypts, cv::Mat& Descriptions, size_t numPoints, std::vector<cv::KeyPoint>& outKeypts, cv::Mat& outDescriptions);

    void selectUniformPoints(std::vector<cv::KeyPoint>& keyPoints, int numRetPoints,
                             cv::Size size, std::vector<cv::KeyPoint>& outputPts, std::vector<int>& indexs);

    void estiTform(std::vector<cv::Point>& prePoints, std::vector<cv::Point>& currPoints, cv::Mat& tform2x3, cv::Mat& inliers, int& status);
};

}  // namespace buildMapping

#endif
