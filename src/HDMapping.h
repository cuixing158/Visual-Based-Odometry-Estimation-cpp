#ifndef HDMAPPING_H
#define HDMAPPING_H

#include <vector>
#include <random>
#include "opencv2/opencv.hpp"
#include "constructWorldMap_types.h"

// DBoW3
#include "DBoW3.h"
#include "DescManip.h"

// #include "spdlog/spdlog.h"

namespace buildMapping {
class HDMapping {
   public:
    HDMapping();
    ~HDMapping();
    void constructWorldMap(const cv::Mat& srcImage);

   public:
    enum matchFeatureMethod { ORB_FEATURES,
                              LK_TRACK_FEATURES,
                              HYBRID_FEATURES };
    HDmap HDmapOutput;
    std::vector<cv::Vec3d> vehiclePoses;  // [x,y,theta]位姿,theta为弧度
    double cumDist;                       // 行驶累计距离，单位：米
    double pixelExtentInWorldXY;          // 每个像素实际物理长度距离，单位：米/像素
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
    matchFeatureMethod method;
    void selectUniformPoints(std::vector<cv::KeyPoint>& keyPoints, int numRetPoints,
                             cv::Size size, std::vector<cv::KeyPoint>& outputPts, std::vector<int>& indexs);

    void estiTform(std::vector<cv::Point2f>& prePoints, std::vector<cv::Point2f>& currPoints, cv::Mat& tform2x3, cv::Mat& inliers, int& status);

    // loop
    std::vector<cv::Mat> m_features;
    DBoW3::Database db;
    void loopDatabase_add_feature(cv::Mat& imageFeature);

    void reset();
    void saveMapData();
    void loadMapData();
};

}  // namespace buildMapping

#endif
