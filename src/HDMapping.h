#ifndef HDMAPPING_H
#define HDMAPPING_H

#include <vector>
#include <random>
#include "opencv2/opencv.hpp"
#include "constructWorldMap_types.h"

// DBoW3
#include "DBoW3.h"
#include "DescManip.h"

//loop
#include "slamPoseGraph.h"
// #include "spdlog/spdlog.h"

namespace buildMapping {

typedef struct imageKptsAndFeatures {
    std::vector<cv::KeyPoint> keyPoints;
    cv::Mat features;
} imageKptsAndFeatures;
class HDMapping {
   public:
    HDMapping();
    ~HDMapping();
    void constructWorldMap(const cv::Mat& srcImage, bool isStopConstructWorldMap = false);

    void localizeWorldMap();

   public:
    enum matchFeatureMethod { ORB_FEATURES,
                              LK_TRACK_FEATURES,
                              HYBRID_FEATURES };

    HDmap m_HDmapOutput;
    std::vector<cv::Vec3d> m_vehiclePoses;  // [x,y,theta]位姿,theta为弧度
    double m_cumDist;                       // 行驶累计距离，单位：米
    double m_pixelExtentInWorldXY;          // 每个像素实际物理长度距离，单位：米/像素
    bool m_isBuildMap;
    double m_buildMapStopFrame;
    bool m_isBuildMapOver;
    bool m_isLocSuccess;
    double m_locVehiclePose[3];

   private:
    // 预定义变量
    cv::Mat m_orbDetectMask, m_BW;  // 分别对应MATLAB中的vehiclePolygon所围成的mask，BW变量，只用于更新的ROI
    std::vector<cv::KeyPoint> m_preKeypts, m_currKeypts;
    cv::Mat m_preDescriptors, m_currDescriptors;
    cv::Mat m_prevImg;
    cv::Mat m_preRelTform, m_relTform, m_previousImgPose;  // 2*3转换矩阵
    cv::Vec3d m_initViclePtPose, m_preViclePtPose;         // [x,y,theta]位姿

    //
    cv::Ptr<cv::Feature2D> m_orbDetector;
    matchFeatureMethod m_method;
    void selectUniformPoints(std::vector<cv::KeyPoint>& keyPoints, int numRetPoints,
                             cv::Size size, std::vector<cv::KeyPoint>& outputPts, std::vector<int>& indexs);

    void estiTform(std::vector<cv::Point2f>& prePoints, std::vector<cv::Point2f>& currPoints, cv::Mat& tform2x3, cv::Mat& inliers, int& status);

    // loop
    std::vector<imageKptsAndFeatures> m_points_features;
    DBoW3::Database m_db;
    SlamGraph2D::myGraph m_instancePtr;
    SlamGraph2D::slamPoseGraph m_pg;  // 索引从1开始
    void detectLoopAndAddGraph();
    void loopDatabaseAddFeaturesAndSave(std::string saveDataBaseYmlGz = "./database.yml.gz");
    DBoW3::QueryResults buildMapping::HDMapping::retrieveImages(cv::Mat queryImage);

    void reset();
    void saveMapData();
    void loadMapData();
};

}  // namespace buildMapping

#endif
