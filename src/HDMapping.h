#pragma once

#include <vector>
#include <random>
#include "opencv2/opencv.hpp"
#include <opencv2/line_descriptor.hpp>
#include "constructWorldMap_types.h"

#include "estimateAffineRigid2D.h"
#include "selectUniform2.h"
#include "coder_array.h"

#include "c_cpp_utils/path.h"

// DBoW3
#include "DBoW3.h"
#include "DescManip.h"

// origin loop
#include "slamPoseGraph.h"
#include "myGraph.h"
#include "poseGraph.h"

// current loop
#include "poseGraphOptimize.h"

// cereal header files
#include "cereal/archives/portable_binary.hpp"
#include "cereal/archives/json.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/vector.hpp"
#include <cereal/types/string.hpp>
#include <fstream>

namespace buildMapping {

// reference: DBOW3/utils/demo_cereal_bench.cpp
template <typename _Tp>
std::vector<_Tp> convertMat2Vector(const cv::Mat& mat) {
    return (std::vector<_Tp>)(mat.reshape(1, 1));  //通道数不变，按行转为一行
}

template <typename _Tp>
cv::Mat convertVector2Mat(std::vector<_Tp> v, int channels, int rows) {
    cv::Mat mat = cv::Mat(v);                            //将vector变成单列的mat
    cv::Mat dest = mat.reshape(channels, rows).clone();  //PS：必须clone()一份，否则返回出错
    return dest;
}

template <class Archive>
void serialize(Archive& ar, cv::Mat& feats) {
    std::vector<uchar> vecFeats = convertMat2Vector<uchar>(feats);
    ar(feats.rows, feats.cols, feats.channels(), vecFeats);
}

template <class Archive>
void serialize(Archive& ar, cerealPoses& pose) {
    ar(pose.x, pose.y, pose.theta);
}

typedef struct keypt {
    float x;
    float y;
    float size;

    template <class Archive>
    void serialize(Archive& ar) {
        ar(x, y, size);
    }
} keypt;

template <class T>
struct imgInfo {
    int rows;
    int cols;
    int channels;
    std::vector<T> datas;

    template <class Archive>
    void serialize(Archive& ar) {
        ar(rows, cols, channels, datas);
    }
};
typedef struct imageKptsAndFeatures {
    std::vector<cv::KeyPoint> keyPoints;
    cv::Mat features;

    template <class Archive>
    void save(Archive& ar) const {
        std::vector<keypt> kpts;
        for (size_t i = 0; i < keyPoints.size(); i++) {
            kpts.push_back(keypt{keyPoints[i].pt.x, keyPoints[i].pt.y, keyPoints[i].size});
        }
        std::vector<uint8_t> data_temp = convertMat2Vector<uint8_t>(features);
        struct imgInfo<uint8_t> feats {
            features.rows, features.cols, features.channels(), data_temp
        };
        ar(kpts, feats);
    }

    template <class Archive>
    void load(Archive& ar) {
        std::vector<keypt> kpts;
        struct imgInfo<uint8_t> feats;
        ar(kpts, feats);
        for (size_t i = 0; i < kpts.size(); i++) {
            keyPoints.push_back(cv::KeyPoint(kpts[i].x, kpts[i].y, kpts[i].size));
        }
        features = convertVector2Mat(feats.datas, feats.channels, feats.rows);
    }
} imageKptsAndFeatures;

template <typename T>
std::vector<int> findItems(std::vector<T> const& v, double greatThanTarget) {
    std::vector<int> indexs;
    auto it = v.begin();
    while ((it = std::find_if(it, v.end(), [&](T const& e) { return e >= greatThanTarget; })) != v.end()) {
        indexs.push_back(std::distance(v.begin(), it));
        it++;
    }
    return indexs;
}

template <typename T>
std::vector<std::vector<T>> nchoosek(std::vector<T> V, int K) {
    int N = V.size();
    std::string bitmask(K, 1);
    bitmask.resize(N, 0);

    std::vector<std::vector<T>> arr;
    do {
        std::vector<T> t_arr(K);
        int j = 0;
        for (int i = 0; i < N; ++i) {
            int ele = V[i];
            if (bitmask[i]) {
                t_arr[j++] = ele;
            }
        }
        arr.push_back(t_arr);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    return arr;
}

template <typename T>
void generateUniformPts(int nums, cv::Mat& mask, std::vector<T>& outPts) {
    outPts.clear();
    cv::Scalar s = cv::sum(mask);
    if (s[0] == 0) {
        return;
    }
    cv::RNG rng(0);
    int times = 0;
    while (times < nums) {
        int x = rng.uniform(0, mask.cols);
        int y = rng.uniform(0, mask.rows);
        T currPt(x, y);
        uchar* pixel = mask.ptr<uchar>(y);
        if (pixel[x]) {
            outPts.push_back(currPt);
            times++;
        }
    }
}

// 对应OpenCV的cv::Mat转MATLAB uint8类型或logical或者double图像
template <typename T>
void convertCVToMatrix(const cv::Mat& srcImg, int rows, int cols, int channels, T* dst) {
    size_t elems = rows * cols;
    if (channels == 3) {
        cv::Mat channels[3];
        cv::split(srcImg.t(), channels);

        memcpy(dst, channels[2].data, elems * sizeof(T));              //copy channel[2] to the red channel
        memcpy(dst + elems, channels[1].data, elems * sizeof(T));      // green
        memcpy(dst + 2 * elems, channels[0].data, elems * sizeof(T));  // blue
    } else {
        cv::Mat dstImg = srcImg.t();
        memcpy(dst, dstImg.data, elems * sizeof(T));
    }
}

// 对应MATLAB uint8类型或者logical图像转cv::Mat，图像在内存中连续
template <typename T>
void convertToMatContinues(const T* inImg, int rows, int cols, int channels, cv::Mat& matBigImg) {
    size_t elems = (size_t)rows * cols;
    T* array = (T*)inImg;
    if (channels == 3) {
        cv::Mat matR = cv::Mat(cols, rows, CV_8UC1, array);  //inImg在内存中必须连续
        cv::Mat matG = cv::Mat(cols, rows, CV_8UC1, array + elems);
        cv::Mat matB = cv::Mat(cols, rows, CV_8UC1, array + 2 * elems);
        std::vector<cv::Mat> matBGR = {matB.t(), matG.t(), matR.t()};
        cv::merge(matBGR, matBigImg);
    } else {
        matBigImg = cv::Mat(cols, rows, CV_8UC1, (T*)inImg);
        matBigImg = matBigImg.t();
    }
}

class HDMapping {
   public:
    typedef enum { ORB_FEATURES,
                   LK_TRACK_FEATURES,
                   HYBRID_FEATURES } matchFeatureMethod;

    typedef enum { BUILD_MAP_SUCCESSFUL,
                   BUILD_MAP_FAILED,
                   BUILD_MAP_PROCESSING } buildMapStatus;

    typedef enum { LOCALIZE_MAP_SUCCESSFUL,
                   LOCALIZE_MAP_FAILED,
                   LOCALIZE_MAP_PROCESSING } localizeMapStatus;

   public:
    HDMapping();
    ~HDMapping();
    buildMapping::HDMapping::buildMapStatus constructWorldMap(const cv::Mat& srcImage, bool isStopConstructWorldMap = false, const char* vocFile = "database.yml.gz", const char* pointsFeatsFile = "pointsFeatures.bin", const char* mapFile = "hdMapCfg.json");

    buildMapping::HDMapping::localizeMapStatus localizeWorldMap(const cv::Mat& srcImage, const char* vocFile = "database.yml.gz", const char* pointsFeatsFile = "pointsFeatures.bin", const char* mapFile = "hdMapCfg.json");

   public:
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
    cv::Mat m_orbDetectMask, m_BW;  // 分别对应MATLAB中的vehiclePolygon所围成的mask，BW变量，分别用于检测ROI，更新的ROI
    std::vector<cv::KeyPoint> m_preKeypts, m_currKeypts;
    cv::Mat m_preDescriptors, m_currDescriptors;
    cv::Mat m_prevImg;
    cv::Mat m_preRelTform, m_relTform, m_previousImgPose;  // 2*3转换矩阵
    cv::Vec3d m_initViclePtPose, m_preViclePtPose;         // [x,y,theta]位姿
    std::vector<cv::Point2f> m_p0;

    //
    cv::Ptr<cv::Feature2D> m_orbDetector;
    matchFeatureMethod m_method;
    void selectUniformPoints(std::vector<cv::KeyPoint>& keyPoints, int numRetPoints,
                             cv::Size size, std::vector<cv::KeyPoint>& outputPts, std::vector<int>& indexs);

    void estiTform(std::vector<cv::Point2f>& prePoints, std::vector<cv::Point2f>& currPoints, cv::Mat& tform2x3, cv::Mat& inliers, int& status);

    void fuseOptimizeHDMap(std::string imageFilesList, std::vector<cv::Vec3d>& updateNodeVehiclePtPoses);

    void detectLanes(cv::Mat& grayImage, std::vector<cv::line_descriptor::KeyLine>& lines);

    // loop
    std::vector<imageKptsAndFeatures> m_points_features;
    DBoW3::Database m_db;
    // SlamGraph2D::slamPoseGraph m_pg;  // 索引从1开始
    void detectLoopAndAddGraph(std::string vocabularyFile = "./database.yml.gz");
    void optimizePoseGraph(cv::Mat& loopIDpairs, std::vector<cv::Vec3d>& relPoses);
    void loopDatabaseAddFeaturesAndSave(std::string saveVocabularyFile);
    DBoW3::QueryResults retrieveImage(cv::Mat queryImage, int topK);

    void reset();
    void saveMapData();
    void loadMapData();
    void drawRoutePath(cv::Mat& drawImage) const;
};

}  // namespace buildMapping
