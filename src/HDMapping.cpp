#include <fstream>
#include <algorithm>
#include "HDMapping.h"

#define DEBUG_SHOW_REALTIME_IMAGE 0
#define DEBUG_SHOW_FUSE_OPTIMIZE_IMAGE 1

namespace buildMapping {
class HDMapping;
}

buildMapping::HDMapping::HDMapping() {
    m_cumDist = 0.0;
    m_pixelExtentInWorldXY = 0.015 / 0.51227;
    m_isBuildMap = true;
    m_buildMapStopFrame = 0;
    m_isBuildMapOver = false;
    m_isLocSuccess = false;
    m_locVehiclePose[0] = 0;
    m_locVehiclePose[1] = 0;
    m_locVehiclePose[2] = 0;

    m_orbDetectMask = 255 * cv::Mat::ones(480, 640, CV_8UC1);
    m_orbDetectMask.rowRange(156, 326).colRange(285, 358) = 0;
    m_BW = 255 * cv::Mat::ones(480, 640, CV_8UC1);
    m_BW.rowRange(480 / 2, 480) = 0;
    m_BW = m_BW &= m_orbDetectMask;
    m_initViclePtPose = cv::Vec3f((285 + 358) / 2.0, (156, 326) / 2.0, 0);
    m_method = matchFeatureMethod::LK_TRACK_FEATURES;

    //第一副图像的像素坐标系为世界坐标系
    m_preRelTform = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                     0, 1, 0);
    m_relTform = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                  0, 1, 0);
    m_previousImgPose = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                         0, 1, 0);
    m_orbDetector = cv::ORB::create(307200);
    m_preViclePtPose = m_initViclePtPose;

    // loop
    // m_pg.init(10000, 5000);
};
buildMapping::HDMapping::~HDMapping() {}

void buildMapping::HDMapping::reset() {
    m_prevImg = cv::Mat();
    m_preKeypts.clear();
    m_preDescriptors = cv::Mat();

    m_cumDist = 0.0;
    m_vehiclePoses.clear();
    m_points_features.clear();

    m_HDmapOutput.bigImg = cv::Mat();
    m_HDmapOutput.ref.XWorldLimits[0] = 0;
    m_HDmapOutput.ref.XWorldLimits[1] = 1;
    m_HDmapOutput.ref.YWorldLimits[0] = 0;
    m_HDmapOutput.ref.YWorldLimits[1] = 1;

    m_isBuildMap = true;
    m_isBuildMapOver = false;
    m_isLocSuccess = false;
    m_locVehiclePose[0] = 0;
    m_locVehiclePose[1] = 0;
    m_locVehiclePose[2] = 0;

    m_preRelTform = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                     0, 1, 0);
    m_relTform = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                  0, 1, 0);
    m_previousImgPose = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                         0, 1, 0);
    m_orbDetector = cv::ORB::create(307200);
}

void buildMapping::HDMapping::saveMapData() {
}

void buildMapping::HDMapping::loadMapData() {
}

buildMapping::HDMapping::buildMapStatus buildMapping::HDMapping::constructWorldMap(const cv::Mat& srcImage,
                                                                                   bool isStopConstructWorldMap, const char* vocFile, const char* pointsFeatsFile, const char* mapFile) {
    if (srcImage.empty()) {
        return buildMapping::HDMapping::buildMapStatus::BUILD_MAP_FAILED;
    }
    static size_t num = 0;

    cv::Mat currImg = srcImage;
    if (currImg.channels() == 3)
        cv::cvtColor(currImg, currImg, cv::COLOR_BGR2GRAY);

    double t1 = cv::getTickCount();
    std::vector<cv::KeyPoint> keyPts;
    m_orbDetector->detect(currImg, keyPts, m_orbDetectMask);
    std::vector<int> indexs;
    selectUniformPoints(keyPts, 2000, cv::Size(currImg.cols, currImg.rows), m_currKeypts, indexs);
    if (m_preKeypts.empty()) {
        m_prevImg = currImg;
        m_preKeypts = m_currKeypts;

        m_HDmapOutput.bigImg = currImg;
        m_HDmapOutput.ref.XWorldLimits[0] = 0;
        m_HDmapOutput.ref.XWorldLimits[1] = currImg.cols - 1;
        m_HDmapOutput.ref.YWorldLimits[0] = 0;
        m_HDmapOutput.ref.YWorldLimits[1] = currImg.rows - 1;
    }
    printf("Detect Elapsed second Time:%.6f\n", (cv::getTickCount() - t1) * 1.0 / cv::getTickFrequency());

    t1 = cv::getTickCount();
    std::vector<cv::DMatch> matches;
    std::vector<cv::Point2f> preP, nextP;
    cv::Mat inliers;
    int status = 0;
    int rigidTformType = 0;
    m_orbDetector->compute(currImg, m_currKeypts, m_currDescriptors);
    if (m_preDescriptors.empty()) {
        m_preDescriptors = m_currDescriptors;
    }
    m_points_features.push_back({m_currKeypts, m_currDescriptors.clone()});

    printf("compute Elapsed second Time:%.6f\n", (cv::getTickCount() - t1) * 1.0 / cv::getTickFrequency());

    // detect lanes
    t1 = cv::getTickCount();
    std::vector<cv::line_descriptor::KeyLine> lines;
    detectLanes(currImg, lines);
    printf("detectLanes Elapsed second Time:%.6f\n",
           (cv::getTickCount() - t1) * 1.0 / cv::getTickFrequency());
    cv::Mat outImg1;
    cv::line_descriptor::drawKeylines(currImg, lines, outImg1, cv::Scalar(0, 0, 255));
    cv::imwrite("keyline.jpg", outImg1);

    t1 = cv::getTickCount();
    if (m_method == buildMapping::HDMapping::matchFeatureMethod::HYBRID_FEATURES) {
        double t2 = cv::getTickCount();
        //match features
        // https://stackoverflow.com/questions/43830849/opencv-use-flann-with-orb-descriptors-to-match-features
        // cv::FlannBasedMatcher matcher;  // = cv::FlannBasedMatcher(cv::makePtr<cv::flann::LshIndexParams>(20, 10, 2));
        cv::BFMatcher matcher(cv::NORM_HAMMING);
        cv::Mat feature1, feature2;
        m_preDescriptors.convertTo(feature1, CV_8U);
        m_currDescriptors.convertTo(feature2, CV_8U);

#if 0
        matcher.match(feature1, feature2, matches);
        //只初步选取匹配前num个较好的特征点
        int numMatches = matches.size() < 10 ? matches.size() : matches.size() / 2;
        std::nth_element(matches.begin(), matches.begin() + numMatches, matches.end());
        matches.erase(matches.begin() + numMatches, matches.end());
#else
        std::vector<std::vector<cv::DMatch>> matchePoints;
        matcher.knnMatch(feature1, feature2, matchePoints, 2);
        for (int i = 0; i < matchePoints.size(); i++) {
            if (matchePoints[i][0].distance < 0.5 * matchePoints[i][1].distance) {
                matches.push_back(matchePoints[i][0]);
            }
        }
#endif
        printf("match Elapsed second Time:%.6f\n", (cv::getTickCount() - t2) * 1.0 / cv::getTickFrequency());
        // cv::Mat drawDstImage;
        // cv::drawMatches(m_prevImg, m_preKeypts, currImg, m_currKeypts, matches, drawDstImage, cv::Scalar(0, 255, 0), -1, std::vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
        // cv::imwrite("drawDstImage.jpg", drawDstImage);
        double t3 = cv::getTickCount();
        // estimate geometry rigid 2D transformation matrix
#ifdef USE_OPENCV_FUNCTION
        std::vector<cv::Point2f>
            points1, points2;
        for (size_t i = 0; i < matches.size(); i++) {
            points1.push_back(preKeypts[matches[i].queryIdx].pt);
            points2.push_back(currKeypts[matches[i].trainIdx].pt);
        }
        cv::Mat inliers = cv::Mat::zeros(matches.size(), 1, CV_8U);
        cv::Mat optimalAffineMat = estimateAffinePartial2D(points1, points2, inliers, cv::RANSAC);
        cv::Mat R = optimalAffineMat.rowRange(0, 2).colRange(0, 2);
        double s = std::sqrt(cv::determinant(R));
        cv::Mat rigidtform2dR = R.mul(1.0 / s);
        cv::hconcat(rigidtform2dR, optimalAffineMat.col(2), m_relTform);
        status = std::abs(s - 1.0) / 1.0 < 0.03 ? 0 : 1;
#else
        preP.resize(matches.size());
        nextP.resize(matches.size());
        for (size_t i = 0; i < matches.size(); i++) {
            preP[i] = m_preKeypts[matches[i].queryIdx].pt;
            nextP[i] = m_currKeypts[matches[i].trainIdx].pt;
        }
        estiTform(preP, nextP, m_relTform, inliers, status);
#endif
        printf("estimate Elapsed second Time:%.6f\n", (cv::getTickCount() - t3) * 1.0 / cv::getTickFrequency());
        // std::cout << "status:" << status << ",sum:" << cv::sum(inliers) << ",matches.size():" << matches.size() << std::endl;
        double t4 = cv::getTickCount();
        bool cond = (status > 0 || cv::sum(inliers)[0] <= 3);
        if (cond) {
            rigidTformType = 1;
            std::vector<uchar> statusOK;
            std::vector<float> err;
            cv::TermCriteria criteria = cv::TermCriteria((cv::TermCriteria::COUNT) +
                                                             (cv::TermCriteria::EPS),
                                                         10, 0.03);
            std::vector<cv::Point2f> p0, p1;
            cv::KeyPoint::convert(m_preKeypts, p0);
            cv::calcOpticalFlowPyrLK(m_prevImg, currImg, p0, p1, statusOK, err, cv::Size(15, 15), 2, criteria);
            preP.clear();
            nextP.clear();
            matches.clear();  // debug
            for (uint i = 0; i < p0.size(); i++) {
                if (statusOK[i] == 1) {
                    matches.push_back(cv::DMatch(preP.size(), nextP.size(), 1));
                    preP.push_back(p0[i]);
                    nextP.push_back(p1[i]);
                }
            }
            estiTform(preP, nextP, m_relTform, inliers, status);
            cond = (status > 0 || cv::sum(inliers)[0] <= 3);
            if (cond) {
                rigidTformType = 2;
                m_relTform = m_preRelTform;
            }
        }
        printf("cond Elapsed second Time:%.6f\n", (cv::getTickCount() - t4) * 1.0 / cv::getTickFrequency());
    } else if (m_method == buildMapping::HDMapping::matchFeatureMethod::LK_TRACK_FEATURES) {
        rigidTformType = 1;
        std::vector<uchar> statusOK;
        std::vector<float> err;
        cv::TermCriteria criteria = cv::TermCriteria((cv::TermCriteria::COUNT) + (cv::TermCriteria::EPS), 10, 0.03);
        std::vector<cv::Point2f> p0, p1;
        cv::KeyPoint::convert(m_preKeypts, p0);
        double t21 = cv::getTickCount();
        cv::calcOpticalFlowPyrLK(m_prevImg, currImg, p0, p1, statusOK, err, cv::Size(15, 15), 2, criteria);
        printf("LK Elapsed second Time:%.6f\n", (cv::getTickCount() - t21) * 1.0 / cv::getTickFrequency());
        preP.clear();
        nextP.clear();
#if DEBUG_SHOW_REALTIME_IMAGE
        matches.clear();  // debug
#endif
        for (uint i = 0; i < p0.size(); i++) {
            if (statusOK[i] == 1) {
#if DEBUG_SHOW_REALTIME_IMAGE
                matches.push_back(cv::DMatch(preP.size(), nextP.size(), 1));
#endif
                preP.push_back(p0[i]);
                nextP.push_back(p1[i]);
            }
        }
        double t22 = cv::getTickCount();

        // cv::Mat inliers = cv::Mat::zeros(matches.size(), 1, CV_8U);
        cv::Mat optimalAffineMat = estimateAffinePartial2D(preP, nextP, inliers, cv::RANSAC, 1.5, 2000, 0.99, 0);
        std::vector<cv::Point2f> refineP0, refineP1;
        for (size_t i = 0; i < inliers.rows; i++) {
            if (inliers.at<uchar>(0, i) == 1) {
                refineP0.push_back(preP[i]);
                refineP1.push_back(nextP[i]);
            }
        }
        cv::Mat inliersIner;
        estiTform(refineP0, refineP1, m_relTform, inliersIner, status);
        int innerIdx = 0;
        for (size_t i = 0; i < inliers.rows; i++) {
            if (inliers.at<uchar>(0, i) == 1) {
                if (inliersIner.at<uchar>(0, innerIdx) == 0) {
                    inliers.at<uchar>(0, i) = 0;
                }
                innerIdx++;
            }
        }

        printf("estiTform Elapsed second Time:%.6f\n", (cv::getTickCount() - t22) * 1.0 / cv::getTickFrequency());
        bool cond = (status > 0 || cv::sum(inliersIner)[0] <= 3);
        if (cond) {
            rigidTformType = 2;
            m_relTform = m_preRelTform;
#if DEBUG_SHOW_REALTIME_IMAGE
            std::ofstream fid("out_may_be_lost_frame.txt", std::ios_base::app);
            fid << num << std::endl;
            fid.close();
#endif
            saveMapData();
            // reset();
        }
    } else {
        //match features
        // cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create("FlannBased");
        cv::FlannBasedMatcher matcher;
        cv::Mat feature1, feature2;
        if (m_preDescriptors.type() != CV_32F) {
            m_preDescriptors.convertTo(feature1, CV_32F);
        }
        if (m_currDescriptors.type() != CV_32F) {
            m_currDescriptors.convertTo(feature2, CV_32F);
        }
#if 0
        matcher.match(feature1,feature2, matches);
        //只初步选取匹配前num个较好的特征点
        int numMatches = matches.size() < 10 ? matches.size() : matches.size() / 2;
        std::nth_element(matches.begin(), matches.begin() + numMatches, matches.end());
        matches.erase(matches.begin() + numMatches, matches.end());
#else
        std::vector<std::vector<cv::DMatch>> matchePoints;
        matcher.knnMatch(feature1, feature2, matchePoints, 2);
        for (int i = 0; i < matchePoints.size(); i++) {
            if (matchePoints[i][0].distance < 0.6 * matchePoints[i][1].distance) {
                matches.push_back(matchePoints[i][0]);
            }
        }
#endif
        // estimate geometry rigid 2D transformation matrix
        preP.resize(matches.size());
        nextP.resize(matches.size());
        for (size_t i = 0; i < matches.size(); i++) {
            preP[i] = m_preKeypts[matches[i].queryIdx].pt;
            nextP[i] = m_currKeypts[matches[i].trainIdx].pt;
        }
        estiTform(preP, nextP, m_relTform, inliers, status);
    }
    printf("  (match+esti+cond)Track_Features Elapsed second Time:%.6f\n",
           (cv::getTickCount() - t1) * 1.0 / cv::getTickFrequency());

    t1 = cv::getTickCount();
    // build map mode
    cv::Mat previousImgPoseA = m_previousImgPose;
    cv::Mat homoMatrix = (cv::Mat_<double>(1, 3) << 0, 0, 1);
    previousImgPoseA.push_back(homoMatrix);
    cv::Mat currImgPose = m_relTform * previousImgPoseA;

    cv::Mat tform,
        tempR, tempT;
    tempR = currImgPose.rowRange(0, 2).colRange(0, 2).t();
    tempT = -tempR * currImgPose.col(2);
    cv::hconcat(tempR, tempT, tform);

    // calcuate vehicle poses
    cv::Mat initViclePtPoseA = (cv::Mat_<double>(3, 3) << 1, 0, m_initViclePtPose[0],
                                0, 1, m_initViclePtPose[1],
                                0, 0, 1);
    cv::Mat currVehiclePoseMat = tform * initViclePtPoseA;
    tempR = (cv::Mat_<double>(3, 3) << currVehiclePoseMat.at<double>(0, 0), currVehiclePoseMat.at<double>(0, 1), 0,
             currVehiclePoseMat.at<double>(1, 0), currVehiclePoseMat.at<double>(1, 1), 0,
             0, 0, 1);
    std::vector<double> angRadius;
    cv::Rodrigues(tempR, angRadius);
    cv::Vec3d currVehiclePtPose = cv::Vec3d(currVehiclePoseMat.at<double>(0, 2), currVehiclePoseMat.at<double>(1, 2), angRadius[2]);
    m_vehiclePoses.push_back(currVehiclePtPose);
    if (m_vehiclePoses.size() >= 2) {
        cv::Vec3d prePose = m_vehiclePoses[m_vehiclePoses.size() - 2];
        cv::Vec3d currPose = m_vehiclePoses[m_vehiclePoses.size() - 1];
        m_cumDist += m_pixelExtentInWorldXY * std::sqrt((currPose[0] - prePose[0]) * (currPose[0] - prePose[0]) + (currPose[1] - prePose[1]) * (currPose[1] - prePose[1]));
    }

    cv::Vec3d meas = currVehiclePtPose - m_preViclePtPose;
    double measurement[3] = {meas[0], meas[1], meas[2]};
    double fromNodeID = num + 1;
    double toNodeID = num + 2;
    // m_pg.addRelativePose2D(measurement, fromNodeID, toNodeID);

    // 开始poseGraph优化建图
    if (isStopConstructWorldMap) {
        m_isBuildMap = false;
        m_isBuildMapOver = true;
        m_buildMapStopFrame = num;
        detectLoopAndAddGraph(std::string(vocFile));

// 融合姿态图输出
#if DEBUG_SHOW_FUSE_OPTIMIZE_IMAGE
        std::string imageFilesList = "imageFilesList.txt";
        fuseOptimizeHDMap(imageFilesList, m_vehiclePoses);
        cv::Mat drawImage;
        drawRoutePath(drawImage);
        cv::imwrite("bigImgCopy_fuseOptimize.jpg", drawImage);
#endif

        // 保存特征和状态相关数据
        std::ofstream os(pointsFeatsFile, std::ios::binary);  // 打开标准输出流
        cereal::PortableBinaryOutputArchive archiveBin(os);   // 构建cereal对象，并用os初始化
        archiveBin(m_points_features);

        std::ofstream file(mapFile);
        cereal::JSONOutputArchive archiveJson(file);
        std::vector<cerealPoses> vehiclePoses;
        for (size_t i = 0; i < m_vehiclePoses.size(); i++) {
            vehiclePoses.push_back(cerealPoses{m_vehiclePoses[i][0], m_vehiclePoses[i][1], m_vehiclePoses[i][2]});
        }
        archiveJson(CEREAL_NVP(m_cumDist), CEREAL_NVP(m_pixelExtentInWorldXY), CEREAL_NVP(m_isBuildMap),
                    CEREAL_NVP(m_isBuildMapOver), CEREAL_NVP(m_buildMapStopFrame), CEREAL_NVP(m_isLocSuccess),
                    CEREAL_NVP(m_locVehiclePose), CEREAL_NVP(vehiclePoses));

        return buildMapStatus::BUILD_MAP_SUCCESSFUL;
    }
    printf("calcuate dist and add pose:%.6f\n", (cv::getTickCount() - t1) * 1.0 / cv::getTickFrequency());
    t1 = cv::getTickCount();

#if DEBUG_SHOW_REALTIME_IMAGE
    cv::Mat currCorner = (cv::Mat_<double>(3, 4) << 0, currImg.cols - 1, currImg.cols - 1, 0,
                          0, 0, currImg.rows - 1, currImg.rows - 1,
                          1, 1, 1, 1);
    currCorner = tform * currCorner;

    double xmin, xmax, ymin, ymax;
    cv::minMaxLoc(currCorner.row(0), &xmin, &xmax);
    cv::minMaxLoc(currCorner.row(1), &ymin, &ymax);
    double pad_left = 0, pad_right = 0, pad_top = 0, pad_down = 0;
    if (xmin < m_HDmapOutput.ref.XWorldLimits[0]) {
        pad_left = m_HDmapOutput.ref.XWorldLimits[0] - xmin;
        m_HDmapOutput.ref.XWorldLimits[0] = xmin;
    }

    if (xmax > m_HDmapOutput.ref.XWorldLimits[1]) {
        pad_right = xmax - m_HDmapOutput.ref.XWorldLimits[1];
        m_HDmapOutput.ref.XWorldLimits[1] = xmax;
    }

    if (ymin < m_HDmapOutput.ref.YWorldLimits[0]) {
        pad_top = m_HDmapOutput.ref.YWorldLimits[0] - ymin;
        m_HDmapOutput.ref.YWorldLimits[0] = ymin;
    }

    if (ymax > m_HDmapOutput.ref.YWorldLimits[1]) {
        pad_down = ymax - m_HDmapOutput.ref.YWorldLimits[1];
        m_HDmapOutput.ref.YWorldLimits[1] = ymax;
    }

    m_HDmapOutput.ref.ImageSize[0] = std::round(m_HDmapOutput.ref.YWorldLimits[1] - m_HDmapOutput.ref.YWorldLimits[0]) + 1;
    m_HDmapOutput.ref.ImageSize[1] = std::round(m_HDmapOutput.ref.XWorldLimits[1] - m_HDmapOutput.ref.XWorldLimits[0]) + 1;

    // 平移到可视化区域
    cv::Mat topImg, bwMask;
    tform.colRange(2, 3) = tform.colRange(2, 3) - (cv::Mat_<double>(2, 1) << m_HDmapOutput.ref.XWorldLimits[0], m_HDmapOutput.ref.YWorldLimits[0]);
    cv::warpAffine(currImg, topImg, tform, cv::Size(m_HDmapOutput.ref.ImageSize[1], m_HDmapOutput.ref.ImageSize[0]));
    cv::warpAffine(m_BW, bwMask, tform, cv::Size(m_HDmapOutput.ref.ImageSize[1], m_HDmapOutput.ref.ImageSize[0]));
    cv::copyMakeBorder(m_HDmapOutput.bigImg, m_HDmapOutput.bigImg, std::round(pad_top), std::round(pad_down), std::round(pad_left), std::round(pad_right), cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::resize(m_HDmapOutput.bigImg, m_HDmapOutput.bigImg, cv::Size(topImg.cols, topImg.rows));
    topImg.copyTo(m_HDmapOutput.bigImg, bwMask);
    std::vector<cv::DMatch> rigidMatches;
    for (size_t i = 0; i < matches.size(); i++) {
        if (inliers.at<uchar>(i, 0)) {
            rigidMatches.push_back(matches[i]);
        }
    }
    std::cout << "matches ratio:" << matches.size() << "/" << m_preKeypts.size() << std::endl;
    std::cout << "rigid estimate matches ratio:" << rigidMatches.size() << "/" << matches.size() << ",rigidTformType:" << rigidTformType << std::endl;

    // assert(rigidMatches.size() > 3);

    // Draw top matches
    cv::Mat imMatches;
    if (rigidTformType == 0) {
        cv::drawMatches(m_prevImg, m_preKeypts, currImg, m_currKeypts, rigidMatches, imMatches, cv::Scalar(0, 255, 255), -1, std::vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
    } else {
        std::vector<cv::KeyPoint> p1_temp, p2_temp;
        cv::KeyPoint::convert(preP, p1_temp, 1, 1, 0, -1);
        cv::KeyPoint::convert(nextP, p2_temp, 1, 1, 0, -1);
        cv::drawMatches(m_prevImg, p1_temp, currImg, p2_temp, rigidMatches, imMatches, cv::Scalar(0, 0, 255), -1, std::vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
    }
    cv::imwrite("matches_" + std::to_string(rigidTformType) + ".jpg", imMatches);
    // cv::Mat dstMat;
    // srcImage.copyTo(dstMat, orbDetectMask);
    cv::Mat drawImage;
    drawRoutePath(drawImage);
    cv::imwrite("bigImgCopy.jpg", drawImage);
#endif

    // update previous state variables
    m_prevImg = currImg;
    m_preKeypts = m_currKeypts;
    m_preDescriptors = m_currDescriptors.clone();
    m_previousImgPose = currImgPose;
    m_preRelTform = m_relTform;
    m_preViclePtPose = currVehiclePtPose;

    num++;
    printf("BUILDMAP Elapsed second Time:%.6f\n\n", (cv::getTickCount() - t1) * 1.0 / cv::getTickFrequency());
    return buildMapStatus::BUILD_MAP_PROCESSING;
}

buildMapping::HDMapping::localizeMapStatus buildMapping::HDMapping::localizeWorldMap(const cv::Mat& srcImage, const char* vocFile, const char* pointsFeatsFile, const char* mapFile) {
    static bool isLoadMapData = false;
    if (srcImage.empty() || !filesystem::path(vocFile).exists() || !filesystem::path(pointsFeatsFile).exists() || !filesystem::path(mapFile).exists()) {
        std::cout << "current path:" << filesystem::path(".").make_absolute() << std::endl;
        return buildMapping::HDMapping::localizeMapStatus::LOCALIZE_MAP_FAILED;
    }

    if (!isLoadMapData) {
        m_db.load(std::string(vocFile));

        std::ifstream is(std::string(pointsFeatsFile), std::ios::binary);
        cereal::PortableBinaryInputArchive iarchive(is);
        iarchive(m_points_features);

        std::ifstream file(mapFile);
        cereal::JSONInputArchive archive1(file);
        std::vector<cerealPoses> vehiclePoses;
        archive1(CEREAL_NVP(m_cumDist), CEREAL_NVP(m_pixelExtentInWorldXY), CEREAL_NVP(m_isBuildMap),
                 CEREAL_NVP(m_isBuildMapOver), CEREAL_NVP(m_buildMapStopFrame), CEREAL_NVP(m_isLocSuccess),
                 CEREAL_NVP(vehiclePoses));  //, CEREAL_NVP(m_vehiclePoses)
        for (size_t i = 0; i < vehiclePoses.size(); i++) {
            m_vehiclePoses.push_back(cv::Vec3d(vehiclePoses[i].x, vehiclePoses[i].y, vehiclePoses[i].theta));
        }
        isLoadMapData = true;
    }

    cv::Mat currImg = srcImage;
    int topK = 20;
    DBoW3::QueryResults result = retrieveImage(currImg, topK);
    DBoW3::QueryResults::const_iterator qit;
    size_t rowIdx = 0;
    int imageIDs[topK] = {0};
    double scores[topK] = {0.0};
    for (qit = result.begin(); qit != result.end(); ++qit) {
        imageIDs[rowIdx] = (double)qit->Id;   // 注意C++索引从0开始
        scores[rowIdx] = (double)qit->Score;  // 注意前面程序用的是topK
        rowIdx++;
    }

    bool isDetectedLoop = false;
    int loopCandidate = 0;
    double minScore = *std::min_element(scores, scores + topK);
    double bestScore = scores[1];
    double ratio = 0.6;           // 根据实际图像调整此值
    double expeThreshold = 0.18;  // 0.2根据数据集经验取得, 见doc/loopClosureDetect.md
    double threshold = std::max({bestScore * ratio, minScore, expeThreshold});

    std::vector<int> validIdxs = findItems(std::vector<double>(scores, scores + topK), threshold);
    std::vector<int> validKeyFrameIds;
    for (size_t validId = 0; validId < validIdxs.size(); validId++) {
        validKeyFrameIds.push_back(imageIDs[validId]);
    }
    if (validKeyFrameIds.size() > 2) {
        std::vector<std::vector<int>> groups = nchoosek(validKeyFrameIds, 3);
        cv::Mat matGroups(groups.size(), groups.at(0).size(), CV_64FC1);
        for (int iRow = 0; iRow < matGroups.rows; ++iRow) {
            double* ele = matGroups.ptr<double>(iRow);
            for (int jCol = 0; jCol < matGroups.cols; ++jCol) {
                ele[jCol] = groups[iRow][jCol];
            }
        }

        for (size_t i_r = 0; i_r < matGroups.rows; i_r++) {
            double minV, maxV;
            cv::minMaxLoc(matGroups.row(i_r), &minV, &maxV, NULL, NULL);
            if (maxV - minV <= 3) {
                cv::Mat consecutiveGroups = matGroups.row(i_r);
                loopCandidate = (int)consecutiveGroups.at<double>(0, 0);
                isDetectedLoop = true;
                break;
            } else {
                isDetectedLoop = false;
            }
        }
    }

    if (isDetectedLoop) {
        std::vector<cv::DMatch> matches;
        std::vector<cv::Point2f> preP, nextP;
        cv::Mat inliers;
        int status = 0;
        cv::BFMatcher matcher(cv::NORM_HAMMING);
        cv::Mat feature1 = m_points_features[loopCandidate].features;
        cv::Mat feature2 = m_currDescriptors;
        std::vector<std::vector<cv::DMatch>> matchePoints;
        matcher.knnMatch(feature1, feature2, matchePoints, 2);
        for (int i = 0; i < matchePoints.size(); i++) {
            if (matchePoints[i][0].distance < 0.6 * matchePoints[i][1].distance) {
                matches.push_back(matchePoints[i][0]);
            }
        }

        preP.resize(matches.size());
        nextP.resize(matches.size());
        for (size_t i = 0; i < matches.size(); i++) {
            preP[i] = m_points_features[loopCandidate].keyPoints[matches[i].queryIdx].pt;
            nextP[i] = m_currKeypts[matches[i].trainIdx].pt;
        }
        estiTform(preP, nextP, m_relTform, inliers, status);
        bool cond = (status > 0 || cv::sum(inliers)[0] <= 3);
        if (cond) {
            return buildMapping::HDMapping::localizeMapStatus::LOCALIZE_MAP_PROCESSING;
        }

        // calcuate localize vehicle absolute pose
        cv::Vec3d localizeBasePose = m_vehiclePoses[loopCandidate];
        cv::Mat localizeBasePoseA = (cv::Mat_<double>(3, 3) << std::cos(localizeBasePose[2]), -std::sin(localizeBasePose[2]), localizeBasePose[0],
                                     std::sin(localizeBasePose[2]), std::cos(localizeBasePose[2]), localizeBasePose[1],
                                     0, 0, 1);
        cv::Mat tform1 = (cv::Mat_<double>(3, 3) << 1, 0, -m_initViclePtPose[0],
                          0, 1, -m_initViclePtPose[1],
                          0, 0, 1);
        cv::Mat localizeBaseImagePoseA = tform1 * localizeBasePoseA;
        cv::Mat relTformA = m_relTform;
        cv::Mat temp = (cv::Mat_<double>(1, 3) << 0, 0, 1);
        relTformA.push_back(temp);
        cv::Mat currImagePoseA = relTformA * localizeBaseImagePoseA;
        cv::Mat tform2 = (cv::Mat_<double>(3, 3) << 1, 0, m_initViclePtPose[0],
                          0, 1, m_initViclePtPose[1],
                          0, 0, 1);
        cv::Mat currLocalizeVehiclePoseA = tform2 * currImagePoseA;

        cv::Mat tempR = (cv::Mat_<double>(3, 3) << currLocalizeVehiclePoseA.at<double>(0, 0), currLocalizeVehiclePoseA.at<double>(0, 1), 0,
                         currLocalizeVehiclePoseA.at<double>(1, 0), currLocalizeVehiclePoseA.at<double>(1, 1), 0,
                         0, 0, 1);
        std::vector<double> angRadius;
        cv::Rodrigues(tempR, angRadius);

        m_locVehiclePose[0] = currLocalizeVehiclePoseA.at<double>(0, 2);
        m_locVehiclePose[1] = currLocalizeVehiclePoseA.at<double>(1, 2);
        m_locVehiclePose[2] = angRadius[2];

        m_isBuildMap = false;
        m_isLocSuccess = true;
        return buildMapping::HDMapping::localizeMapStatus::LOCALIZE_MAP_SUCCESSFUL;
    } else {
        m_isLocSuccess = false;
        return buildMapping::HDMapping::localizeMapStatus::LOCALIZE_MAP_PROCESSING;
    }
}

void buildMapping::HDMapping::estiTform(std::vector<cv::Point2f>& prePoints, std::vector<cv::Point2f>& currPoints,
                                        cv::Mat& tform2x3, cv::Mat& inliers, int& status) {
    coder::array<double, 2U> pts1_tmp, pts2_tmp;
    coder::array<boolean_T, 2U> inlierIndex;
    double tform2x3_[6];

    pts1_tmp.set_size(prePoints.size(), 2);
    pts2_tmp.set_size(currPoints.size(), 2);
    for (size_t i = 0; i < currPoints.size(); i++) {
        pts1_tmp[i] = prePoints[i].x;
        pts1_tmp[i + currPoints.size()] = prePoints[i].y;

        pts2_tmp[i] = currPoints[i].x;
        pts2_tmp[i + currPoints.size()] = currPoints[i].y;
    }
    estimateAffineRigid2D::estimateAffineRigid2D(pts1_tmp, pts2_tmp, tform2x3_, inlierIndex, &status);
    inliers = cv::Mat(prePoints.size(), 1, CV_8U, inlierIndex.data()).clone();
    tform2x3 = (cv::Mat_<double>(2, 3) << tform2x3_[0], tform2x3_[2], tform2x3_[4], tform2x3_[1], tform2x3_[3], tform2x3_[5]);
}

void buildMapping::HDMapping::selectUniformPoints(std::vector<cv::KeyPoint>& keyPoints, int numRetPoints,
                                                  cv::Size size, std::vector<cv::KeyPoint>& outputPts, std::vector<int>& indexs) {
    outputPts.clear();
    indexs.clear();
    if (numRetPoints >= keyPoints.size()) {
        outputPts = keyPoints;
        indexs.reserve(keyPoints.size());
        std::iota(indexs.begin(), indexs.end(), 0);
        return;
    }

    coder::array<double, 2U> points;
    coder::array<double, 2U> pointsOut;
    coder::array<double, 1U> b_index;
    coder::array<double, 1U> responses;
    double dv[2];

    points.set_size(keyPoints.size(), 2);
    responses.set_size(keyPoints.size());
    for (size_t i = 0; i < points.size(0); i++) {
        points[i] = keyPoints[i].pt.x;
        points[i + points.size(0)] = keyPoints[i].pt.y;
        responses[i] = keyPoints[i].response;
    }

    dv[0] = size.height;
    dv[1] = size.width;
    selectUniform2::selectUniform2(points, responses, (double)numRetPoints, dv, pointsOut, b_index);

    for (size_t i = 0; i < pointsOut.size(0); i++) {
        indexs.push_back(b_index[i] - 1);
        outputPts.push_back(keyPoints[(int)(b_index[i]) - 1]);
    }
}

void buildMapping::HDMapping::loopDatabaseAddFeaturesAndSave(std::string saveDataBaseYmlGz) {
    // branching factor and depth levels
    const int k = 10;
    const int L = 4;
    const DBoW3::WeightingType weight = DBoW3::TF_IDF;
    const DBoW3::ScoringType score = DBoW3::L1_NORM;

    std::cout << "From features,Create vocabulary,please wait ..." << std::endl;
    DBoW3::Vocabulary voc(k, L, weight, score);
    std::vector<cv::Mat> features;
    for (size_t i = 0; i < m_points_features.size(); i++) {
        features.push_back(m_points_features[i].features);
    }
    voc.create(features);

    m_db.setVocabulary(voc, false, 0);  // false = do not use direct index
    // (so ignore the last param)
    // The direct index is useful if we want to retrieve the features that
    // belong to some vocabulary node.
    // db creates a copy of the vocabulary, we may get rid of "voc" now

    // add images to the database
    for (size_t i = 0; i < features.size(); i++)
        m_db.add(features[i]);

    std::cout << "add features per image done!" << std::endl;

    std::string databaseFile = saveDataBaseYmlGz;
    std::cout << "Vocabulary information: " << std::endl
              << voc << std::endl
              << "have saved this path:" << databaseFile << std::endl;
    m_db.save(databaseFile);
}

DBoW3::QueryResults buildMapping::HDMapping::retrieveImage(cv::Mat queryImage, int topK) {
    // add images to the database
    cv::Ptr<cv::Feature2D> fdetector = cv::ORB::create(480 * 640);
    std::vector<cv::KeyPoint> keypoints, outkeypoints;
    cv::Mat descriptors;
    if (queryImage.channels() == 3) {
        cv::cvtColor(queryImage, queryImage, cv::COLOR_BGR2GRAY);
    }

    if (queryImage.empty()) throw std::runtime_error("Could not open image");
    fdetector->detect(queryImage, keypoints, m_orbDetectMask);
    std::vector<int> indexs;
    selectUniformPoints(keypoints, 2000, cv::Size(queryImage.cols, queryImage.rows), outkeypoints, indexs);
    fdetector->compute(queryImage, outkeypoints, descriptors);

    m_currDescriptors = descriptors;
    m_currKeypts = outkeypoints;

    DBoW3::QueryResults ret;
    m_db.query(descriptors, ret, topK);  // m_db是已经加入了每幅图像的特征的database，此处选取的是top 20

    return ret;
}

void buildMapping::HDMapping::detectLoopAndAddGraph(std::string vocabularyFile) {
    std::string saveDataBaseYmlGz = vocabularyFile;
    if (filesystem::path(saveDataBaseYmlGz).exists()) {
        std::cout << "loading database,please wait... ,name:" << saveDataBaseYmlGz << std::endl;
        m_db.load(saveDataBaseYmlGz);
        std::cout << "load finished!" << std::endl;
    } else {
        loopDatabaseAddFeaturesAndSave(saveDataBaseYmlGz);
    }

    // multiple loop
    int buildMapStopFrame = m_points_features.size();
    int startDetectLoopIndex = std::max(buildMapStopFrame - 500, 0);
    int endDetectLoopIndex = buildMapStopFrame;
    int topK = 20;

    for (int idx = startDetectLoopIndex; idx < endDetectLoopIndex; idx = idx + 10) {
        DBoW3::QueryResults result;
        m_db.query(m_points_features[idx].features, result, topK);

        DBoW3::QueryResults::const_iterator qit;
        size_t rowIdx = 0;
        int imageIDs[topK] = {0};
        double scores[topK] = {0.0};
        for (qit = result.begin(); qit != result.end(); ++qit) {
            imageIDs[rowIdx] = (double)qit->Id;   // 注意C++索引从0开始
            scores[rowIdx] = (double)qit->Score;  // 注意前面程序用的是topK
            rowIdx++;
        }

        std::vector<int> loopKeyFrameIds, validKeyFrameIds;
        std::vector<double> keyFrameScores, validKeyFrameScores;
        int nearestIDsStartIdx = std::max(0, idx - 100);  // 建图模式，位于当前帧至少100帧前才开始回环检测
        for (int iq = 0; iq < topK; iq++) {
            if (imageIDs[iq] < nearestIDsStartIdx) {
                loopKeyFrameIds.push_back(imageIDs[iq]);
                keyFrameScores.push_back(scores[iq]);
            }
        }

        bool isDetectedLoop = false;
        int loopCandidate = 0;
        if (!loopKeyFrameIds.empty()) {
            double minScore = *std::min_element(scores, scores + topK);
            double bestScore = keyFrameScores[0];
            double ratio = 0.6;                                               // 根据实际图像调整此值
            double threshold = std::max({bestScore * ratio, minScore, 0.2});  // 0.2根据数据集经验取得, 见doc/loopClosureDetect.md
            std::vector<int> validIdxs = findItems(keyFrameScores, threshold);
            for (size_t validId = 0; validId < validIdxs.size(); validId++) {
                validKeyFrameIds.push_back(loopKeyFrameIds[validId]);
            }
            if (validKeyFrameIds.size() > 2) {
                std::vector<std::vector<int>> groups = nchoosek(validKeyFrameIds, 3);
                cv::Mat matGroups(groups.size(), groups.at(0).size(), CV_64FC1);
                for (int iRow = 0; iRow < matGroups.rows; ++iRow) {
                    for (int jCol = 0; jCol < matGroups.cols; ++jCol) {
                        matGroups.at<double>(iRow, jCol) = groups.at(iRow).at(jCol);
                    }
                }

                for (size_t i_r = 0; i_r < matGroups.rows; i_r++) {
                    double minV, maxV;
                    cv::minMaxLoc(matGroups.row(i_r), &minV, &maxV, NULL, NULL);
                    if (maxV - minV <= 3) {
                        cv::Mat consecutiveGroups = matGroups.row(i_r);
                        loopCandidate = (int)consecutiveGroups.at<double>(0, 0);
                        isDetectedLoop = true;
                        break;
                    } else {
                        isDetectedLoop = false;
                    }
                }
            }
        }

        if (isDetectedLoop) {
            std::vector<cv::Point2f> preP, nextP;
            std::vector<cv::DMatch> matches;
            cv::BFMatcher matcher(cv::NORM_HAMMING);
            std::vector<std::vector<cv::DMatch>> matchePoints;
            matcher.knnMatch(m_points_features[loopCandidate].features, m_points_features[idx].features, matchePoints, 2);
            for (int i = 0; i < matchePoints.size(); i++) {
                if (matchePoints[i][0].distance < 0.5 * matchePoints[i][1].distance) {
                    matches.push_back(matchePoints[i][0]);
                }
            }
            preP.resize(matches.size());
            nextP.resize(matches.size());
            for (size_t i = 0; i < matches.size(); i++) {
                preP[i] = m_points_features[loopCandidate].keyPoints[matches[i].queryIdx].pt;
                nextP[i] = m_points_features[idx].keyPoints[matches[i].trainIdx].pt;
            }
            cv::Mat loopTform, inliers;
            int status = 0;
            estiTform(preP, nextP, loopTform, inliers, status);
            if ((status > 0) || (cv::sum(inliers)[0] <= 3)) {
                continue;
            }

            cv::Mat invertR = loopTform.rowRange(0, 2).colRange(0, 2).t();
            cv::Mat invertLoopTform;
            cv::hconcat(invertR, -invertR * loopTform.col(2), invertLoopTform);
            cv::Mat tempR = (cv::Mat_<double>(3, 3) << invertLoopTform.at<double>(0, 0), invertLoopTform.at<double>(0, 1), 0,
                             invertLoopTform.at<double>(1, 0), invertLoopTform.at<double>(1, 1), 0,
                             0, 0, 1);
            std::vector<double> angRadius;
            cv::Rodrigues(tempR, angRadius);
            double measurement[3] =
                {invertLoopTform.at<double>(0, 2), invertLoopTform.at<double>(1, 2), angRadius[2]};  // 注意角度，是相对值, 单位为弧度

            // optimize loop
            cv::Mat loopIDpairs = (cv::Mat_<double>(1, 2) << idx, loopCandidate);
            std::vector<cv::Vec3d> relPoses{cv::Vec3d(measurement)};
            optimizePoseGraph(loopIDpairs, relPoses);
            break;  // 仅使用一次闭环即可
        }
    }
}

void buildMapping::HDMapping::fuseOptimizeHDMap(std::string imageFilesList, std::vector<cv::Vec3d>& updateNodeVehiclePtPoses) {
    int numNodes = updateNodeVehiclePtPoses.size();
    cv::Mat relTformA = (cv::Mat_<double>(3, 3) << 1, 0, -updateNodeVehiclePtPoses[0][0],
                         0, 1, -updateNodeVehiclePtPoses[0][1],
                         0, 0, 1);

    m_HDmapOutput.ref.XWorldLimits[0] = 0.0;
    m_HDmapOutput.ref.XWorldLimits[1] = 0.5;
    m_HDmapOutput.ref.YWorldLimits[0] = 0.0;
    m_HDmapOutput.ref.YWorldLimits[1] = 0.5;
    m_HDmapOutput.bigImg.release();

    int idx = 0;
    std::ifstream fid(imageFilesList);
    std::string line;
    while (std::getline(fid, line) && (idx < numNodes)) {
        cv::Mat currImg = cv::imread(line);
        double theta = updateNodeVehiclePtPoses[idx][2];  // in radian
        cv::Mat currNodeVehiclePtPoseA = (cv::Mat_<double>(3, 3) << std::cos(theta), -std::sin(theta), updateNodeVehiclePtPoses[idx][0],
                                          std::sin(theta), std::cos(theta), updateNodeVehiclePtPoses[idx][1],
                                          0, 0, 1);
        cv::Mat imageTargetPoseA = currNodeVehiclePtPoseA * relTformA;
        cv::Mat tform = imageTargetPoseA.rowRange(0, 2);
        cv::Mat currCorner = (cv::Mat_<double>(3, 4) << 0, currImg.cols - 1, currImg.cols - 1, 0,
                              0, 0, currImg.rows - 1, currImg.rows - 1,
                              1, 1, 1, 1);
        currCorner = tform * currCorner;

        double xmin, xmax, ymin, ymax;
        cv::minMaxLoc(currCorner.row(0), &xmin, &xmax);
        cv::minMaxLoc(currCorner.row(1), &ymin, &ymax);
        double pad_left = 0, pad_right = 0, pad_top = 0, pad_down = 0;
        if (xmin < m_HDmapOutput.ref.XWorldLimits[0]) {
            pad_left = m_HDmapOutput.ref.XWorldLimits[0] - xmin;
            m_HDmapOutput.ref.XWorldLimits[0] = xmin;
        }

        if (xmax > m_HDmapOutput.ref.XWorldLimits[1]) {
            pad_right = xmax - m_HDmapOutput.ref.XWorldLimits[1];
            m_HDmapOutput.ref.XWorldLimits[1] = xmax;
        }

        if (ymin < m_HDmapOutput.ref.YWorldLimits[0]) {
            pad_top = m_HDmapOutput.ref.YWorldLimits[0] - ymin;
            m_HDmapOutput.ref.YWorldLimits[0] = ymin;
        }

        if (ymax > m_HDmapOutput.ref.YWorldLimits[1]) {
            pad_down = ymax - m_HDmapOutput.ref.YWorldLimits[1];
            m_HDmapOutput.ref.YWorldLimits[1] = ymax;
        }

        m_HDmapOutput.ref.ImageSize[0] = std::round(m_HDmapOutput.ref.YWorldLimits[1] - m_HDmapOutput.ref.YWorldLimits[0]) + 1;
        m_HDmapOutput.ref.ImageSize[1] = std::round(m_HDmapOutput.ref.XWorldLimits[1] - m_HDmapOutput.ref.XWorldLimits[0]) + 1;

        // 平移到可视化区域
        cv::Mat topImg, bwMask;
        tform.colRange(2, 3) = tform.colRange(2, 3) - (cv::Mat_<double>(2, 1) << m_HDmapOutput.ref.XWorldLimits[0], m_HDmapOutput.ref.YWorldLimits[0]);
        cv::warpAffine(currImg, topImg, tform, cv::Size(m_HDmapOutput.ref.ImageSize[1], m_HDmapOutput.ref.ImageSize[0]));
        if (idx == 0) {
            cv::Mat temp = cv::Mat::ones(m_BW.rows, m_BW.cols, m_BW.type());
            cv::warpAffine(temp, bwMask, tform, cv::Size(m_HDmapOutput.ref.ImageSize[1], m_HDmapOutput.ref.ImageSize[0]));
            m_HDmapOutput.bigImg = currImg;
        } else {
            cv::warpAffine(m_BW, bwMask, tform, cv::Size(m_HDmapOutput.ref.ImageSize[1], m_HDmapOutput.ref.ImageSize[0]));
        }
        cv::copyMakeBorder(m_HDmapOutput.bigImg, m_HDmapOutput.bigImg, std::round(pad_top), std::round(pad_down),
                           std::round(pad_left), std::round(pad_right), cv::BORDER_CONSTANT, cv::Scalar::all(0));
        cv::resize(m_HDmapOutput.bigImg, m_HDmapOutput.bigImg, cv::Size(topImg.cols, topImg.rows));
        topImg.copyTo(m_HDmapOutput.bigImg, bwMask);

        idx = idx + 1;
    }
    fid.close();
}

void buildMapping::HDMapping::drawRoutePath(cv::Mat& drawImage) const {
    drawImage.release();
    if (m_HDmapOutput.bigImg.channels() == 3) {
        drawImage = m_HDmapOutput.bigImg;
    } else {
        cv::cvtColor(m_HDmapOutput.bigImg, drawImage, cv::COLOR_GRAY2BGR);
    }

    int lenPoses = m_vehiclePoses.size();
    if (lenPoses > 1) {
        for (size_t i = 0; i < lenPoses - 1; i++) {
            cv::Point pt1, pt2;
            cv::Vec3d prePose, nextPose;
            prePose = m_vehiclePoses[i];
            nextPose = m_vehiclePoses[i + 1];
            pt1 = cv::Point(prePose[0] - m_HDmapOutput.ref.XWorldLimits[0], prePose[1] - m_HDmapOutput.ref.YWorldLimits[0]);
            pt2 = cv::Point(nextPose[0] - m_HDmapOutput.ref.XWorldLimits[0], nextPose[1] - m_HDmapOutput.ref.YWorldLimits[0]);
            cv::line(drawImage, pt1, pt2, cv::Scalar(0, 0, 255), 4);
        }
    }
}

void buildMapping::HDMapping::optimizePoseGraph(cv::Mat& loopIDpairs, std::vector<cv::Vec3d>& relPoses) {
    coder::array<double, 2U> absposes_tmp, rel_poses;
    coder::array<double, 2U> loopNodePairs;
    coder::array<double, 2U> updatedPoses;

    absposes_tmp.set_size(m_vehiclePoses.size(), 3);
    for (size_t i = 0; i < m_vehiclePoses.size(); i++) {
        absposes_tmp[i] = m_vehiclePoses[i][0];
        absposes_tmp[i + m_vehiclePoses.size()] = m_vehiclePoses[i][1];
        absposes_tmp[i + 2 * m_vehiclePoses.size()] = m_vehiclePoses[i][2];
    }
    rel_poses.set_size(relPoses.size(), 3);
    for (size_t i = 0; i < relPoses.size(); i++) {
        rel_poses[i] = relPoses[i][0];
        rel_poses[i + relPoses.size()] = relPoses[i][1];
        rel_poses[i + 2 * relPoses.size()] = relPoses[i][2];
    }
    loopNodePairs.set_size(loopIDpairs.rows, 2);
    for (size_t i = 0; i < loopIDpairs.rows; i++) {
        double* numel = loopIDpairs.ptr<double>(i);
        loopNodePairs[i] = numel[0];
        loopNodePairs[i + loopIDpairs.rows] = numel[1];
    }
    poseGraphOptimize::poseGraphOptimize(absposes_tmp, loopNodePairs,
                                         rel_poses, updatedPoses);
    m_vehiclePoses.clear();
    for (size_t i = 0; i < updatedPoses.size(0); i++) {
        m_vehiclePoses.push_back(cv::Vec3d(updatedPoses[i], updatedPoses[i + updatedPoses.size(0)], updatedPoses[i + 2 * updatedPoses.size(0)]) +
                                 m_initViclePtPose);
    }
    poseGraphOptimize::poseGraphOptimize_terminate();
}

void buildMapping::HDMapping::detectLanes(cv::Mat& grayImage, std::vector<cv::line_descriptor::KeyLine>& lines) {
    if (grayImage.channels() == 3) {
        cv::cvtColor(grayImage, grayImage, cv::COLOR_BGR2GRAY);
    }

    coder::array<detectLaneMarkerRidge::struct0_T, 2U> mat_lines_struct;
    coder::array<unsigned char, 2U> bevImage_tmp, mask_tmp;
    double approxLaneWidthPixels_tmp = 4;
    double sensitive = 0.25;

    bevImage_tmp.set_size(grayImage.rows, grayImage.cols);
    convertCVToMatrix(grayImage, grayImage.rows, grayImage.cols, grayImage.channels(), bevImage_tmp.data());

    cv::Mat mask = cv::Mat::ones(grayImage.rows, grayImage.cols, CV_8UC1);
    mask_tmp.set_size(grayImage.rows, grayImage.cols);
    convertCVToMatrix(mask, mask.rows, mask.cols, mask.channels(), mask_tmp.data());

    detectLaneMarkerRidge::detectLaneMarkerRidge(
        bevImage_tmp, approxLaneWidthPixels_tmp, mask_tmp,
        sensitive, mat_lines_struct);

    lines.clear();
    for (size_t i = 0; i < mat_lines_struct.size(1); i++) {
        detectLaneMarkerRidge::struct0_T currLine = mat_lines_struct[i];
        cv::line_descriptor::KeyLine kl;
        kl.startPointX = currLine.point1[0];
        kl.startPointY = currLine.point1[1];
        kl.endPointX = currLine.point2[0];
        kl.endPointY = currLine.point2[1];
        kl.sPointInOctaveX = currLine.point1[0];
        kl.sPointInOctaveY = currLine.point1[1];
        kl.ePointInOctaveX = currLine.point2[0];
        kl.ePointInOctaveY = currLine.point2[1];
        kl.lineLength = (float)sqrt(pow(currLine.point1[0] - currLine.point2[0], 2) +
                                    pow(currLine.point1[1] - currLine.point2[1], 2));
        cv::LineIterator li(grayImage, cv::Point2f(currLine.point1[0], currLine.point1[1]),
                            cv::Point2f(currLine.point2[0], currLine.point2[1]));
        kl.numOfPixels = li.count;

        kl.angle = atan2((kl.endPointY - kl.startPointY), (kl.endPointX - kl.startPointX));
        kl.class_id = i;  // or zero ?
        kl.octave = 0;
        kl.size = (kl.endPointX - kl.startPointX) * (kl.endPointY - kl.startPointY);
        kl.response = kl.lineLength / std::max(grayImage.cols, grayImage.rows);
        kl.pt = cv::Point2f((kl.endPointX + kl.startPointX) / 2, (kl.endPointY + kl.startPointY) / 2);
        lines.push_back(kl);
    }
    detectLaneMarkerRidge::detectLaneMarkerRidge_terminate();
}