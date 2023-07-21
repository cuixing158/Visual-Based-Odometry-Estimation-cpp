#include <fstream>
#include "HDMapping.h"

#include "selectUniform2.h"
#include "coder_array.h"
#include "estimateAffineRigid2D.h"
#include "rt_nonfinite.h"

#define DEBUG_SHOW_IMAGE 1
// #define USE_SSC_UNIFORM
namespace buildMapping {
class HDMapping;
}

buildMapping::HDMapping::HDMapping() {
    cumDist = 0.0;
    pixelExtentInWorldXY = 0.015 / 0.51227;
    isBuildMap = true;
    buildMapStopFrame = 1200;
    isBuildMapOver = false;
    isLocSuccess = false;
    locVehiclePose[3] = {0};

    orbDetectMask = 255 * cv::Mat::ones(480, 640, CV_8UC1);
    orbDetectMask.rowRange(156, 326).colRange(285, 358) = 0;
    BW = 255 * cv::Mat::ones(480, 640, CV_8UC1);
    BW.rowRange(480 / 2, 480) = 0;
    BW = BW &= orbDetectMask;
    initViclePtPose = cv::Vec3f((285 + 358) / 2.0, (156, 326) / 2.0, 0);
    method = matchFeatureMethod::LK_TRACK_FEATURES;

    //第一副图像的像素坐标系为世界坐标系
    preRelTform = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                   0, 1, 0);
    relTform = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                0, 1, 0);
    previousImgPose = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                       0, 1, 0);
    orbDetector = cv::ORB::create(307200);
};
buildMapping::HDMapping::~HDMapping() {}

void buildMapping::HDMapping::reset() {
    prevImg = cv::Mat();
    preKeypts.clear();
    preDescriptors = cv::Mat();

    cumDist = 0.0;
    vehiclePoses.clear();

    HDmapOutput.bigImg = cv::Mat();
    HDmapOutput.ref.XWorldLimits[0] = 0;
    HDmapOutput.ref.XWorldLimits[1] = 1;
    HDmapOutput.ref.YWorldLimits[0] = 0;
    HDmapOutput.ref.YWorldLimits[1] = 1;

    preRelTform = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                   0, 1, 0);
    relTform = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                0, 1, 0);
    previousImgPose = (cv::Mat_<double>(2, 3) << 1, 0, 0,
                       0, 1, 0);
    orbDetector = cv::ORB::create(307200);
}

void buildMapping::HDMapping::saveMapData() {
}

void buildMapping::HDMapping::loadMapData() {
}

void buildMapping::HDMapping::constructWorldMap(const cv::Mat& srcImage) {
    static size_t num = 0;
    num++;

    cv::Mat currImg = srcImage;
    if (currImg.channels() == 3)
        cv::cvtColor(currImg, currImg, cv::COLOR_BGR2GRAY);

    if (currImg.empty()) throw std::runtime_error("currImg is empty!");

    std::vector<cv::KeyPoint> keyPts;
    orbDetector->detect(currImg, keyPts, orbDetectMask);
    std::vector<int> indexs;
    selectUniformPoints(keyPts, 2000, cv::Size(currImg.cols, currImg.rows), currKeypts, indexs);
    if (preKeypts.empty()) {
        prevImg = currImg;
        preKeypts = currKeypts;

        HDmapOutput.bigImg = currImg;
        HDmapOutput.ref.XWorldLimits[0] = 0;
        HDmapOutput.ref.XWorldLimits[1] = currImg.cols - 1;
        HDmapOutput.ref.YWorldLimits[0] = 0;
        HDmapOutput.ref.YWorldLimits[1] = currImg.rows - 1;
    }

    std::vector<cv::DMatch> matches;
    std::vector<cv::Point2f> preP, nextP;
    cv::Mat inliers;
    int status = 0;
    int rigidTformType = 0;
    if (method == buildMapping::HDMapping::matchFeatureMethod::HYBRID_FEATURES) {
        orbDetector->compute(currImg, currKeypts, currDescriptors);
        if (preDescriptors.empty()) {
            preDescriptors = currDescriptors;
        }
        //match features
        cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create("FlannBased");
        if (
            preDescriptors.type() != CV_32F) {
            preDescriptors.convertTo(preDescriptors, CV_32F);
        }
        if (currDescriptors.type() != CV_32F) {
            currDescriptors.convertTo(currDescriptors, CV_32F);
        }

#ifdef USE_TOPK_BEST_MATCHES
        matcher->match(preDescriptors, currDescriptors, matches);
        //只初步选取匹配前num个较好的特征点
        int numMatches = matches.size() < 10 ? matches.size() : matches.size() / 2;
        std::nth_element(matches.begin(), matches.begin() + numMatches, matches.end());
        matches.erase(matches.begin() + numMatches, matches.end());
#else
        std::vector<std::vector<cv::DMatch> > matchePoints;
        matcher->knnMatch(preDescriptors, currDescriptors, matchePoints, 2);
        for (int i = 0; i < matchePoints.size(); i++) {
            if (matchePoints[i][0].distance < 0.6 * matchePoints[i][1].distance) {
                matches.push_back(matchePoints[i][0]);
            }
        }
#endif

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
        cv::hconcat(rigidtform2dR, optimalAffineMat.col(2), relTform);
#else
        preP.resize(matches.size());
        nextP.resize(matches.size());
        for (size_t i = 0; i < matches.size(); i++) {
            preP[i] = preKeypts[matches[i].queryIdx].pt;
            nextP[i] = currKeypts[matches[i].trainIdx].pt;
        }
        estiTform(preP, nextP, relTform, inliers, status);
#endif
        bool cond = (status > 0 || cv::sum(inliers)[0] <= 3);
        if (cond) {
            rigidTformType = 1;
            std::vector<uchar> statusOK;
            std::vector<float> err;
            cv::TermCriteria criteria = cv::TermCriteria((cv::TermCriteria::COUNT) + (cv::TermCriteria::EPS), 10, 0.03);
            std::vector<cv::Point2f> p0, p1;
            cv::KeyPoint::convert(preKeypts, p0);
            cv::calcOpticalFlowPyrLK(prevImg, currImg, p0, p1, statusOK, err, cv::Size(15, 15), 2, criteria);
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
            estiTform(preP, nextP, relTform, inliers, status);
            cond = (status > 0 || cv::sum(inliers)[0] <= 3);
            if (cond) {
                rigidTformType = 2;
                relTform = preRelTform;
            }
        }
    } else if (method == buildMapping::HDMapping::matchFeatureMethod::LK_TRACK_FEATURES) {
        rigidTformType = 1;
        std::vector<uchar> statusOK;
        std::vector<float> err;
        cv::TermCriteria criteria = cv::TermCriteria((cv::TermCriteria::COUNT) + (cv::TermCriteria::EPS), 10, 0.03);
        std::vector<cv::Point2f> p0, p1;
        cv::KeyPoint::convert(preKeypts, p0);
        cv::calcOpticalFlowPyrLK(prevImg, currImg, p0, p1, statusOK, err, cv::Size(15, 15), 2, criteria);
        preP.clear();
        nextP.clear();
#if DEBUG_SHOW_IMAGE
        matches.clear();  // debug
#endif
        for (uint i = 0; i < p0.size(); i++) {
            if (statusOK[i] == 1) {
#if DEBUG_SHOW_IMAGE
                matches.push_back(cv::DMatch(preP.size(), nextP.size(), 1));
#endif
                preP.push_back(p0[i]);
                nextP.push_back(p1[i]);
            }
        }
        estiTform(preP, nextP, relTform, inliers, status);
        bool cond = (status > 0 || cv::sum(inliers)[0] <= 3);
        if (cond) {
            rigidTformType = 2;
            relTform = preRelTform;
#if DEBUG_SHOW_IMAGE
            std::ofstream fid("out_may_be_lost_frame.txt", std::ios_base::app);
            fid << num << std::endl;
            fid.close();
#endif
            saveMapData();
            reset();
        }
    } else {
        orbDetector->compute(currImg, currKeypts, currDescriptors);
        if (preDescriptors.empty()) {
            preDescriptors = currDescriptors;
        }
        //match features
        cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create("FlannBased");
        if (
            preDescriptors.type() != CV_32F) {
            preDescriptors.convertTo(preDescriptors, CV_32F);
        }
        if (currDescriptors.type() != CV_32F) {
            currDescriptors.convertTo(currDescriptors, CV_32F);
        }

#ifdef USE_TOPK_BEST_MATCHES
        matcher->match(preDescriptors, currDescriptors, matches);
        //只初步选取匹配前num个较好的特征点
        int numMatches = matches.size() < 10 ? matches.size() : matches.size() / 2;
        std::nth_element(matches.begin(), matches.begin() + numMatches, matches.end());
        matches.erase(matches.begin() + numMatches, matches.end());
#else
        std::vector<std::vector<cv::DMatch> > matchePoints;
        matcher->knnMatch(preDescriptors, currDescriptors, matchePoints, 2);
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
            preP[i] = preKeypts[matches[i].queryIdx].pt;
            nextP[i] = currKeypts[matches[i].trainIdx].pt;
        }
        estiTform(preP, nextP, relTform, inliers, status);
    }

    // build map mode
    cv::Mat previousImgPoseA = previousImgPose;
    cv::Mat homoMatrix = (cv::Mat_<double>(1, 3) << 0, 0, 1);
    previousImgPoseA.push_back(homoMatrix);
    cv::Mat currImgPose = relTform * previousImgPoseA;

    cv::Mat tform,
        tempR, tempT;
    tempR = currImgPose.rowRange(0, 2).colRange(0, 2).t();
    tempT = -tempR * currImgPose.col(2);
    cv::hconcat(tempR, tempT, tform);

    // calcuate vehicle poses
    cv::Mat initViclePtPoseA = (cv::Mat_<double>(3, 3) << 1, 0, initViclePtPose[0],
                                0, 1, initViclePtPose[1],
                                0, 0, 1);
    cv::Mat currVehiclePose = tform * initViclePtPoseA;
    tempR = (cv::Mat_<double>(3, 3) << currVehiclePose.at<double>(0, 0), currVehiclePose.at<double>(0, 1), 0,
             currVehiclePose.at<double>(1, 0), currVehiclePose.at<double>(1, 1), 0,
             0, 0, 1);
    std::vector<double> angRadius;
    cv::Rodrigues(tempR, angRadius);
    vehiclePoses.push_back(cv::Vec3d(currVehiclePose.at<double>(0, 2), currVehiclePose.at<double>(1, 2), angRadius[2]));
    if (vehiclePoses.size() >= 2) {
        cv::Vec3d prePose = vehiclePoses[vehiclePoses.size() - 2];
        cv::Vec3d currPose = vehiclePoses[vehiclePoses.size() - 1];
        cumDist += pixelExtentInWorldXY * std::sqrt((currPose[0] - prePose[0]) * (currPose[0] - prePose[0]) + (currPose[1] - prePose[1]) * (currPose[1] - prePose[1]));
    }

#if DEBUG_SHOW_IMAGE
    cv::Mat currCorner = (cv::Mat_<double>(3, 4) << 0, currImg.cols - 1, currImg.cols - 1, 0,
                          0, 0, currImg.rows - 1, currImg.rows - 1,
                          1, 1, 1, 1);
    currCorner = tform * currCorner;

    double xmin, xmax, ymin, ymax;
    cv::minMaxLoc(currCorner.row(0), &xmin, &xmax);
    cv::minMaxLoc(currCorner.row(1), &ymin, &ymax);
    double pad_left = 0, pad_right = 0, pad_top = 0, pad_down = 0;
    if (xmin < HDmapOutput.ref.XWorldLimits[0]) {
        pad_left = HDmapOutput.ref.XWorldLimits[0] - xmin;
        HDmapOutput.ref.XWorldLimits[0] = xmin;
    }

    if (xmax > HDmapOutput.ref.XWorldLimits[1]) {
        pad_right = xmax - HDmapOutput.ref.XWorldLimits[1];
        HDmapOutput.ref.XWorldLimits[1] = xmax;
    }

    if (ymin < HDmapOutput.ref.YWorldLimits[0]) {
        pad_top = HDmapOutput.ref.YWorldLimits[0] - ymin;
        HDmapOutput.ref.YWorldLimits[0] = ymin;
    }

    if (ymax > HDmapOutput.ref.YWorldLimits[1]) {
        pad_down = ymax - HDmapOutput.ref.YWorldLimits[1];
        HDmapOutput.ref.YWorldLimits[1] = ymax;
    }

    HDmapOutput.ref.ImageSize[0] = std::round(HDmapOutput.ref.YWorldLimits[1] - HDmapOutput.ref.YWorldLimits[0]) + 1;
    HDmapOutput.ref.ImageSize[1] = std::round(HDmapOutput.ref.XWorldLimits[1] - HDmapOutput.ref.XWorldLimits[0]) + 1;

    // 平移到可视化区域
    cv::Mat topImg, bwMask;
    tform.colRange(2, 3) = tform.colRange(2, 3) - (cv::Mat_<double>(2, 1) << HDmapOutput.ref.XWorldLimits[0], HDmapOutput.ref.YWorldLimits[0]);
    cv::warpAffine(currImg, topImg, tform, cv::Size(HDmapOutput.ref.ImageSize[1], HDmapOutput.ref.ImageSize[0]));
    cv::warpAffine(BW, bwMask, tform, cv::Size(HDmapOutput.ref.ImageSize[1], HDmapOutput.ref.ImageSize[0]));
    cv::copyMakeBorder(HDmapOutput.bigImg, HDmapOutput.bigImg, std::round(pad_top), std::round(pad_down), std::round(pad_left), std::round(pad_right), cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::resize(HDmapOutput.bigImg, HDmapOutput.bigImg, cv::Size(topImg.cols, topImg.rows));
    topImg.copyTo(HDmapOutput.bigImg, bwMask);
    std::vector<cv::DMatch> rigidMatches;
    for (size_t i = 0; i < matches.size(); i++) {
        if (inliers.at<uchar>(i, 0)) {
            rigidMatches.push_back(matches[i]);
        }
    }
    std::cout << "matches ratio:" << matches.size() << "/" << preKeypts.size() << std::endl;
    std::cout << "rigid estimate matches ratio:" << rigidMatches.size() << "/" << matches.size() << std::endl;

    // assert(rigidMatches.size() > 3);

    // Draw top matches
    cv::Mat imMatches;
    if (rigidTformType == 0) {
        cv::drawMatches(prevImg, preKeypts, currImg, currKeypts, rigidMatches, imMatches, cv::Scalar(0, 255, 255), -1, std::vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
    } else {
        std::vector<cv::KeyPoint> p1_temp, p2_temp;
        cv::KeyPoint::convert(preP, p1_temp, 1, 1, 0, -1);
        cv::KeyPoint::convert(nextP, p2_temp, 1, 1, 0, -1);
        cv::drawMatches(prevImg, p1_temp, currImg, p2_temp, rigidMatches, imMatches, cv::Scalar(0, 0, 255), -1, std::vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
    }
    cv::imwrite("matches_" + std::to_string(rigidTformType) + ".jpg", imMatches);
    // cv::Mat dstMat;
    // srcImage.copyTo(dstMat, orbDetectMask);
    cv::Mat drawImage;
    cv::cvtColor(HDmapOutput.bigImg, drawImage, cv::COLOR_GRAY2BGR);
    int lenPoses = vehiclePoses.size();
    if (lenPoses > 1) {
        for (size_t i = 0; i < lenPoses - 1; i++) {
            cv::Point pt1, pt2;
            cv::Vec3d prePose, nextPose;
            prePose = vehiclePoses[i];
            nextPose = vehiclePoses[i + 1];
            pt1 = cv::Point(prePose[0] - HDmapOutput.ref.XWorldLimits[0], prePose[1] - HDmapOutput.ref.YWorldLimits[0]);
            pt2 = cv::Point(nextPose[0] - HDmapOutput.ref.XWorldLimits[0], nextPose[1] - HDmapOutput.ref.YWorldLimits[0]);
            cv::line(drawImage, pt1, pt2, cv::Scalar(0, 0, 255), 4);
        }
    }

    cv::imwrite("bigImgCopy.png", drawImage);
#endif

    // update previous state variables
    prevImg = currImg;
    preKeypts = currKeypts;
    preDescriptors = currDescriptors;
    previousImgPose = currImgPose;
    preRelTform = relTform;
}

void buildMapping::HDMapping::estiTform(std::vector<cv::Point2f>& prePoints, std::vector<cv::Point2f>& currPoints, cv::Mat& tform2x3, cv::Mat& inliers, int& status) {
    coder::array<double, 2U>
        pts1_tmp, pts2_tmp;
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
    estimateAffineRigid2D::estimateAffineRigid2D(pts1_tmp, pts2_tmp, tform2x3_,
                                                 inlierIndex, &status);
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