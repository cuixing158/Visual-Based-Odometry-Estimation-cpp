#include "HDMapping.h"

namespace buildMapping {
class HDMapping;
}

buildMapping::HDMapping::HDMapping() {
    double cumDist = 0.0;
    double pixelExtentInWorldXY = 0.0;
    bool isBuildMap = true;
    double buildMapStopFrame = 1200;
    bool isBuildMapOver = false;
    bool isLocSuccess = false;
    double locVehiclePose[3] = {0, 0, 0};

    orbDetectMask = 255 * cv::Mat::ones(480, 640, CV_8UC1);
    orbDetectMask.rowRange(156, 326).colRange(285, 358) = 0;
    BW = 255 * cv::Mat::ones(480, 640, CV_8UC1);
    BW.rowRange(480 / 2, 480) = 0;
    BW = BW &= orbDetectMask;

    HDmapOutput.ref.XWorldLimits[0] = 0;
    HDmapOutput.ref.XWorldLimits[1] = 0.5;
    HDmapOutput.ref.YWorldLimits[0] = 0;
    HDmapOutput.ref.YWorldLimits[1] = 0.5;

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

void buildMapping::HDMapping::constructWorldMap(const cv::Mat& srcImage) {
    static size_t num = 0;
    num++;

    cv::Mat currImg = srcImage;
    if (currImg.channels() == 3)
        cv::cvtColor(currImg, currImg, cv::COLOR_BGR2GRAY);

    if (currImg.empty()) throw std::runtime_error("currImg is empty!");

    std::vector<cv::KeyPoint> keyPts;
    cv::Mat Descriptions;
    orbDetector->detectAndCompute(currImg, orbDetectMask, keyPts, Descriptions);
    selectUniform(keyPts, Descriptions, 2000, currKeypts, currDescriptors);
    if (preKeypts.empty()) {
        prevImg = currImg;
        preKeypts = currKeypts;
        preDescriptors = currDescriptors;
    }
    //match features
    std::vector<cv::DMatch> matches, good_matches;
    cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create("FlannBased");
    matcher->match(preDescriptors, currDescriptors, matches);
    std::vector<cv::Point2f> points1, points2;
    for (size_t i = 0; i < matches.size(); i++) {
        points1.push_back(preKeypts[matches[i].queryIdx].pt);
        points2.push_back(currKeypts[matches[i].trainIdx].pt);
    }

    // estimate geometry rigid 2D transformation matrix
    cv::Mat inliers = cv::Mat::zeros(preKeypts.size(), 1, CV_8U);
    cv::Mat optimalAffineMat = estimateAffinePartial2D(points1, points2, inliers, cv::RANSAC);
    cv::Mat R = optimalAffineMat.rowRange(0, 2).colRange(0, 2);
    double s = std::sqrt(cv::determinant(R));
    cv::Mat rigidtform2dR = R.mul(1.0 / s);
    cv::hconcat(rigidtform2dR, optimalAffineMat.col(2), relTform);

    // build map mode
    cv::Mat previousImgPoseA = previousImgPose;
    cv::Mat homoMatrix = (cv::Mat_<double>(1, 3) << 0, 0, 1);
    previousImgPoseA.push_back(homoMatrix);
    cv::Mat currImgPose = relTform * previousImgPoseA;

    cv::Mat tform, tempR, tempT;
    tempR = currImgPose.rowRange(0, 2).colRange(0, 2).t();
    tempT = -tempR * currImgPose.col(2);
    cv::hconcat(tempR, tempT, tform);

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
    // if (num > 52) {
    //     std::cout << std::round(pad_top) << "," << std::round(pad_down) << "," << std::round(pad_left) << "," << std::round(pad_right) << std::endl;
    //     std::cout << "HDmapOutput.bigImg.size:" << HDmapOutput.bigImg.size << ",topImg.size:" << topImg.size << ",BW.size:" << BW.size << std::endl;
    //     cv::imwrite("bigImg" + std::to_string(num) + ".jpg", HDmapOutput.bigImg);
    //     cv::imwrite("BW" + std::to_string(num) + ".jpg", BW);
    //     cv::imwrite("topImg" + std::to_string(num) + ".jpg", topImg);
    // }

    topImg.copyTo(HDmapOutput.bigImg, bwMask);
    for (size_t i = 0; i < matches.size(); i++) {
        if (inliers.at<uchar>(i, 0)) {
            good_matches.push_back(matches[i]);
        }
    }

#if DEBUG_SHOW_IMAGE
    // Draw top matches
    cv::Mat imMatches;
    cv::drawMatches(prevImg, preKeypts, currImg, currKeypts, good_matches, imMatches, cv::Scalar(0, 255, 255), -1, std::vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
    cv::imwrite("matches.jpg", imMatches);
    // cv::Mat dstMat;
    // srcImage.copyTo(dstMat, orbDetectMask);
    cv::imwrite("bigImgCopy" + std::to_string(num) + ".jpg", HDmapOutput.bigImg);
    //
#endif

    // update previous state variables
    prevImg = currImg;
    preKeypts = currKeypts;
    preDescriptors = currDescriptors;
    previousImgPose = currImgPose;
}

void buildMapping::HDMapping::selectUniform(std::vector<cv::KeyPoint>& keyPts, cv::Mat& Descriptions, size_t numPoints, std::vector<cv::KeyPoint>& outKeypts, cv::Mat& outDescriptions) {
    outKeypts.clear();
    outDescriptions = cv::Mat();

    std::vector<size_t> uniformIdxs(keyPts.size());
    std::iota(uniformIdxs.begin(), uniformIdxs.end(), 0);
    std::random_shuffle(uniformIdxs.begin(), uniformIdxs.end());

    if (uniformIdxs.size() > numPoints) {
        for (size_t i = 0; i < numPoints; i++) {
            size_t idx = uniformIdxs[i];
            outKeypts.push_back(keyPts[idx]);
            outDescriptions.push_back(Descriptions.row(idx));
        }
    } else {
        outKeypts = keyPts;
        outDescriptions = Descriptions;
    }
}