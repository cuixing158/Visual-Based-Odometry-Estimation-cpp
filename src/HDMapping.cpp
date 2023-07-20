#include <fstream>
#include "HDMapping.h"

#include "selectUniform2.h"
#include "coder_array.h"
#include "estimateAffineRigid2D.h"
#include "rt_nonfinite.h"

#define DEBUG_SHOW_IMAGE 1
namespace buildMapping {
class HDMapping;
}

void ssc(std::vector<cv::KeyPoint> keyPoints, int numRetPoints,
         float tolerance, int cols, int rows, std::vector<cv::KeyPoint>& outputPts, std::vector<int>& indexs) {
    // several temp expression variables to simplify solution equation
    int exp1 = rows + cols + 2 * numRetPoints;
    long long exp2 =
        ((long long)4 * cols + (long long)4 * numRetPoints +
         (long long)4 * rows * numRetPoints + (long long)rows * rows +
         (long long)cols * cols - (long long)2 * rows * cols +
         (long long)4 * rows * cols * numRetPoints);
    double exp3 = sqrt(exp2);
    double exp4 = numRetPoints - 1;

    double sol1 = -round((exp1 + exp3) / exp4);  // first solution
    double sol2 = -round((exp1 - exp3) / exp4);  // second solution

    // binary search range initialization with positive solution
    int high = (sol1 > sol2) ? sol1 : sol2;
    int low = floor(sqrt((double)keyPoints.size() / numRetPoints));
    low = std::max(1, low);

    int width;
    int prevWidth = -1;

    std::vector<int> ResultVec;
    bool complete = false;
    unsigned int K = numRetPoints;
    unsigned int Kmin = round(K - (K * tolerance));
    unsigned int Kmax = round(K + (K * tolerance));

    std::vector<int> result;
    result.reserve(keyPoints.size());
    while (!complete) {
        width = low + (high - low) / 2;
        if (width == prevWidth ||
            low >
                high) {          // needed to reassure the same radius is not repeated again
            ResultVec = result;  // return the keypoints from the previous iteration
            break;
        }
        result.clear();
        double c = (double)width / 2.0;  // initializing Grid
        int numCellCols = floor(cols / c);
        int numCellRows = floor(rows / c);
        std::vector<std::vector<bool>> coveredVec(numCellRows + 1,
                                                  std::vector<bool>(numCellCols + 1, false));

        for (unsigned int i = 0; i < keyPoints.size(); ++i) {
            int row =
                floor(keyPoints[i].pt.y /
                      c);  // get position of the cell current point is located at
            int col = floor(keyPoints[i].pt.x / c);
            if (coveredVec[row][col] == false) {  // if the cell is not covered
                result.push_back(i);
                int rowMin = ((row - floor(width / c)) >= 0)
                                 ? (row - floor(width / c))
                                 : 0;  // get range which current radius is covering
                int rowMax = ((row + floor(width / c)) <= numCellRows)
                                 ? (row + floor(width / c))
                                 : numCellRows;
                int colMin =
                    ((col - floor(width / c)) >= 0) ? (col - floor(width / c)) : 0;
                int colMax = ((col + floor(width / c)) <= numCellCols)
                                 ? (col + floor(width / c))
                                 : numCellCols;
                for (int rowToCov = rowMin; rowToCov <= rowMax; ++rowToCov) {
                    for (int colToCov = colMin; colToCov <= colMax; ++colToCov) {
                        if (!coveredVec[rowToCov][colToCov])
                            coveredVec[rowToCov][colToCov] =
                                true;  // cover cells within the square bounding box with width
                                       // w
                    }
                }
            }
        }

        if (result.size() >= Kmin && result.size() <= Kmax) {  // solution found
            ResultVec = result;
            complete = true;
        } else if (result.size() < Kmin)
            high = width - 1;  // update binary search range
        else
            low = width + 1;
        prevWidth = width;
    }
    // retrieve final keypoints
    outputPts.clear();
    indexs.clear();
    for (unsigned int i = 0; i < ResultVec.size(); i++)
        outputPts.push_back(keyPoints[ResultVec[i]]);

    indexs = ResultVec;
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
    initViclePtPose = cv::Vec3f((285 + 358) / 2.0, (156, 326) / 2.0, 0);

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
    orbDetector->detect(currImg, keyPts, orbDetectMask);
    std::vector<int> indexs;
#ifdef USE_SSC_UNIFORM
    // selectUniform(keyPts, Descriptions, 2000, currKeypts, currDescriptors);
    currKeypts.clear();
    ssc(keyPts, 2000, 0.1, currImg.cols, currImg.rows, currKeypts, indexs);

#else
    selectUniformPoints(keyPts, 2000, cv::Size(currImg.cols, currImg.rows), currKeypts, indexs);
#endif
    orbDetector->compute(currImg, currKeypts, currDescriptors);
    if (preKeypts.empty()) {
        prevImg = currImg;
        preKeypts = currKeypts;
        preDescriptors = currDescriptors;

        HDmapOutput.bigImg = currImg;
        HDmapOutput.ref.XWorldLimits[0] = 0;
        HDmapOutput.ref.XWorldLimits[1] = currImg.cols - 1;
        HDmapOutput.ref.YWorldLimits[0] = 0;
        HDmapOutput.ref.YWorldLimits[1] = currImg.rows - 1;
    }
    //match features
    std::vector<cv::DMatch> matches;
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
    std::vector<std::vector<cv::DMatch>> matchePoints;
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

//     cv::Mat tempShowImg;
// // cv::drawKeypoints(currImg, currKeypts, tempShowImg);
// tempShowImg = currImg;
// cv::RNG rng(time(0));
// for (size_t i = 0; i < points2.size(); i++) {
//     cv::circle(tempShowImg, points2[i], 2, cv::Scalar::all(255));
// }
// cv::imwrite("tempShow.jpg", tempShowImg);
#else
    std::vector<cv::Point> preP, nextP;
    preP.resize(matches.size());
    nextP.resize(matches.size());
    for (size_t i = 0; i < matches.size(); i++) {
        preP[i] = preKeypts[matches[i].queryIdx].pt;
        nextP[i] = currKeypts[matches[i].trainIdx].pt;
    }
    cv::Mat inliers = cv::Mat::zeros(preP.size(), 1, CV_8UC1);
    int status = 0;
    estiTform(preP, nextP, relTform, inliers, status);
#endif
    bool cond = (status > 0 || cv::sum(inliers)[0] <= 3);
    // std::cout << "inliers:" << inliers << std::endl;
    // cv::Scalar ss = cv::sum(inliers);
    // std::cout << __LINE__ << ",cv::sum(inliers)[0]:" << ss << "," << cv::sum(inliers)[0] << std::endl;
    if (cond) {
        std::vector<uchar> statusOK;
        std::vector<float> err;
        cv::TermCriteria criteria = cv::TermCriteria((cv::TermCriteria::COUNT) + (cv::TermCriteria::EPS), 10, 0.03);
        std::vector<cv::Point2f> p0, p1;
        cv::KeyPoint::convert(preKeypts, p0);
        cv::calcOpticalFlowPyrLK(prevImg, currImg, p0, p1, statusOK, err, cv::Size(15, 15), 2, criteria);
        std::vector<cv::Point2f> good_new;
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
        if (cond) {  // may be lost
            relTform = preRelTform;

            std::cout << "may be lost" << std::endl;
            // // debug output
            // cv::Mat imgMatches;
            // std::vector<cv::KeyPoint> temp1, temp2;
            // matches.clear();  // debug
            // for (size_t i = 0; i < preP.size(); i++) {
            //     temp1.push_back(cv::KeyPoint(preP[i], 1.f));
            //     temp2.push_back(cv::KeyPoint(nextP[i], 1.f));
            //     matches.push_back(cv::DMatch(i, i, 1));
            // }
            // cv::drawMatches(prevImg, temp1, currImg, temp2, matches, imgMatches, cv::Scalar(0, 255, 255), -1, std::vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
            // cv::imwrite("matchesForce" + std::to_string(num) + ".jpg", imgMatches);
        }
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
    std::cout << "matches ratio:" << matches.size() << "/" << matchePoints.size() << std::endl;
    std::cout << "rigid estimate matches ratio:" << rigidMatches.size() << "/" << matches.size() << std::endl;

    assert(rigidMatches.size() > 3);

    // Draw top matches
    cv::Mat imMatches;
    cv::drawMatches(prevImg, preKeypts, currImg, currKeypts, rigidMatches, imMatches, cv::Scalar(0, 255, 255), -1, std::vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
    cv::imwrite("matches.jpg", imMatches);
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
    cv::imwrite("bigImgCopy.jpg", drawImage);
#endif

    // update previous state variables
    prevImg = currImg;
    preKeypts = currKeypts;
    preDescriptors = currDescriptors;
    previousImgPose = currImgPose;
    preRelTform = relTform;
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

void buildMapping::HDMapping::estiTform(std::vector<cv::Point>& prePoints, std::vector<cv::Point>& currPoints, cv::Mat& tform2x3, cv::Mat& inliers, int& status) {
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