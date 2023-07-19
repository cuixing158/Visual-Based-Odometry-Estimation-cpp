#include "HDMapping.h"

#include "selectUniform2.h"
#include "coder_array.h"
#include "estimateAffineRigid2D.h"
#include "rt_nonfinite.h"
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

void selectUniformPoints(std::vector<cv::KeyPoint> keyPoints, int numRetPoints,
                         cv::Size size, std::vector<cv::KeyPoint>& outputPts, std::vector<int>& indexs) {
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
    outputPts.clear();
    indexs.clear();
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
    std::vector<cv::DMatch> matches, good_matches;
    cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create("FlannBased");
    if (
        preDescriptors.type() != CV_32F) {
        preDescriptors.convertTo(preDescriptors, CV_32F);
    }
    if (currDescriptors.type() != CV_32F) {
        currDescriptors.convertTo(currDescriptors, CV_32F);
    }

    matcher->match(preDescriptors, currDescriptors, matches);
    std::vector<cv::Point2f> points1, points2;
    for (size_t i = 0; i < matches.size(); i++) {
        points1.push_back(preKeypts[matches[i].queryIdx].pt);
        points2.push_back(currKeypts[matches[i].trainIdx].pt);
    }
// cv::Mat tempShowImg;
// // cv::drawKeypoints(currImg, currKeypts, tempShowImg);
// tempShowImg = currImg;
// cv::RNG rng(time(0));
// for (size_t i = 0; i < points2.size(); i++) {
//     cv::circle(tempShowImg, points2[i], 2, cv::Scalar::all(255));
// }
// cv::imwrite("tempShow.jpg", tempShowImg);

// estimate geometry rigid 2D transformation matrix
#ifdef USE_OPENCV_FUNCTION
    cv::Mat inliers = cv::Mat::zeros(matches.size(), 1, CV_8U);
    cv::Mat optimalAffineMat = estimateAffinePartial2D(points1, points2, inliers, cv::RANSAC);
    cv::Mat R = optimalAffineMat.rowRange(0, 2).colRange(0, 2);
    double s = std::sqrt(cv::determinant(R));
    cv::Mat rigidtform2dR = R.mul(1.0 / s);
    cv::hconcat(rigidtform2dR, optimalAffineMat.col(2), relTform);
#else
    coder::array<double, 2U> pts1_tmp, pts2_tmp;
    coder::array<boolean_T, 2U> inlierIndex;
    double tform2x3[6];
    int status;

    pts1_tmp.set_size(matches.size(), 2);
    pts2_tmp.set_size(matches.size(), 2);
    for (size_t i = 0; i < matches.size(); i++) {
        pts1_tmp[i] = points1[i].x;
        pts1_tmp[i + matches.size()] = points1[i].y;

        pts2_tmp[i] = points2[i].x;
        pts2_tmp[i + matches.size()] = points2[i].y;
    }
    estimateAffineRigid2D::estimateAffineRigid2D(pts1_tmp, pts2_tmp, tform2x3,
                                                 inlierIndex, &status);
    cv::Mat inliers = cv::Mat(matches.size(), 1, CV_8U, inlierIndex.data());
    relTform = (cv::Mat_<double>(2, 3) << tform2x3[0], tform2x3[2], tform2x3[4], tform2x3[1], tform2x3[3], tform2x3[5]);
#endif

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

    // cv::imwrite("bigImgCopy1.jpg", HDmapOutput.bigImg);
    cv::copyMakeBorder(HDmapOutput.bigImg, HDmapOutput.bigImg, std::round(pad_top), std::round(pad_down), std::round(pad_left), std::round(pad_right), cv::BORDER_CONSTANT, cv::Scalar::all(0));
    // cv::imwrite("bigImgCopy2.jpg", HDmapOutput.bigImg);

    cv::resize(HDmapOutput.bigImg, HDmapOutput.bigImg, cv::Size(topImg.cols, topImg.rows));
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
    cv::imwrite("bigImgCopy.jpg", HDmapOutput.bigImg);
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