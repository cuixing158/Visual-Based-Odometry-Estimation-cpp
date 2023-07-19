#include <iostream>
// #include <math>
#include "opencv2/opencv.hpp"

void subsetPoints(std::vector<cv::KeyPoint> keyPts, int N, std::vector<cv::KeyPoint> &outputPts) {
    outputPts.clear();
    if (keyPts.size() > N) {
        for (size_t i = 0; i < N; i++) {
            outputPts.push_back(keyPts[i]);
        }
    }
}

int main(int, char **) {
    cv::Mat srcImg = cv::imread("/opt_disk2/rd22946/MATLAB/R2023a/toolbox/matlab/imagesci/peppers.png", 0);
    cv::Ptr<cv::FeatureDetector> orbDetector = cv::ORB::create(srcImg.rows * srcImg.cols);

    cv::Mat currDescriptors;
    std::vector<cv::KeyPoint> keyPts, currKeypts;
    orbDetector->detect(srcImg, keyPts);
    // selectUniform(keyPts, Descriptions, 2000, currKeypts, currDescriptors);
    std::vector<int> indexs;
    // currKeypts.clear();
    // currDescriptors = cv::Mat();
    // ssc(keyPts, 2000, 0.1, currImg.cols, currImg.rows, currKeypts, indexs);
    // for (size_t i = 0; i < indexs.size(); i++) {
    //     currDescriptors.push_back(Descriptions.row(indexs[i]));
    // }
    subsetPoints(keyPts, 20, currKeypts);

    // for (size_t i = 0; i < 10; i++) {
    //     std::cout << currKeypts[i].pt << std::endl;
    // }

    orbDetector->compute(srcImg, currKeypts, currDescriptors);
    std::cout << currDescriptors.size << std::endl;
}