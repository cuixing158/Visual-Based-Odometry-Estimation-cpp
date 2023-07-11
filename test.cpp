#include <iostream>
// #include <math>
#include "opencv2/opencv.hpp"

int main(int, char **) {
    cv::Mat srcImg = cv::imread("/opt_disk2/rd22946/MATLAB/R2023a/toolbox/matlab/imagesci/peppers.png");
    double theta = 3.141592653;
    cv::Mat tform = (cv::Mat_<double>(2, 3) << std::cos(theta / 6), -std::sin(theta / 6), 100, std::sin(theta / 6), std::cos(theta / 6), 200);
    cv::Mat dstImg;
    cv::warpAffine(srcImg, dstImg, tform, cv::Size(512, 384));
    cv::imwrite("dstImg.jpg", dstImg);
}