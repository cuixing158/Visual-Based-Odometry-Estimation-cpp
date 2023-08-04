// https://stackoverflow.com/questions/43830849/opencv-use-flann-with-orb-descriptors-to-match-features
#include <iostream>
// #include <math>
#include "opencv2/opencv.hpp"
cv::Mat src;

void myfun(cv::Mat& out) {
        src = (cv::Mat_<double>(2, 3) << 1, 2, 3, 4, 5, 9);
    out = cv::Mat(2, 3, CV_64F, src.data);
    src.release();
}
int main(int argc, char** argv) {
    // Read both images.
    cv::Mat image1 = cv::imread("/opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/map_R_new_undistort/map_R/1.jpg", cv::IMREAD_GRAYSCALE);
    if (image1.empty()) {
        std::cerr << "Couldn't read image in " << argv[1] << std::endl;
        return 1;
    }
    cv::Mat image2 = cv::imread("/opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/map_R_new_undistort/map_R/2.jpg", cv::IMREAD_GRAYSCALE);
    if (image2.empty()) {
        std::cerr << "Couldn't read image in " << argv[2] << std::endl;
        return 1;
    }

    cv::Mat image3;
    myfun(image3);
    std::cout << image3 << std::endl;
    // cv::imshow("Matches", image_matches);
    return 0;
}