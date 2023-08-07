// https://stackoverflow.com/questions/43830849/opencv-use-flann-with-orb-descriptors-to-match-features
#include <iostream>
#include "DBoW3.h"
#include "c_cpp_utils/path.h"
// #include <math>
#include "opencv2/opencv.hpp"

using namespace DBoW3;

void myfun(cv::Mat& out) {
    cv::Mat src = (cv::Mat_<double>(2, 3) << 1, 2, 3, 4, 5, 9);
    out = cv::Mat(2, 3, CV_64F, src.data);
    src.release();
}

int main(int argc, char** argv) {
    std::string imagePath = "/opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/map_R_new_undistort/map_R";  // opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/116_new_undistort/116";  //"/opt_disk2/rd22946/AllDataAndModels/from_tongwenchao/map_R_new_undistort/map_R";
    std::vector<std::string> imagePaths;
    size_t numImgs = getFullNames(filesystem::path(imagePath), imagePaths, ".jpg");
    // std::sort(imagePaths.begin(), imagePaths.end(),
    //           [](std::string p1, std::string p2) { return atoi(filesystem::path(p1).filenameNoExt().substr(6).c_str()) < atoi(filesystem::path(p2).filenameNoExt().substr(6).c_str()); });
    std::sort(imagePaths.begin(), imagePaths.end(),
              [](std::string p1, std::string p2) { return atoi(filesystem::path(p1).filenameNoExt().c_str()) < atoi(filesystem::path(p2).filenameNoExt().c_str()); });

    std::vector<cv::KeyPoint> keypts;
    cv::Mat feature;
    std::vector<cv::Mat> features;
    cv::Ptr<cv::ORB> orbDetector = cv::ORB::create();
    for (size_t i = 0; i < 200; i++) {
        cv::Mat srcImg = cv::imread(imagePaths[i], cv::IMREAD_COLOR);
        cv::cvtColor(srcImg, srcImg, cv::COLOR_BGR2GRAY);

        orbDetector->detectAndCompute(srcImg, cv::noArray(), keypts, feature);
        features.push_back(feature);
    }
    std::string saveDataBaseYmlGz = "./database.yml.gz";

    // branching factor and depth levels
    const int k = 10;
    const int L = 4;
    const DBoW3::WeightingType weight = DBoW3::TF_IDF;
    const DBoW3::ScoringType score = DBoW3::L1_NORM;

    std::cout << "From features,Create vocabulary,please wait ..." << std::endl;
    DBoW3::Vocabulary voc(k, L, weight, score);

    voc.create(features);
    Database db;
    db.setVocabulary(voc, false, 0);  // false = do not use direct index
    // (so ignore the last param)
    // The direct index is useful if we want to retrieve the features that
    // belong to some vocabulary node.
    // db creates a copy of the vocabulary, we may get rid of "voc" now

    // add images to the database
    for (size_t i = 0; i < features.size(); i++)
        db.add(features[i]);

    std::cout << "add features per image done!" << std::endl;

    std::string databaseFile = saveDataBaseYmlGz;
    std::cout << "Vocabulary information: " << std::endl
              << voc << std::endl
              << "have saved this path:" << databaseFile << std::endl;
    db.save(databaseFile);

    db.load(databaseFile);

    return 0;
}