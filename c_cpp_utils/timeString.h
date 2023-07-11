/**
 * @file           : timeString.h
 * @target         : windows or linux or mac os
 * @details        : get current timestamp string inux,Windows,MacOS
 * @author         : cuixingxing
 * @email          : xingxing.cui@long-horn.com
 * @date           : 17-Mar-2022 12:02:52
 * @version        : V1.0.0
 * @copyright      : Copyright (C) 2022 cuixingxing.All rights reserved.
 */

#pragma once
#include <time.h>

#include <iostream>

#include "opencv2/opencv.hpp"
// 0月日时分秒毫秒
// 1月日时分
// 2月日

inline std::string timeStamp(int onlyMDH = 0) {
    static int countN = 0;
    time_t curtime = time(0);
    tm tim = *localtime(&curtime);
    int year, month, day, hour, minute, second, msecond;
    year = 1900 + tim.tm_year;
    month = tim.tm_mon + 1;
    day = tim.tm_mday;
    hour = tim.tm_hour;
    minute = tim.tm_min;
    second = tim.tm_sec;
    long long time_ = (long long)cv::getTickCount();
    // msecond = time_ % (1000000);
    msecond = ++countN;
    if (onlyMDH == 1) {
        return "_M" + std::to_string(month) + "D" + std::to_string(day) + "H" +
               std::to_string(hour) + "m" + std::to_string(minute) + "_";
    }
    if (onlyMDH == 2) {
        return "_M" + std::to_string(month) + "D" + std::to_string(day) + "_";
    }
    return "_Y" + std::to_string(year) + "M" + std::to_string(month) + "D" +
           std::to_string(day) + "H" + std::to_string(hour) + "m" +
           std::to_string(minute) + "s" + std::to_string(second) + "m" +
           std::to_string(msecond) + "_";
}

inline char* gettm(int64 timestamp) {
    int64 milli =
        timestamp +
        (int64)8 * 60 * 60 *
            1000;  //此处转化为东八区北京时间，如果是其它时区需要按需求修改
    auto mTime = std::chrono::milliseconds(milli);
    auto tp = std::chrono::time_point<std::chrono::system_clock,
                                      std::chrono::milliseconds>(mTime);
    auto tt = std::chrono::system_clock::to_time_t(tp);
    std::tm* now = std::gmtime(&tt);
    char strtime[255];
    sprintf(strtime, "%4d_%02d_%02d__%02d_%02d_%02d_\n", now->tm_year + 1900,
            now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min,
            now->tm_sec);
    return strtime;
}

inline std::time_t getTimeStamp() {
    std::chrono::time_point<std::chrono::system_clock, std::chrono::milliseconds>
        tp = std::chrono::time_point_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now());
    auto tmp = std::chrono::duration_cast<std::chrono::milliseconds>(
        tp.time_since_epoch());
    std::time_t timestamp = tmp.count();
    // std::time_t timestamp = std::chrono::system_clock::to_time_t(tp);
    return timestamp;  // 毫秒级时间戳
}