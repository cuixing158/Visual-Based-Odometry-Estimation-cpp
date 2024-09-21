///
/// @file           : Sim3Helpers.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "Sim3Helpers.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_data.h"
#include <cmath>
#include <cstring>

/// Function Definitions
///
/// @fn             : multiplyLogSim3
/// @brief          :
/// @param          : const double S1[16]
///                   const double S2[16]
///                   const double S3[16]
///                   double e[7]
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
void Sim3Helpers::multiplyLogSim3(const double S1[16], const double S2[16],
                                  const double S3[16], double e[7]) {
    double S12[16];
    double S123[16];
    double omega2[9];
    double w[9];
    double omega1[3];
    double a;
    double a1;
    double b1;
    double c;
    double d;
    double sigma;
    int b_omega2_tmp;
    int omega2_tmp;
    int r1;
    int r2;
    int r3;
    int rtemp;
    for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
        a1 = 0.0;
        for (int k{0}; k < 3; k++) {
            rtemp = k << 2;
            omega2[omega2_tmp + 3 * k] =
                (S1[omega2_tmp] * S2[rtemp] + S1[omega2_tmp + 4] * S2[rtemp + 1]) +
                S1[omega2_tmp + 8] * S2[rtemp + 2];
            a1 += S1[15] * S1[omega2_tmp + rtemp] * S2[k + 12];
        }
        omega1[omega2_tmp] = a1 + S1[omega2_tmp + 12];
    }
    for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
        rtemp = omega2_tmp << 2;
        S12[rtemp] = omega2[3 * omega2_tmp];
        S12[rtemp + 1] = omega2[3 * omega2_tmp + 1];
        S12[rtemp + 2] = omega2[3 * omega2_tmp + 2];
        S12[omega2_tmp + 12] = omega1[omega2_tmp];
    }
    S12[3] = 0.0;
    S12[7] = 0.0;
    S12[11] = 0.0;
    S12[15] = S1[15] * S2[15];
    for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
        a1 = 0.0;
        for (int k{0}; k < 3; k++) {
            rtemp = k << 2;
            omega2[omega2_tmp + 3 * k] =
                (S12[omega2_tmp] * S3[rtemp] + S12[omega2_tmp + 4] * S3[rtemp + 1]) +
                S12[omega2_tmp + 8] * S3[rtemp + 2];
            a1 += S12[15] * S12[omega2_tmp + rtemp] * S3[k + 12];
        }
        omega1[omega2_tmp] = a1 + S12[omega2_tmp + 12];
    }
    for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
        rtemp = omega2_tmp << 2;
        S123[rtemp] = omega2[3 * omega2_tmp];
        S123[rtemp + 1] = omega2[3 * omega2_tmp + 1];
        S123[rtemp + 2] = omega2[3 * omega2_tmp + 2];
        S123[omega2_tmp + 12] = omega1[omega2_tmp];
    }
    S123[3] = 0.0;
    S123[7] = 0.0;
    S123[11] = 0.0;
    S123[15] = S12[15] * S3[15];
    sigma = std::log(S123[15]);
    d = 0.5 * (((S123[0] + S123[5]) + S123[10]) - 1.0);
    if (std::abs(sigma) < 1.0E-5) {
        c = 1.0 - sigma / 2.0;
        if (d > 0.99999) {
            for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
                b_omega2_tmp = omega2_tmp << 2;
                omega2[3 * omega2_tmp] = S123[b_omega2_tmp] - S123[omega2_tmp];
                omega2[3 * omega2_tmp + 1] =
                    S123[b_omega2_tmp + 1] - S123[omega2_tmp + 4];
                omega2[3 * omega2_tmp + 2] =
                    S123[b_omega2_tmp + 2] - S123[omega2_tmp + 8];
            }
            omega1[0] = 0.5 * omega2[5];
            omega1[1] = 0.5 * omega2[6];
            omega1[2] = 0.5 * omega2[1];
            omega2[0] = 0.0;
            omega2[3] = -omega1[2];
            omega2[6] = omega1[1];
            omega2[1] = omega1[2];
            omega2[4] = 0.0;
            omega2[7] = -omega1[0];
            omega2[2] = -omega1[1];
            omega2[5] = omega1[0];
            omega2[8] = 0.0;
            a = 0.5;
            a1 = 0.16666666666666666;
        } else {
            double th;
            double thSquare;
            th = std::acos(d);
            thSquare = th * th;
            a = th / (2.0 * std::sqrt(1.0 - d * d));
            for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
                b_omega2_tmp = omega2_tmp << 2;
                omega2[3 * omega2_tmp] = S123[b_omega2_tmp] - S123[omega2_tmp];
                omega2[3 * omega2_tmp + 1] =
                    S123[b_omega2_tmp + 1] - S123[omega2_tmp + 4];
                omega2[3 * omega2_tmp + 2] =
                    S123[b_omega2_tmp + 2] - S123[omega2_tmp + 8];
            }
            omega1[0] = a * omega2[5];
            omega1[1] = a * omega2[6];
            omega1[2] = a * omega2[1];
            omega2[0] = 0.0;
            omega2[3] = -omega1[2];
            omega2[6] = omega1[1];
            omega2[1] = omega1[2];
            omega2[4] = 0.0;
            omega2[7] = -omega1[0];
            omega2[2] = -omega1[1];
            omega2[5] = omega1[0];
            omega2[8] = 0.0;
            a = (1.0 - std::cos(th)) / thSquare;
            a1 = (th - std::sin(th)) / (thSquare * th);
        }
    } else {
        c = (S123[15] - 1.0) / sigma;
        if (d > 0.99999) {
            a1 = sigma * sigma;
            for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
                b_omega2_tmp = omega2_tmp << 2;
                omega2[3 * omega2_tmp] = S123[b_omega2_tmp] - S123[omega2_tmp];
                omega2[3 * omega2_tmp + 1] =
                    S123[b_omega2_tmp + 1] - S123[omega2_tmp + 4];
                omega2[3 * omega2_tmp + 2] =
                    S123[b_omega2_tmp + 2] - S123[omega2_tmp + 8];
            }
            omega1[0] = 0.5 * omega2[5];
            omega1[1] = 0.5 * omega2[6];
            omega1[2] = 0.5 * omega2[1];
            omega2[0] = 0.0;
            omega2[3] = -omega1[2];
            omega2[6] = omega1[1];
            omega2[1] = omega1[2];
            omega2[4] = 0.0;
            omega2[7] = -omega1[0];
            omega2[2] = -omega1[1];
            omega2[5] = omega1[0];
            omega2[8] = 0.0;
            a = ((sigma - 1.0) * S123[15] + 1.0) / a1;
            a1 = (((0.5 * a1 - sigma) + 1.0) * S123[15] - 1.0) / (a1 * sigma);
        } else {
            double th;
            th = std::acos(d);
            a = th / (2.0 * std::sqrt(1.0 - d * d));
            for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
                b_omega2_tmp = omega2_tmp << 2;
                omega2[3 * omega2_tmp] = S123[b_omega2_tmp] - S123[omega2_tmp];
                omega2[3 * omega2_tmp + 1] =
                    S123[b_omega2_tmp + 1] - S123[omega2_tmp + 4];
                omega2[3 * omega2_tmp + 2] =
                    S123[b_omega2_tmp + 2] - S123[omega2_tmp + 8];
            }
            double thSquare;
            omega1[0] = a * omega2[5];
            omega1[1] = a * omega2[6];
            omega1[2] = a * omega2[1];
            omega2[0] = 0.0;
            omega2[3] = -omega1[2];
            omega2[6] = omega1[1];
            omega2[1] = omega1[2];
            omega2[4] = 0.0;
            omega2[7] = -omega1[0];
            omega2[2] = -omega1[1];
            omega2[5] = omega1[0];
            omega2[8] = 0.0;
            thSquare = th * th;
            a1 = S123[15] * std::sin(th);
            b1 = S123[15] * std::cos(th);
            d = thSquare + sigma * sigma;
            a = (a1 * sigma + (1.0 - b1) * th) / (th * d);
            a1 = (c - ((b1 - 1.0) * sigma + a1 * th) / d) / thSquare;
        }
    }
    for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
        for (int k{0}; k < 3; k++) {
            rtemp = omega2_tmp + 3 * k;
            w[rtemp] = (a * omega2[rtemp] +
                        ((a1 * omega2[omega2_tmp] * omega2[3 * k] +
                          a1 * omega2[omega2_tmp + 3] * omega2[3 * k + 1]) +
                         a1 * omega2[omega2_tmp + 6] * omega2[3 * k + 2])) +
                       c * static_cast<double>(iv[rtemp]);
        }
    }
    r1 = 0;
    r2 = 1;
    r3 = 2;
    a1 = std::abs(w[0]);
    d = std::abs(w[1]);
    if (d > a1) {
        a1 = d;
        r1 = 1;
        r2 = 0;
    }
    if (std::abs(w[2]) > a1) {
        r1 = 2;
        r2 = 1;
        r3 = 0;
    }
    w[r2] /= w[r1];
    w[r3] /= w[r1];
    w[r2 + 3] -= w[r2] * w[r1 + 3];
    w[r3 + 3] -= w[r3] * w[r1 + 3];
    w[r2 + 6] -= w[r2] * w[r1 + 6];
    w[r3 + 6] -= w[r3] * w[r1 + 6];
    if (std::abs(w[r3 + 3]) > std::abs(w[r2 + 3])) {
        rtemp = r2;
        r2 = r3;
        r3 = rtemp;
    }
    w[r3 + 3] /= w[r2 + 3];
    w[r3 + 6] -= w[r3 + 3] * w[r2 + 6];
    a1 = S123[12];
    d = S123[13];
    b1 = S123[14];
    for (int k{0}; k < 3; k++) {
        b_omega2_tmp = k + 3 * r1;
        omega2[b_omega2_tmp] = static_cast<double>(iv[k]) / w[r1];
        rtemp = k + 3 * r2;
        omega2[rtemp] =
            static_cast<double>(iv[k + 3]) - omega2[b_omega2_tmp] * w[r1 + 3];
        omega2_tmp = k + 3 * r3;
        omega2[omega2_tmp] =
            static_cast<double>(iv[k + 6]) - omega2[b_omega2_tmp] * w[r1 + 6];
        omega2[rtemp] /= w[r2 + 3];
        omega2[omega2_tmp] -= omega2[rtemp] * w[r2 + 6];
        omega2[omega2_tmp] /= w[r3 + 6];
        omega2[rtemp] -= omega2[omega2_tmp] * w[r3 + 3];
        omega2[b_omega2_tmp] -= omega2[omega2_tmp] * w[r3];
        omega2[b_omega2_tmp] -= omega2[rtemp] * w[r2];
        e[k] = (omega2[k] * a1 + omega2[k + 3] * d) + omega2[k + 6] * b1;
        e[k + 3] = omega1[k];
    }
    e[6] = sigma;
}

///
/// @fn             : sim3ToSform
/// @brief          :
/// @param          : const double minVecSim3[7]
///                   double S[16]
/// @return         : void
///
void Sim3Helpers::sim3ToSform(const double minVecSim3[7], double S[16]) {
    double R[9];
    double w[9];
    double wSquare[9];
    double a_tmp;
    double absxk;
    double c;
    double s;
    double scale;
    double t;
    double th;
    int S_tmp;
    w[0] = 0.0;
    w[3] = -minVecSim3[5];
    w[6] = minVecSim3[4];
    w[1] = minVecSim3[5];
    w[4] = 0.0;
    w[7] = -minVecSim3[3];
    w[2] = -minVecSim3[4];
    w[5] = minVecSim3[3];
    w[8] = 0.0;
    scale = 3.3121686421112381E-170;
    absxk = std::abs(minVecSim3[3]);
    if (absxk > 3.3121686421112381E-170) {
        th = 1.0;
        scale = absxk;
    } else {
        t = absxk / 3.3121686421112381E-170;
        th = t * t;
    }
    absxk = std::abs(minVecSim3[4]);
    if (absxk > scale) {
        t = scale / absxk;
        th = th * t * t + 1.0;
        scale = absxk;
    } else {
        t = absxk / scale;
        th += t * t;
    }
    absxk = std::abs(minVecSim3[5]);
    if (absxk > scale) {
        t = scale / absxk;
        th = th * t * t + 1.0;
        scale = absxk;
    } else {
        t = absxk / scale;
        th += t * t;
    }
    th = scale * std::sqrt(th);
    s = std::exp(minVecSim3[6]);
    for (int i{0}; i < 3; i++) {
        for (S_tmp = 0; S_tmp < 3; S_tmp++) {
            wSquare[i + 3 * S_tmp] =
                (w[i] * w[3 * S_tmp] + w[i + 3] * w[3 * S_tmp + 1]) +
                w[i + 6] * w[3 * S_tmp + 2];
        }
    }
    if (std::abs(minVecSim3[6]) < 1.0E-5) {
        c = 1.0;
        if (th < 1.0E-5) {
            a_tmp = 0.5;
            t = 0.16666666666666666;
            std::memset(&R[0], 0, 9U * sizeof(double));
            R[0] = 1.0;
            R[4] = 1.0;
            R[8] = 1.0;
            for (int i{0}; i < 9; i++) {
                R[i] = (R[i] + w[i]) + wSquare[i] / 2.0;
            }
        } else {
            double thSquare;
            thSquare = th * th;
            a_tmp = (1.0 - std::cos(th)) / thSquare;
            scale = std::sin(th);
            t = (th - scale) / (thSquare * th);
            absxk = scale / th;
            std::memset(&R[0], 0, 9U * sizeof(double));
            R[0] = 1.0;
            R[4] = 1.0;
            R[8] = 1.0;
            for (int i{0}; i < 9; i++) {
                R[i] = (R[i] + absxk * w[i]) + a_tmp * wSquare[i];
            }
        }
    } else {
        c = (s - 1.0) / minVecSim3[6];
        if (th < 1.0E-5) {
            scale = minVecSim3[6] * minVecSim3[6];
            a_tmp = ((minVecSim3[6] - 1.0) * s + 1.0) / scale;
            t = (((0.5 * scale - minVecSim3[6]) + 1.0) * s - 1.0) /
                (scale * minVecSim3[6]);
            std::memset(&R[0], 0, 9U * sizeof(double));
            R[0] = 1.0;
            R[4] = 1.0;
            R[8] = 1.0;
            for (int i{0}; i < 9; i++) {
                R[i] = (R[i] + w[i]) + wSquare[i] / 2.0;
            }
        } else {
            double a1_tmp;
            double b1_tmp;
            double thSquare;
            a1_tmp = std::sin(th);
            scale = s * a1_tmp;
            b1_tmp = std::cos(th);
            absxk = s * b1_tmp;
            thSquare = th * th;
            t = thSquare + minVecSim3[6] * minVecSim3[6];
            a_tmp = (scale * minVecSim3[6] + (1.0 - absxk) * th) / (th * t);
            t = (c - ((absxk - 1.0) * minVecSim3[6] + scale * th) / t) / thSquare;
            absxk = a1_tmp / th;
            scale = (1.0 - b1_tmp) / thSquare;
            std::memset(&R[0], 0, 9U * sizeof(double));
            R[0] = 1.0;
            R[4] = 1.0;
            R[8] = 1.0;
            for (int i{0}; i < 9; i++) {
                R[i] = (R[i] + absxk * w[i]) + scale * wSquare[i];
            }
        }
    }
    for (int i{0}; i < 9; i++) {
        w[i] = (a_tmp * w[i] + t * wSquare[i]) + c * static_cast<double>(iv[i]);
    }
    scale = minVecSim3[0];
    absxk = minVecSim3[1];
    t = minVecSim3[2];
    for (int i{0}; i < 3; i++) {
        S_tmp = i << 2;
        S[S_tmp] = R[3 * i];
        S[S_tmp + 1] = R[3 * i + 1];
        S[S_tmp + 2] = R[3 * i + 2];
        S[i + 12] = (w[i] * scale + w[i + 3] * absxk) + w[i + 6] * t;
        S[S_tmp + 3] = 0.0;
    }
    S[15] = s;
}

}  // namespace internal
}  // namespace core
}  // namespace robotics
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for Sim3Helpers.cpp
///
/// [EOF]
///
