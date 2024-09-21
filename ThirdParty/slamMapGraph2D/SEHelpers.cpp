///
/// @file           : SEHelpers.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "SEHelpers.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_rtwutil.h"
#include "sqrt.h"
#include "svd.h"
#include <algorithm>
#include <cmath>
#include <cstring>

/// Function Definitions
///
/// @fn             : expSE3hat
/// @brief          :
/// @param          : const double e[6]
///                   double T[16]
/// @return         : void
///
namespace SlamGraph2D {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
void SEHelpers::expSE3hat(const double e[6], double T[16]) {
    double Sphi[9];
    double V[9];
    double absxk;
    double c;
    double c_tmp;
    double d;
    double d1;
    double d2;
    double scale;
    double t;
    double theta;
    int T_tmp;
    Sphi[0] = 0.0;
    Sphi[3] = -e[5];
    Sphi[6] = e[4];
    Sphi[1] = e[5];
    Sphi[4] = 0.0;
    Sphi[7] = -e[3];
    Sphi[2] = -e[4];
    Sphi[5] = e[3];
    Sphi[8] = 0.0;
    scale = 3.3121686421112381E-170;
    absxk = std::abs(e[3]);
    if (absxk > 3.3121686421112381E-170) {
        theta = 1.0;
        scale = absxk;
    } else {
        t = absxk / 3.3121686421112381E-170;
        theta = t * t;
    }
    absxk = std::abs(e[4]);
    if (absxk > scale) {
        t = scale / absxk;
        theta = theta * t * t + 1.0;
        scale = absxk;
    } else {
        t = absxk / scale;
        theta += t * t;
    }
    absxk = std::abs(e[5]);
    if (absxk > scale) {
        t = scale / absxk;
        theta = theta * t * t + 1.0;
        scale = absxk;
    } else {
        t = absxk / scale;
        theta += t * t;
    }
    theta = scale * std::sqrt(theta);
    scale = theta * theta;
    absxk = 1.0 - std::cos(theta);
    t = absxk / scale;
    c_tmp = std::sin(theta);
    c = (theta - c_tmp) / (scale * theta);
    if ((t < 2.2204460492503131E-16) || (std::isinf(t) || std::isnan(t))) {
        std::memset(&V[0], 0, 9U * sizeof(double));
        V[0] = 1.0;
        V[4] = 1.0;
        V[8] = 1.0;
        std::copy(&V[0], &V[9], &Sphi[0]);
    } else {
        double b_Sphi[9];
        signed char V_tmp[9];
        for (int i{0}; i < 9; i++) {
            V_tmp[i] = 0;
        }
        V_tmp[0] = 1;
        V_tmp[4] = 1;
        V_tmp[8] = 1;
        for (int i{0}; i < 3; i++) {
            d = Sphi[i + 3];
            d1 = Sphi[i + 6];
            for (T_tmp = 0; T_tmp < 3; T_tmp++) {
                double d3;
                double d4;
                int b_V_tmp;
                d2 = Sphi[3 * T_tmp];
                d3 = Sphi[3 * T_tmp + 1];
                d4 = Sphi[3 * T_tmp + 2];
                b_V_tmp = i + 3 * T_tmp;
                V[b_V_tmp] = (static_cast<double>(V_tmp[b_V_tmp]) + t * Sphi[b_V_tmp]) +
                             ((c * Sphi[i] * d2 + c * d * d3) + c * d1 * d4);
                b_Sphi[b_V_tmp] = ((Sphi[i] * d2 + d * d3) + d1 * d4) / scale;
            }
        }
        for (int i{0}; i < 9; i++) {
            Sphi[i] = (static_cast<double>(V_tmp[i]) + Sphi[i] * c_tmp / theta) +
                      b_Sphi[i] * absxk;
        }
    }
    d = e[0];
    d1 = e[1];
    d2 = e[2];
    for (int i{0}; i < 3; i++) {
        T_tmp = i << 2;
        T[T_tmp] = Sphi[3 * i];
        T[T_tmp + 1] = Sphi[3 * i + 1];
        T[T_tmp + 2] = Sphi[3 * i + 2];
        T[i + 12] = (V[i] * d + V[i + 3] * d1) + V[i + 6] * d2;
    }
    T[3] = 0.0;
    T[7] = 0.0;
    T[11] = 0.0;
    T[15] = 1.0;
}

///
/// @fn             : veelogmSE3
/// @brief          :
/// @param          : const double T[16]
///                   double vec[6]
/// @return         : void
///
void SEHelpers::veelogmSE3(const double T[16], double vec[6]) {
    creal_T u;
    double Vinv[9];
    double b_I[9];
    double wx[9];
    double wv[3];
    double a;
    double theta;
    double thetaSq;
    int wx_tmp;
    bool guard1{false};
    bool y;
    a = 0.5 * (((T[0] + T[5]) + T[10]) - 1.0);
    if (!(std::abs(a) > 1.0)) {
        u.re = std::acos(a);
    } else {
        creal_T v;
        v.re = a + 1.0;
        v.im = 0.0;
        ::SlamGraph2D::coder::internal::scalar::b_sqrt(v);
        u.re = 1.0 - a;
        u.im = 0.0;
        ::SlamGraph2D::coder::internal::scalar::b_sqrt(u);
        a = u.re;
        u.re = 2.0 * rt_atan2d_snf(a, v.re);
    }
    a = u.re / std::sin(u.re);
    for (int i{0}; i < 3; i++) {
        wx_tmp = i << 2;
        wx[3 * i] = T[wx_tmp] - T[i];
        wx[3 * i + 1] = T[wx_tmp + 1] - T[i + 4];
        wx[3 * i + 2] = T[wx_tmp + 2] - T[i + 8];
    }
    wv[0] = wx[5];
    wv[1] = wx[6];
    wv[2] = wx[1];
    guard1 = false;
    if ((!std::isinf(a)) && (!std::isnan(a))) {
        bool exitg1;
        y = true;
        wx_tmp = 0;
        exitg1 = false;
        while ((!exitg1) && (wx_tmp < 3)) {
            if (!(wv[wx_tmp] == 0.0)) {
                y = false;
                exitg1 = true;
            } else {
                wx_tmp++;
            }
        }
        if (!y) {
            wv[0] = wx[5] * a / 2.0;
            wv[1] = wx[6] * a / 2.0;
            wv[2] = wx[1] * a / 2.0;
        } else {
            guard1 = true;
        }
    } else {
        guard1 = true;
    }
    if (guard1) {
        std::memset(&b_I[0], 0, 9U * sizeof(double));
        b_I[0] = 1.0;
        b_I[4] = 1.0;
        b_I[8] = 1.0;
        for (int i{0}; i < 3; i++) {
            int I_tmp;
            wx_tmp = i << 2;
            b_I[3 * i] -= T[wx_tmp];
            I_tmp = 3 * i + 1;
            b_I[I_tmp] -= T[wx_tmp + 1];
            I_tmp = 3 * i + 2;
            b_I[I_tmp] -= T[wx_tmp + 2];
        }
        y = true;
        for (wx_tmp = 0; wx_tmp < 9; wx_tmp++) {
            if (y) {
                a = b_I[wx_tmp];
                if (std::isinf(a) || std::isnan(a)) {
                    y = false;
                }
            } else {
                y = false;
            }
        }
        if (y) {
            ::SlamGraph2D::coder::internal::svd(b_I, Vinv, wv, wx);
        } else {
            for (int i{0}; i < 9; i++) {
                wx[i] = rtNaN;
            }
        }
        a = 1.0 / std::sqrt((wx[6] * wx[6] + wx[7] * wx[7]) + wx[8] * wx[8]);
        wv[0] = wx[6] * a * u.re;
        wv[1] = wx[7] * a * u.re;
        wv[2] = wx[8] * a * u.re;
    }
    theta = std::sqrt((wv[0] * wv[0] + wv[1] * wv[1]) + wv[2] * wv[2]);
    thetaSq = theta * theta;
    a = (1.0 - std::cos(theta)) / thetaSq;
    if ((a < 2.2204460492503131E-16) || (std::isinf(a) || std::isnan(a))) {
        std::memset(&Vinv[0], 0, 9U * sizeof(double));
        Vinv[0] = 1.0;
        Vinv[4] = 1.0;
        Vinv[8] = 1.0;
    } else {
        wx[0] = 0.0;
        wx[3] = -wv[2];
        wx[6] = wv[1];
        wx[1] = wv[2];
        wx[4] = 0.0;
        wx[7] = -wv[0];
        wx[2] = -wv[1];
        wx[5] = wv[0];
        wx[8] = 0.0;
        a = 1.0 / thetaSq * (1.0 - std::sin(theta) / theta / (2.0 * a));
        std::memset(&b_I[0], 0, 9U * sizeof(double));
        for (wx_tmp = 0; wx_tmp < 3; wx_tmp++) {
            b_I[wx_tmp + 3 * wx_tmp] = 1.0;
            for (int i{0}; i < 3; i++) {
                Vinv[wx_tmp + 3 * i] =
                    (wx[wx_tmp] * wx[3 * i] + wx[wx_tmp + 3] * wx[3 * i + 1]) +
                    wx[wx_tmp + 6] * wx[3 * i + 2];
            }
        }
        for (int i{0}; i < 9; i++) {
            Vinv[i] = (b_I[i] - 0.5 * wx[i]) + a * Vinv[i];
        }
    }
    a = T[12];
    theta = T[13];
    thetaSq = T[14];
    for (int i{0}; i < 3; i++) {
        vec[i] = (Vinv[i] * a + Vinv[i + 3] * theta) + Vinv[i + 6] * thetaSq;
        vec[i + 3] = wv[i];
    }
}

}  // namespace internal
}  // namespace core
}  // namespace robotics
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for SEHelpers.cpp
///
/// [EOF]
///
