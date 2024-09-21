///
/// @file           : PoseGraphHelpers.cpp
/// @target         : Texas Instruments->C6000
/// @details        : pose graph algorithms
/// @author         : cuixingxing
/// @email          : cuixingxing150@gmail.com
/// @date           : 26-Jul-2023 07:45:22
/// @version        : V0.1.0
/// @copyright      : Copyright (C) 2023 TheMatrix Inc.All rights reserved.
///

/// @include file    : Include Files
#include "PoseGraphHelpers.h"
#include "BlockInserter2.h"
#include "BlockMatrix.h"
#include "SEHelpers.h"
#include "Sim3Helpers.h"
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "slamMapGraph2D_data.h"
#include "slamMapGraph2D_rtwutil.h"
#include "sparse.h"
#include "sparse1.h"
#include "coder_array.h"
#include <algorithm>
#include <cmath>
#include <cstring>

/// Function Declarations
namespace SlamGraph2D {
static int div_s32(int numerator, int denominator);

}

/// Function Definitions
///
/// @fn             : costBetweenTwoNodes
/// @brief          :
/// @param          : const ::coder::array<double, 2U> &Toi
///                   const ::coder::array<double, 2U> &Toj
///                   const ::coder::array<double, 2U> &measurement
///                   const ::coder::array<double, 2U> &Omega
///                   bool nodejIsLandmark
///                   double gradi_data[]
///                   int &gradi_size
///                   double gradj_data[]
///                   int &gradj_size
///                   double hessii_data[]
///                   int hessii_size[2]
///                   double hessij_data[]
///                   int hessij_size[2]
///                   double hessji_data[]
///                   int hessji_size[2]
///                   double hessjj_data[]
///                   int hessjj_size[2]
/// @return         : double
///
namespace SlamGraph2D {
namespace coder {
namespace nav {
namespace algs {
namespace internal {
double PoseGraphHelpers::costBetweenTwoNodes(
    const ::coder::array<double, 2U> &Toi,
    const ::coder::array<double, 2U> &Toj,
    const ::coder::array<double, 2U> &measurement,
    const ::coder::array<double, 2U> &Omega, bool nodejIsLandmark,
    double gradi_data[], int &gradi_size, double gradj_data[], int &gradj_size,
    double hessii_data[], int hessii_size[2], double hessij_data[],
    int hessij_size[2], double hessji_data[], int hessji_size[2],
    double hessjj_data[], int hessjj_size[2]) {
    static const signed char b_N[36]{1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                     0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
    double Jaci_data[49];
    double Jacj_data[49];
    double Jaci[36];
    double Jacj[36];
    double b_y_tmp_data[21];
    double y_tmp_data[21];
    double dv[18];
    double f_R[9];
    double e_data[7];
    double e_R[3];
    double cost;
    int Jaci_size[2];
    int Jacj_size[2];
    int b_y_tmp_size[2];
    int y_tmp_size[2];
    int aoffset;
    int b_i;
    int boffset;
    int coffset;
    int i;
    int j;
    if (Omega.size(0) == 2) {
        double R[9];
        double N[6];
        double cosTheta;
        double sinTheta;
        for (i = 0; i < 3; i++) {
            R[3 * i] = measurement[measurement.size(0) * i];
            R[3 * i + 1] = measurement[measurement.size(0) * i + 1];
            R[3 * i + 2] = measurement[measurement.size(0) * i + 2];
        }
        double b_measurement;
        double d;
        double q;
        double s;
        double sinThetaDy;
        d = rt_atan2d_snf(Toi[1], Toi[0]);
        sinTheta = std::sin(d);
        cosTheta = std::cos(d);
        q = Toj[Toj.size(0) * 2 + 1] - Toi[Toi.size(0) * 2 + 1];
        sinThetaDy = q * sinTheta;
        b_measurement = Toj[Toj.size(0) * 2] - Toi[Toi.size(0) * 2];
        s = b_measurement * cosTheta;
        b_i = 2;
        e_data[0] = (sinThetaDy + s) - R[6];
        q = q * cosTheta - b_measurement * sinTheta;
        e_data[1] = q - R[7];
        N[0] = -cosTheta;
        N[2] = -sinTheta;
        N[4] = q;
        N[1] = sinTheta;
        N[3] = -cosTheta;
        N[5] = -sinThetaDy - s;
        Jaci_size[0] = 2;
        Jaci_size[1] = 3;
        for (i = 0; i < 6; i++) {
            Jaci_data[i] = N[i];
        }
        Jacj_size[0] = 2;
        Jacj_size[1] = 2;
        Jacj_data[0] = cosTheta;
        Jacj_data[1] = -sinTheta;
        Jacj_data[2] = sinTheta;
        Jacj_data[3] = cosTheta;
    } else if (Omega.size(0) == 3) {
        if (nodejIsLandmark) {
            double Sio[16];
            double Tio[16];
            double deltatform[16];
            double R[9];
            double b_R[9];
            double t[3];
            double d;
            double d1;
            double d2;
            for (i = 0; i < 4; i++) {
                coffset = i << 2;
                Sio[coffset] = Toi[Toi.size(0) * i];
                deltatform[coffset] = measurement[measurement.size(0) * i];
                Sio[coffset + 1] = Toi[Toi.size(0) * i + 1];
                deltatform[coffset + 1] = measurement[measurement.size(0) * i + 1];
                Sio[coffset + 2] = Toi[Toi.size(0) * i + 2];
                deltatform[coffset + 2] = measurement[measurement.size(0) * i + 2];
                Sio[coffset + 3] = Toi[Toi.size(0) * i + 3];
                deltatform[coffset + 3] = measurement[measurement.size(0) * i + 3];
            }
            for (i = 0; i < 3; i++) {
                R[3 * i] = Sio[i];
                R[3 * i + 1] = Sio[i + 4];
                R[3 * i + 2] = Sio[i + 8];
            }
            for (i = 0; i < 9; i++) {
                b_R[i] = -R[i];
            }
            d = Sio[12];
            d1 = Sio[13];
            d2 = Sio[14];
            for (i = 0; i < 3; i++) {
                coffset = i << 2;
                Tio[coffset] = R[3 * i];
                Tio[coffset + 1] = R[3 * i + 1];
                Tio[coffset + 2] = R[3 * i + 2];
                Tio[i + 12] = (b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2;
            }
            Tio[3] = 0.0;
            Tio[7] = 0.0;
            Tio[11] = 0.0;
            Tio[15] = 1.0;
            for (i = 0; i < 3; i++) {
                f_R[3 * i] = -deltatform[i];
                f_R[3 * i + 1] = -deltatform[i + 4];
                f_R[3 * i + 2] = -deltatform[i + 8];
                t[i] = ((Tio[i] * Toj[Toj.size(0) * 3] +
                         Tio[i + 4] * Toj[Toj.size(0) * 3 + 1]) +
                        Tio[i + 8] * Toj[Toj.size(0) * 3 + 2]) +
                       Tio[i + 12];
            }
            b_i = 3;
            d = deltatform[12];
            d1 = deltatform[13];
            d2 = deltatform[14];
            for (i = 0; i < 3; i++) {
                dv[3 * i] = iv[3 * i];
                j = 3 * i + 1;
                dv[j] = iv[j];
                j = 3 * i + 2;
                dv[j] = iv[j];
                e_data[i] = t[i] + ((f_R[i] * d + f_R[i + 3] * d1) + f_R[i + 6] * d2);
            }
            dv[9] = -0.0;
            dv[12] = t[2];
            dv[15] = -t[1];
            dv[10] = -t[2];
            dv[13] = -0.0;
            dv[16] = t[0];
            dv[11] = t[1];
            dv[14] = -t[0];
            dv[17] = -0.0;
            Jaci_size[0] = 3;
            Jaci_size[1] = 6;
            std::copy(&dv[0], &dv[18], &Jaci_data[0]);
            Jacj_size[0] = 3;
            Jacj_size[1] = 3;
            for (i = 0; i < 3; i++) {
                coffset = i << 2;
                Jacj_data[3 * i] = Tio[coffset];
                Jacj_data[3 * i + 1] = Tio[coffset + 1];
                Jacj_data[3 * i + 2] = Tio[coffset + 2];
            }
        } else {
            double R[9];
            double b_R[9];
            double c_R[9];
            double dRoidtheta[4];
            double d_R[4];
            double y[4];
            double y_tmp[4];
            double b_y_tmp[2];
            double e_tmp[2];
            double b_measurement;
            double cosTheta;
            double d;
            double d1;
            double d2;
            double d3;
            double d4;
            double q;
            double sinTheta;
            for (i = 0; i < 3; i++) {
                R[3 * i] = Toi[Toi.size(0) * i];
                b_R[3 * i] = Toj[Toj.size(0) * i];
                c_R[3 * i] = measurement[measurement.size(0) * i];
                boffset = 3 * i + 1;
                R[boffset] = Toi[Toi.size(0) * i + 1];
                b_R[boffset] = Toj[Toj.size(0) * i + 1];
                c_R[boffset] = measurement[measurement.size(0) * i + 1];
                boffset = 3 * i + 2;
                R[boffset] = Toi[Toi.size(0) * i + 2];
                b_R[boffset] = Toj[Toj.size(0) * i + 2];
                c_R[boffset] = measurement[measurement.size(0) * i + 2];
            }
            for (i = 0; i < 2; i++) {
                d = c_R[3 * i];
                d1 = c_R[3 * i + 1];
                d2 = R[0] * d + R[3] * d1;
                d = R[1] * d + R[4] * d1;
                y[i] = d2 * b_R[0] + d * b_R[1];
                y[i + 2] = d2 * b_R[3] + d * b_R[4];
                e_tmp[i] = b_R[i + 6] - R[i + 6];
            }
            b_measurement = rt_atan2d_snf(y[1], y[0]);
            d_R[0] = R[0];
            d_R[1] = R[3];
            d_R[2] = R[1];
            d_R[3] = R[4];
            if (std::isnan(b_measurement + 3.1415926535897931) ||
                std::isinf(b_measurement + 3.1415926535897931)) {
                cosTheta = rtNaN;
            } else if (b_measurement + 3.1415926535897931 == 0.0) {
                cosTheta = 0.0;
            } else {
                bool rEQ0;
                cosTheta =
                    std::fmod(b_measurement + 3.1415926535897931, 6.2831853071795862);
                rEQ0 = (cosTheta == 0.0);
                if (!rEQ0) {
                    q = std::abs((b_measurement + 3.1415926535897931) /
                                 6.2831853071795862);
                    rEQ0 =
                        !(std::abs(q - std::floor(q + 0.5)) > 2.2204460492503131E-16 * q);
                }
                if (rEQ0) {
                    cosTheta = 0.0;
                } else if (b_measurement + 3.1415926535897931 < 0.0) {
                    cosTheta += 6.2831853071795862;
                }
            }
            dRoidtheta[0] = R[3];
            dRoidtheta[2] = -R[0];
            dRoidtheta[1] = R[0];
            dRoidtheta[3] = R[3];
            d = c_R[0];
            d1 = c_R[1];
            d2 = c_R[3];
            d3 = c_R[4];
            for (j = 0; j < 2; j++) {
                coffset = j << 1;
                b_measurement = dRoidtheta[j + 2];
                d4 = dRoidtheta[j];
                y[coffset] = d * d4 + d1 * b_measurement;
                y[coffset + 1] = d2 * d4 + d3 * b_measurement;
            }
            d = e_tmp[0];
            d1 = e_tmp[1];
            for (j = 0; j < 2; j++) {
                coffset = j << 1;
                b_measurement = R[j % 2];
                q = R[(j + 2) % 2 + 3];
                dRoidtheta[coffset] = c_R[0] * b_measurement + c_R[1] * q;
                dRoidtheta[coffset + 1] = c_R[3] * b_measurement + c_R[4] * q;
                b_y_tmp[j] = (d_R[j] * d + d_R[j + 2] * d1) - c_R[j + 6];
            }
            b_i = 3;
            e_data[0] = c_R[0] * b_y_tmp[0] + b_y_tmp[1] * c_R[1];
            e_data[1] = b_y_tmp[0] * c_R[3] + b_y_tmp[1] * c_R[4];
            e_data[2] = cosTheta - 3.1415926535897931;
            y_tmp[0] = -c_R[0];
            y_tmp[1] = -c_R[3];
            y_tmp[2] = -c_R[1];
            y_tmp[3] = -c_R[4];
            d = R[0];
            d1 = R[3];
            d2 = R[1];
            d3 = R[4];
            d4 = e_tmp[0];
            sinTheta = e_tmp[1];
            for (i = 0; i < 2; i++) {
                q = y_tmp[i + 2];
                b_measurement = y_tmp[i];
                d_R[i] = b_measurement * d + q * d1;
                d_R[i + 2] = b_measurement * d2 + q * d3;
                b_y_tmp[i] = y[i] * d4 + y[i + 2] * sinTheta;
            }
            f_R[0] = d_R[0];
            f_R[1] = d_R[1];
            f_R[6] = b_y_tmp[0];
            f_R[3] = d_R[2];
            f_R[4] = d_R[3];
            f_R[7] = b_y_tmp[1];
            f_R[2] = 0.0;
            f_R[5] = 0.0;
            f_R[8] = -1.0;
            Jaci_size[0] = 3;
            Jaci_size[1] = 3;
            std::copy(&f_R[0], &f_R[9], &Jaci_data[0]);
            f_R[0] = dRoidtheta[0];
            f_R[1] = dRoidtheta[1];
            f_R[6] = 0.0;
            f_R[3] = dRoidtheta[2];
            f_R[4] = dRoidtheta[3];
            f_R[7] = 0.0;
            f_R[2] = 0.0;
            f_R[5] = 0.0;
            f_R[8] = 1.0;
            Jacj_size[0] = 3;
            Jacj_size[1] = 3;
            std::copy(&f_R[0], &f_R[9], &Jacj_data[0]);
        }
    } else if (Omega.size(0) == 6) {
        double Sio[16];
        double Sji[16];
        double Tio[16];
        double Tji[16];
        double Tjo[16];
        double Tjop[16];
        double deltatform[16];
        double R[9];
        double b_R[9];
        double dv1[6];
        double d;
        double d1;
        double d2;
        double d3;
        for (i = 0; i < 4; i++) {
            coffset = i << 2;
            Sio[coffset] = Toi[Toi.size(0) * i];
            Sji[coffset] = Toj[Toj.size(0) * i];
            deltatform[coffset] = measurement[measurement.size(0) * i];
            Sio[coffset + 1] = Toi[Toi.size(0) * i + 1];
            Sji[coffset + 1] = Toj[Toj.size(0) * i + 1];
            deltatform[coffset + 1] = measurement[measurement.size(0) * i + 1];
            Sio[coffset + 2] = Toi[Toi.size(0) * i + 2];
            Sji[coffset + 2] = Toj[Toj.size(0) * i + 2];
            deltatform[coffset + 2] = measurement[measurement.size(0) * i + 2];
            Sio[coffset + 3] = Toi[Toi.size(0) * i + 3];
            Sji[coffset + 3] = Toj[Toj.size(0) * i + 3];
            deltatform[coffset + 3] = measurement[measurement.size(0) * i + 3];
        }
        for (i = 0; i < 3; i++) {
            R[3 * i] = Sio[i];
            R[3 * i + 1] = Sio[i + 4];
            R[3 * i + 2] = Sio[i + 8];
        }
        for (i = 0; i < 9; i++) {
            b_R[i] = -R[i];
        }
        d = Sio[12];
        d1 = Sio[13];
        d2 = Sio[14];
        for (i = 0; i < 3; i++) {
            coffset = i << 2;
            Tio[coffset] = R[3 * i];
            aoffset = 3 * i + 1;
            Tio[coffset + 1] = R[aoffset];
            boffset = 3 * i + 2;
            Tio[coffset + 2] = R[boffset];
            Tio[i + 12] = (b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2;
            R[3 * i] = deltatform[i];
            R[aoffset] = deltatform[i + 4];
            R[boffset] = deltatform[i + 8];
        }
        Tio[3] = 0.0;
        Tio[7] = 0.0;
        Tio[11] = 0.0;
        Tio[15] = 1.0;
        for (i = 0; i < 9; i++) {
            b_R[i] = -R[i];
        }
        d = deltatform[12];
        d1 = deltatform[13];
        d2 = deltatform[14];
        for (i = 0; i < 3; i++) {
            coffset = i << 2;
            Tji[coffset] = R[3 * i];
            Tji[coffset + 1] = R[3 * i + 1];
            Tji[coffset + 2] = R[3 * i + 2];
            Tji[i + 12] = (b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2;
        }
        Tji[3] = 0.0;
        Tji[7] = 0.0;
        Tji[11] = 0.0;
        Tji[15] = 1.0;
        for (i = 0; i < 4; i++) {
            d = Tji[i];
            d1 = Tji[i + 4];
            d2 = Tji[i + 8];
            d3 = Tji[i + 12];
            for (j = 0; j < 4; j++) {
                coffset = j << 2;
                deltatform[i + coffset] = ((d * Tio[coffset] + d1 * Tio[coffset + 1]) +
                                           d2 * Tio[coffset + 2]) +
                                          d3 * Tio[coffset + 3];
            }
        }
        for (i = 0; i < 3; i++) {
            R[3 * i] = Sji[i];
            R[3 * i + 1] = Sji[i + 4];
            R[3 * i + 2] = Sji[i + 8];
        }
        for (i = 0; i < 9; i++) {
            b_R[i] = -R[i];
        }
        d = Sji[12];
        d1 = Sji[13];
        d2 = Sji[14];
        for (i = 0; i < 3; i++) {
            coffset = i << 2;
            Tjo[coffset] = R[3 * i];
            Tjo[coffset + 1] = R[3 * i + 1];
            Tjo[coffset + 2] = R[3 * i + 2];
            Tjo[i + 12] = (b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2;
        }
        Tjo[3] = 0.0;
        Tjo[7] = 0.0;
        Tjo[11] = 0.0;
        Tjo[15] = 1.0;
        for (b_i = 0; b_i < 6; b_i++) {
            double Tm[16];
            double g_R[16];
            double c_R[9];
            double N[6];
            for (i = 0; i < 6; i++) {
                N[i] = static_cast<double>(b_N[i + 6 * b_i]) * 1.0E-5;
            }
            robotics::core::internal::SEHelpers::expSE3hat(N, Sio);
            for (i = 0; i < 6; i++) {
                N[i] = -static_cast<double>(b_N[i + 6 * b_i]) * 1.0E-5;
            }
            robotics::core::internal::SEHelpers::expSE3hat(N, Tm);
            for (i = 0; i < 4; i++) {
                d = Sio[i];
                d1 = Sio[i + 4];
                d2 = Sio[i + 8];
                d3 = Sio[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    Tjop[i + coffset] = ((d * Tio[coffset] + d1 * Tio[coffset + 1]) +
                                         d2 * Tio[coffset + 2]) +
                                        d3 * Tio[coffset + 3];
                }
            }
            for (i = 0; i < 4; i++) {
                d = Tji[i];
                d1 = Tji[i + 4];
                d2 = Tji[i + 8];
                d3 = Tji[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    g_R[i + coffset] = ((d * Tjop[coffset] + d1 * Tjop[coffset + 1]) +
                                        d2 * Tjop[coffset + 2]) +
                                       d3 * Tjop[coffset + 3];
                }
            }
            for (i = 0; i < 4; i++) {
                d = g_R[i];
                d1 = g_R[i + 4];
                d2 = g_R[i + 8];
                d3 = g_R[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    Tjop[i + coffset] = ((d * Sji[coffset] + d1 * Sji[coffset + 1]) +
                                         d2 * Sji[coffset + 2]) +
                                        d3 * Sji[coffset + 3];
                }
            }
            robotics::core::internal::SEHelpers::veelogmSE3(Tjop, dv1);
            for (i = 0; i < 4; i++) {
                d = Tm[i];
                d1 = Tm[i + 4];
                d2 = Tm[i + 8];
                d3 = Tm[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    Tjop[i + coffset] = ((d * Tio[coffset] + d1 * Tio[coffset + 1]) +
                                         d2 * Tio[coffset + 2]) +
                                        d3 * Tio[coffset + 3];
                }
            }
            for (i = 0; i < 4; i++) {
                d = Tji[i];
                d1 = Tji[i + 4];
                d2 = Tji[i + 8];
                d3 = Tji[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    g_R[i + coffset] = ((d * Tjop[coffset] + d1 * Tjop[coffset + 1]) +
                                        d2 * Tjop[coffset + 2]) +
                                       d3 * Tjop[coffset + 3];
                }
            }
            for (i = 0; i < 4; i++) {
                d = g_R[i];
                d1 = g_R[i + 4];
                d2 = g_R[i + 8];
                d3 = g_R[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    Tjop[i + coffset] = ((d * Sji[coffset] + d1 * Sji[coffset + 1]) +
                                         d2 * Sji[coffset + 2]) +
                                        d3 * Sji[coffset + 3];
                }
            }
            robotics::core::internal::SEHelpers::veelogmSE3(Tjop, N);
            for (i = 0; i < 6; i++) {
                Jaci[i + 6 * b_i] = (dv1[i] - N[i]) / 2.0E-5;
            }
            for (i = 0; i < 4; i++) {
                d = Sio[i];
                d1 = Sio[i + 4];
                d2 = Sio[i + 8];
                d3 = Sio[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    Tjop[i + coffset] = ((d * Tjo[coffset] + d1 * Tjo[coffset + 1]) +
                                         d2 * Tjo[coffset + 2]) +
                                        d3 * Tjo[coffset + 3];
                }
                d = Tm[i];
                d1 = Tm[i + 4];
                d2 = Tm[i + 8];
                d3 = Tm[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    Sio[i + coffset] = ((d * Tjo[coffset] + d1 * Tjo[coffset + 1]) +
                                        d2 * Tjo[coffset + 2]) +
                                       d3 * Tjo[coffset + 3];
                }
            }
            for (i = 0; i < 3; i++) {
                R[3 * i] = Tjop[i];
                b_R[3 * i] = Sio[i];
                boffset = 3 * i + 1;
                R[boffset] = Tjop[i + 4];
                b_R[boffset] = Sio[i + 4];
                boffset = 3 * i + 2;
                R[boffset] = Tjop[i + 8];
                b_R[boffset] = Sio[i + 8];
            }
            for (i = 0; i < 9; i++) {
                c_R[i] = -R[i];
            }
            d = Tjop[12];
            d1 = Tjop[13];
            d2 = Tjop[14];
            for (i = 0; i < 3; i++) {
                boffset = i << 2;
                g_R[boffset] = R[3 * i];
                g_R[boffset + 1] = R[3 * i + 1];
                g_R[boffset + 2] = R[3 * i + 2];
                g_R[i + 12] = (c_R[i] * d + c_R[i + 3] * d1) + c_R[i + 6] * d2;
            }
            g_R[3] = 0.0;
            g_R[7] = 0.0;
            g_R[11] = 0.0;
            g_R[15] = 1.0;
            for (i = 0; i < 4; i++) {
                d = deltatform[i];
                d1 = deltatform[i + 4];
                d2 = deltatform[i + 8];
                d3 = deltatform[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    Tjop[i + coffset] = ((d * g_R[coffset] + d1 * g_R[coffset + 1]) +
                                         d2 * g_R[coffset + 2]) +
                                        d3 * g_R[coffset + 3];
                }
            }
            robotics::core::internal::SEHelpers::veelogmSE3(Tjop, dv1);
            for (i = 0; i < 9; i++) {
                R[i] = -b_R[i];
            }
            d = Sio[12];
            d1 = Sio[13];
            d2 = Sio[14];
            for (i = 0; i < 3; i++) {
                boffset = i << 2;
                g_R[boffset] = b_R[3 * i];
                g_R[boffset + 1] = b_R[3 * i + 1];
                g_R[boffset + 2] = b_R[3 * i + 2];
                g_R[i + 12] = (R[i] * d + R[i + 3] * d1) + R[i + 6] * d2;
            }
            g_R[3] = 0.0;
            g_R[7] = 0.0;
            g_R[11] = 0.0;
            g_R[15] = 1.0;
            for (i = 0; i < 4; i++) {
                d = deltatform[i];
                d1 = deltatform[i + 4];
                d2 = deltatform[i + 8];
                d3 = deltatform[i + 12];
                for (j = 0; j < 4; j++) {
                    coffset = j << 2;
                    Tjop[i + coffset] = ((d * g_R[coffset] + d1 * g_R[coffset + 1]) +
                                         d2 * g_R[coffset + 2]) +
                                        d3 * g_R[coffset + 3];
                }
            }
            robotics::core::internal::SEHelpers::veelogmSE3(Tjop, N);
            for (i = 0; i < 6; i++) {
                Jacj[i + 6 * b_i] = (dv1[i] - N[i]) / 2.0E-5;
            }
        }
        for (i = 0; i < 4; i++) {
            d = deltatform[i];
            d1 = deltatform[i + 4];
            d2 = deltatform[i + 8];
            d3 = deltatform[i + 12];
            for (j = 0; j < 4; j++) {
                coffset = j << 2;
                Tjop[i + coffset] = ((d * Sji[coffset] + d1 * Sji[coffset + 1]) +
                                     d2 * Sji[coffset + 2]) +
                                    d3 * Sji[coffset + 3];
            }
        }
        robotics::core::internal::SEHelpers::veelogmSE3(Tjop, dv1);
        b_i = 6;
        for (i = 0; i < 6; i++) {
            e_data[i] = dv1[i];
        }
        Jaci_size[0] = 6;
        Jaci_size[1] = 6;
        Jacj_size[0] = 6;
        Jacj_size[1] = 6;
        std::copy(&Jaci[0], &Jaci[36], &Jaci_data[0]);
        std::copy(&Jacj[0], &Jacj[36], &Jacj_data[0]);
    } else {
        double Sio[16];
        double Sji[16];
        double Tji[16];
        double Tjop[16];
        double g_R[16];
        double R[9];
        double b_R[9];
        double deltavec[7];
        double b_measurement;
        double cosTheta;
        double d;
        double d1;
        double d2;
        double q;
        double sinThetaDy;
        for (i = 0; i < 4; i++) {
            coffset = i << 2;
            Tjop[coffset] = Toi[Toi.size(0) * i];
            Tji[coffset] = Toj[Toj.size(0) * i];
            Sio[coffset] = measurement[measurement.size(0) * i];
            Tjop[coffset + 1] = Toi[Toi.size(0) * i + 1];
            Tji[coffset + 1] = Toj[Toj.size(0) * i + 1];
            Sio[coffset + 1] = measurement[measurement.size(0) * i + 1];
            Tjop[coffset + 2] = Toi[Toi.size(0) * i + 2];
            Tji[coffset + 2] = Toj[Toj.size(0) * i + 2];
            Sio[coffset + 2] = measurement[measurement.size(0) * i + 2];
            Tjop[coffset + 3] = Toi[Toi.size(0) * i + 3];
            Tji[coffset + 3] = Toj[Toj.size(0) * i + 3];
            Sio[coffset + 3] = measurement[measurement.size(0) * i + 3];
        }
        for (i = 0; i < 3; i++) {
            R[3 * i] = Tjop[i];
            b_R[3 * i] = Sio[i];
            boffset = 3 * i + 1;
            R[boffset] = Tjop[i + 4];
            b_R[boffset] = Sio[i + 4];
            boffset = 3 * i + 2;
            R[boffset] = Tjop[i + 8];
            b_R[boffset] = Sio[i + 8];
        }
        b_measurement = measurement[measurement.size(0) * 3 + 3];
        d = Sio[12];
        d1 = Sio[13];
        d2 = Sio[14];
        for (i = 0; i < 3; i++) {
            coffset = i << 2;
            Sji[coffset] = b_R[3 * i];
            Sji[coffset + 1] = b_R[3 * i + 1];
            Sji[coffset + 2] = b_R[3 * i + 2];
            Sji[i + 12] =
                -((b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2) / b_measurement;
        }
        Sji[3] = 0.0;
        Sji[7] = 0.0;
        Sji[11] = 0.0;
        Sji[15] = 1.0 / measurement[measurement.size(0) * 3 + 3];
        for (i = 0; i < 3; i++) {
            b_R[3 * i] = Tjop[i];
            b_R[3 * i + 1] = Tjop[i + 4];
            b_R[3 * i + 2] = Tjop[i + 8];
        }
        q = Toi[Toi.size(0) * 3 + 3];
        d = Tjop[12];
        d1 = Tjop[13];
        d2 = Tjop[14];
        for (i = 0; i < 3; i++) {
            coffset = i << 2;
            Sio[coffset] = b_R[3 * i];
            Sio[coffset + 1] = b_R[3 * i + 1];
            Sio[coffset + 2] = b_R[3 * i + 2];
            Sio[i + 12] = -((b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2) / q;
        }
        Sio[3] = 0.0;
        Sio[7] = 0.0;
        Sio[11] = 0.0;
        Sio[15] = 1.0 / Toi[Toi.size(0) * 3 + 3];
        q = Toi[Toi.size(0) * 3 + 3];
        b_measurement = Toj[Toj.size(0) * 3 + 3];
        cosTheta = Toi[Toi.size(0) * 3 + 3];
        sinThetaDy = Toj[Toj.size(0) * 3 + 3];
        d = Toi[Toi.size(0) * 3 + 3];
        d1 = Toj[Toj.size(0) * 3 + 3];
        for (int k{0}; k < 7; k++) {
            double Tjo[16];
            double Tm[16];
            double deltatform[16];
            double c_R[9];
            double t[3];
            double d3;
            double d4;
            double s;
            double sinTheta;
            for (b_i = 0; b_i < 7; b_i++) {
                deltavec[b_i] = 0.0;
            }
            deltavec[k] = 1.0E-9;
            robotics::core::internal::Sim3Helpers::sim3ToSform(deltavec, Tjo);
            for (i = 0; i < 3; i++) {
                b_R[3 * i] = Tjo[i];
                b_R[3 * i + 1] = Tjo[i + 4];
                b_R[3 * i + 2] = Tjo[i + 8];
            }
            d2 = Tjo[12];
            d3 = Tjo[13];
            d4 = Tjo[14];
            sinTheta = Tjo[15];
            for (i = 0; i < 3; i++) {
                coffset = i << 2;
                deltatform[coffset] = b_R[3 * i];
                deltatform[coffset + 1] = b_R[3 * i + 1];
                deltatform[coffset + 2] = b_R[3 * i + 2];
                deltatform[i + 12] =
                    -((b_R[i] * d2 + b_R[i + 3] * d3) + b_R[i + 6] * d4) / sinTheta;
            }
            deltatform[3] = 0.0;
            deltatform[7] = 0.0;
            deltatform[11] = 0.0;
            deltatform[15] = 1.0 / Tjo[15];
            for (i = 0; i < 3; i++) {
                d2 = 0.0;
                for (j = 0; j < 3; j++) {
                    coffset = j << 2;
                    f_R[i + 3 * j] = (Tjop[i] * deltatform[coffset] +
                                      Tjop[i + 4] * deltatform[coffset + 1]) +
                                     Tjop[i + 8] * deltatform[coffset + 2];
                    d2 += q * Tjop[i + coffset] * deltatform[j + 12];
                }
                e_R[i] = d2 + Tjop[i + 12];
            }
            for (i = 0; i < 3; i++) {
                coffset = i << 2;
                Tjo[coffset] = f_R[3 * i];
                Tjo[coffset + 1] = f_R[3 * i + 1];
                Tjo[coffset + 2] = f_R[3 * i + 2];
                Tjo[i + 12] = e_R[i];
            }
            Tjo[3] = 0.0;
            Tjo[7] = 0.0;
            Tjo[11] = 0.0;
            Tjo[15] = d * deltatform[15];
            s = d1 * deltatform[15];
            deltavec[k] = -1.0E-9;
            robotics::core::internal::Sim3Helpers::sim3ToSform(deltavec, Tm);
            for (i = 0; i < 3; i++) {
                d2 = 0.0;
                for (j = 0; j < 3; j++) {
                    boffset = j << 2;
                    coffset = i + boffset;
                    aoffset = j + 3 * i;
                    b_R[aoffset] = Tjo[coffset];
                    c_R[i + 3 * j] = (Tji[i] * deltatform[boffset] +
                                      Tji[i + 4] * deltatform[boffset + 1]) +
                                     Tji[i + 8] * deltatform[boffset + 2];
                    d2 += b_measurement * Tji[coffset] * deltatform[j + 12];
                    f_R[aoffset] = Tm[coffset];
                }
                t[i] = d2 + Tji[i + 12];
            }
            d2 = Tm[12];
            d3 = Tm[13];
            d4 = Tm[14];
            sinTheta = Tm[15];
            for (i = 0; i < 3; i++) {
                coffset = i << 2;
                deltatform[coffset] = f_R[3 * i];
                deltatform[coffset + 1] = f_R[3 * i + 1];
                deltatform[coffset + 2] = f_R[3 * i + 2];
                deltatform[i + 12] =
                    -((f_R[i] * d2 + f_R[i + 3] * d3) + f_R[i + 6] * d4) / sinTheta;
            }
            deltatform[3] = 0.0;
            deltatform[7] = 0.0;
            deltatform[11] = 0.0;
            deltatform[15] = 1.0 / Tm[15];
            for (i = 0; i < 3; i++) {
                d2 = 0.0;
                for (j = 0; j < 3; j++) {
                    coffset = j << 2;
                    f_R[i + 3 * j] = (Tjop[i] * deltatform[coffset] +
                                      Tjop[i + 4] * deltatform[coffset + 1]) +
                                     Tjop[i + 8] * deltatform[coffset + 2];
                    d2 += cosTheta * Tjop[i + coffset] * deltatform[j + 12];
                }
                e_R[i] = d2 + Tjop[i + 12];
            }
            for (i = 0; i < 3; i++) {
                coffset = i << 2;
                Tm[coffset] = f_R[3 * i];
                Tm[coffset + 1] = f_R[3 * i + 1];
                Tm[coffset + 2] = f_R[3 * i + 2];
                Tm[i + 12] = e_R[i];
            }
            Tm[3] = 0.0;
            Tm[7] = 0.0;
            Tm[11] = 0.0;
            Tm[15] = d * deltatform[15];
            d2 = Tjo[12];
            d3 = Tjo[13];
            d4 = Tjo[14];
            sinTheta = Tjo[15];
            for (i = 0; i < 3; i++) {
                f_R[3 * i] = Tm[i];
                boffset = i << 2;
                g_R[boffset] = b_R[3 * i];
                coffset = 3 * i + 1;
                f_R[coffset] = Tm[i + 4];
                g_R[boffset + 1] = b_R[coffset];
                coffset = 3 * i + 2;
                f_R[coffset] = Tm[i + 8];
                g_R[boffset + 2] = b_R[coffset];
                g_R[i + 12] =
                    -((b_R[i] * d2 + b_R[i + 3] * d3) + b_R[i + 6] * d4) / sinTheta;
            }
            g_R[3] = 0.0;
            g_R[7] = 0.0;
            g_R[11] = 0.0;
            g_R[15] = 1.0 / Tjo[15];
            robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, g_R, Tji,
                                                                   deltavec);
            d2 = Tm[12];
            d3 = Tm[13];
            d4 = Tm[14];
            sinTheta = Tm[15];
            for (i = 0; i < 3; i++) {
                boffset = i << 2;
                g_R[boffset] = f_R[3 * i];
                g_R[boffset + 1] = f_R[3 * i + 1];
                g_R[boffset + 2] = f_R[3 * i + 2];
                g_R[i + 12] =
                    -((f_R[i] * d2 + f_R[i + 3] * d3) + f_R[i + 6] * d4) / sinTheta;
            }
            g_R[3] = 0.0;
            g_R[7] = 0.0;
            g_R[11] = 0.0;
            g_R[15] = 1.0 / Tm[15];
            robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, g_R, Tji,
                                                                   e_data);
            for (i = 0; i < 7; i++) {
                Jaci_data[i + 7 * k] =
                    (deltavec[i] - e_data[i]) * 4.9999999999999994E+8;
            }
            for (i = 0; i < 3; i++) {
                boffset = i << 2;
                g_R[boffset] = c_R[3 * i];
                g_R[boffset + 1] = c_R[3 * i + 1];
                g_R[boffset + 2] = c_R[3 * i + 2];
                g_R[i + 12] = t[i];
            }
            g_R[3] = 0.0;
            g_R[7] = 0.0;
            g_R[11] = 0.0;
            g_R[15] = s;
            robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, Sio, g_R,
                                                                   deltavec);
            for (i = 0; i < 3; i++) {
                d2 = 0.0;
                for (j = 0; j < 3; j++) {
                    coffset = j << 2;
                    f_R[i + 3 * j] = (Tji[i] * deltatform[coffset] +
                                      Tji[i + 4] * deltatform[coffset + 1]) +
                                     Tji[i + 8] * deltatform[coffset + 2];
                    d2 += sinThetaDy * Tji[i + coffset] * deltatform[j + 12];
                }
                e_R[i] = d2 + Tji[i + 12];
            }
            for (i = 0; i < 3; i++) {
                boffset = i << 2;
                g_R[boffset] = f_R[3 * i];
                g_R[boffset + 1] = f_R[3 * i + 1];
                g_R[boffset + 2] = f_R[3 * i + 2];
                g_R[i + 12] = e_R[i];
            }
            g_R[3] = 0.0;
            g_R[7] = 0.0;
            g_R[11] = 0.0;
            g_R[15] = d1 * deltatform[15];
            robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, Sio, g_R,
                                                                   e_data);
            for (i = 0; i < 7; i++) {
                Jacj_data[i + 7 * k] =
                    (deltavec[i] - e_data[i]) * 4.9999999999999994E+8;
            }
        }
        q = Toi[Toi.size(0) * 3 + 3];
        d = Tjop[12];
        d1 = Tjop[13];
        d2 = Tjop[14];
        for (i = 0; i < 3; i++) {
            boffset = i << 2;
            g_R[boffset] = R[3 * i];
            g_R[boffset + 1] = R[3 * i + 1];
            g_R[boffset + 2] = R[3 * i + 2];
            g_R[i + 12] = -((R[i] * d + R[i + 3] * d1) + R[i + 6] * d2) / q;
        }
        g_R[3] = 0.0;
        g_R[7] = 0.0;
        g_R[11] = 0.0;
        g_R[15] = 1.0 / Toi[Toi.size(0) * 3 + 3];
        robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, g_R, Tji,
                                                               deltavec);
        b_i = 7;
        for (i = 0; i < 7; i++) {
            e_data[i] = deltavec[i];
        }
        Jaci_size[0] = 7;
        Jaci_size[1] = 7;
        Jacj_size[0] = 7;
        Jacj_size[1] = 7;
    }
    coffset = Omega.size(1);
    aoffset = Omega.size(1);
    for (j = 0; j < coffset; j++) {
        boffset = j * Omega.size(0);
        e_R[j] = 0.0;
        for (int k{0}; k < b_i; k++) {
            e_R[j] += e_data[k] * Omega[boffset + k];
        }
    }
    cost = 0.0;
    for (i = 0; i < aoffset; i++) {
        cost += e_R[i] * e_data[i];
    }
    ::SlamGraph2D::coder::internal::blas::mtimes(Jaci_data, Jaci_size, Omega,
                                                 y_tmp_data, y_tmp_size);
    coffset = y_tmp_size[0] - 1;
    i = y_tmp_size[1];
    gradi_size = y_tmp_size[0];
    std::memset(&gradi_data[0], 0,
                static_cast<unsigned int>(coffset + 1) * sizeof(double));
    for (int k{0}; k < i; k++) {
        aoffset = k * y_tmp_size[0];
        for (b_i = 0; b_i <= coffset; b_i++) {
            gradi_data[b_i] += y_tmp_data[aoffset + b_i] * e_data[k];
        }
    }
    ::SlamGraph2D::coder::internal::blas::mtimes(Jacj_data, Jacj_size, Omega,
                                                 b_y_tmp_data, b_y_tmp_size);
    coffset = b_y_tmp_size[0] - 1;
    i = b_y_tmp_size[1];
    gradj_size = b_y_tmp_size[0];
    std::memset(&gradj_data[0], 0,
                static_cast<unsigned int>(coffset + 1) * sizeof(double));
    for (int k{0}; k < i; k++) {
        aoffset = k * b_y_tmp_size[0];
        for (b_i = 0; b_i <= coffset; b_i++) {
            gradj_data[b_i] += b_y_tmp_data[aoffset + b_i] * e_data[k];
        }
    }
    ::SlamGraph2D::coder::internal::blas::mtimes(
        y_tmp_data, y_tmp_size, Jaci_data, Jaci_size, hessii_data, hessii_size);
    ::SlamGraph2D::coder::internal::blas::mtimes(
        y_tmp_data, y_tmp_size, Jacj_data, Jacj_size, hessij_data, hessij_size);
    ::SlamGraph2D::coder::internal::blas::mtimes(b_y_tmp_data, b_y_tmp_size,
                                                 Jaci_data, Jaci_size,
                                                 hessji_data, hessji_size);
    ::SlamGraph2D::coder::internal::blas::mtimes(b_y_tmp_data, b_y_tmp_size,
                                                 Jacj_data, Jacj_size,
                                                 hessjj_data, hessjj_size);
    return cost;
}

///
/// @fn             : div_s32
/// @brief          :
/// @param          : int numerator
///                   int denominator
/// @return         : int
///
}  // namespace internal
}  // namespace algs
}  // namespace nav
}  // namespace coder
static int div_s32(int numerator, int denominator) {
    int quotient;
    if (denominator == 0) {
        if (numerator >= 0) {
            quotient = MAX_int32_T;
        } else {
            quotient = MIN_int32_T;
        }
    } else {
        unsigned int tempAbsQuotient;
        unsigned int u;
        if (numerator < 0) {
            tempAbsQuotient = ~static_cast<unsigned int>(numerator) + 1U;
        } else {
            tempAbsQuotient = static_cast<unsigned int>(numerator);
        }
        if (denominator < 0) {
            u = ~static_cast<unsigned int>(denominator) + 1U;
        } else {
            u = static_cast<unsigned int>(denominator);
        }
        tempAbsQuotient /= u;
        if ((numerator < 0) != (denominator < 0)) {
            quotient = -static_cast<int>(tempAbsQuotient);
        } else {
            quotient = static_cast<int>(tempAbsQuotient);
        }
    }
    return quotient;
}

///
/// @fn             : poseGraphCost
/// @brief          :
/// @param          : const ::coder::array<double, 2U> &posesMat
///                   const ::coder::array<double, 2U> &args_edgeNodePairs
///                   const ::coder::array<double, 2U> &args_edgeMeasurements
///                   const ::coder::array<double, 2U> &args_edgeInfoMats
///                   const double args_tformSize[2]
///                   const double args_infoMatSize[2]
///                   double args_poseDeltaLength
///                   const ::coder::array<double, 1U> &args_nodeMap
///                   const ::coder::array<double, 1U> &args_nodeDims
///                   const ::coder::array<bool, 1U> &args_IsLandmarkNode
///                   ::coder::array<double, 1U> &gradient
///                   sparse &hessian
/// @return         : double
///
namespace coder {
namespace nav {
namespace algs {
namespace internal {
double PoseGraphHelpers::poseGraphCost(
    const ::coder::array<double, 2U> &posesMat,
    const ::coder::array<double, 2U> &args_edgeNodePairs,
    const ::coder::array<double, 2U> &args_edgeMeasurements,
    const ::coder::array<double, 2U> &args_edgeInfoMats,
    const double args_tformSize[2], const double args_infoMatSize[2],
    double args_poseDeltaLength, const ::coder::array<double, 1U> &args_nodeMap,
    const ::coder::array<double, 1U> &args_nodeDims,
    const ::coder::array<bool, 1U> &args_IsLandmarkNode,
    ::coder::array<double, 1U> &gradient, sparse &hessian) {
    BlockInserter2 bi;
    robotics::core::internal::BlockMatrix edgeInfoMats;
    robotics::core::internal::BlockMatrix edgeMeasurements;
    robotics::core::internal::BlockMatrix poses;
    ::coder::array<double, 2U> OmegaIn;
    ::coder::array<double, 2U> Tij;
    ::coder::array<double, 2U> Toi;
    ::coder::array<double, 2U> Toj;
    ::coder::array<double, 2U> varargin_1;
    ::coder::array<double, 2U> y;
    ::coder::array<double, 1U> b_bi;
    ::coder::array<double, 1U> c_bi;
    ::coder::array<double, 1U> d_bi;
    ::coder::array<int, 2U> b_m1;
    ::coder::array<int, 2U> m2;
    ::coder::array<int, 2U> m3;
    ::coder::array<int, 2U> m4;
    ::coder::array<int, 2U> n1;
    ::coder::array<int, 2U> n2;
    ::coder::array<int, 2U> n3;
    ::coder::array<int, 2U> v1;
    ::coder::array<int, 2U> vk;
    ::coder::array<signed char, 2U> b_I;
    double H_data[196];
    double hessii_data[49];
    double hessij_data[49];
    double hessji_data[49];
    double hessjj_data[49];
    double gradi_data[7];
    double gradj_data[7];
    double cost;
    double maxNodeDim;
    double numEntries1;
    int i;
    int i1;
    int i2;
    int idx;
    int k;
    int last;
    int loop_ub;
    edgeMeasurements.Matrix.set_size(args_edgeMeasurements.size(0), 3);
    loop_ub = args_edgeMeasurements.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        edgeMeasurements.Matrix[i] = args_edgeMeasurements[i];
    }
    edgeMeasurements.BlockSize[0] = args_tformSize[0];
    edgeMeasurements.BlockSize[1] = args_tformSize[1];
    edgeMeasurements.NumRowBlocks =
        static_cast<double>(args_edgeMeasurements.size(0)) / args_tformSize[0];
    edgeMeasurements.NumColBlocks = 3.0 / args_tformSize[1];
    edgeInfoMats.Matrix.set_size(args_edgeInfoMats.size(0), 3);
    loop_ub = args_edgeInfoMats.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        edgeInfoMats.Matrix[i] = args_edgeInfoMats[i];
    }
    edgeInfoMats.BlockSize[0] = args_infoMatSize[0];
    edgeInfoMats.BlockSize[1] = args_infoMatSize[1];
    edgeInfoMats.NumRowBlocks =
        static_cast<double>(args_edgeInfoMats.size(0)) / args_infoMatSize[0];
    edgeInfoMats.NumColBlocks = 3.0 / args_infoMatSize[1];
    cost = 0.0;
    poses.Matrix.set_size(posesMat.size(0), 3);
    loop_ub = posesMat.size(0) * 3;
    for (i = 0; i < loop_ub; i++) {
        poses.Matrix[i] = posesMat[i];
    }
    poses.BlockSize[0] = args_tformSize[0];
    poses.BlockSize[1] = args_tformSize[1];
    poses.NumRowBlocks =
        static_cast<double>(posesMat.size(0)) / args_tformSize[0];
    poses.NumColBlocks = 3.0 / args_tformSize[1];
    last = static_cast<int>(poses.NumRowBlocks * args_poseDeltaLength);
    bi.Gradient.set_size(last);
    for (i = 0; i < last; i++) {
        bi.Gradient[i] = 0.0;
    }
    bi.NodeDims.set_size(args_nodeDims.size(0));
    loop_ub = args_nodeDims.size(0);
    for (i = 0; i < loop_ub; i++) {
        bi.NodeDims[i] = args_nodeDims[i];
    }
    bi.NodeMap.set_size(args_nodeMap.size(0));
    loop_ub = args_nodeMap.size(0);
    for (i = 0; i < loop_ub; i++) {
        bi.NodeMap[i] = args_nodeMap[i];
    }
    last = args_nodeDims.size(0);
    if (args_nodeDims.size(0) <= 2) {
        if (args_nodeDims.size(0) == 1) {
            maxNodeDim = args_nodeDims[0];
        } else {
            maxNodeDim = args_nodeDims[args_nodeDims.size(0) - 1];
            if ((!(args_nodeDims[0] < maxNodeDim)) &&
                ((!std::isnan(args_nodeDims[0])) || std::isnan(maxNodeDim))) {
                maxNodeDim = args_nodeDims[0];
            }
        }
    } else {
        if (!std::isnan(args_nodeDims[0])) {
            idx = 1;
        } else {
            bool exitg1;
            idx = 0;
            k = 2;
            exitg1 = false;
            while ((!exitg1) && (k <= last)) {
                if (!std::isnan(args_nodeDims[k - 1])) {
                    idx = k;
                    exitg1 = true;
                } else {
                    k++;
                }
            }
        }
        if (idx == 0) {
            maxNodeDim = args_nodeDims[0];
        } else {
            maxNodeDim = args_nodeDims[idx - 1];
            i = idx + 1;
            for (k = i; k <= last; k++) {
                numEntries1 = args_nodeDims[k - 1];
                if (maxNodeDim < numEntries1) {
                    maxNodeDim = numEntries1;
                }
            }
        }
    }
    i = static_cast<int>(
        (4.0 * static_cast<double>(args_edgeNodePairs.size(0)) + 1.0) *
        maxNodeDim * maxNodeDim);
    bi.HessianCSC.set_size(i, 3);
    loop_ub = i * 3;
    for (i = 0; i < loop_ub; i++) {
        bi.HessianCSC[i] = 0.0;
    }
    bi.HessianCSCCount = 1.0;
    i = args_edgeNodePairs.size(0);
    for (k = 0; k < i; k++) {
        double b_i;
        double j;
        double numEntries2;
        double numEntries2_tmp;
        int hessii_size[2];
        int hessij_size[2];
        int hessji_size[2];
        int hessjj_size[2];
        int b_loop_ub;
        int c_loop_ub;
        int m1;
        int m2_idx_0;
        int m3_idx_0;
        int varargin_1_tmp;
        signed char i3;
        bool b;
        b_i = args_edgeNodePairs[k];
        j = args_edgeNodePairs[k + args_edgeNodePairs.size(0)];
        edgeMeasurements.extractBlock(static_cast<double>(k) + 1.0, Tij);
        edgeInfoMats.extractBlock(static_cast<double>(k) + 1.0, OmegaIn);
        b = args_IsLandmarkNode[static_cast<int>(j) - 1];
        if (b) {
            maxNodeDim = args_nodeDims[static_cast<int>(j) - 1];
            if (maxNodeDim < 1.0) {
                last = 0;
            } else {
                last = static_cast<int>(maxNodeDim);
            }
            for (i1 = 0; i1 < last; i1++) {
                for (i2 = 0; i2 < last; i2++) {
                    OmegaIn[i2 + last * i1] = OmegaIn[i2 + OmegaIn.size(0) * i1];
                }
            }
            OmegaIn.set_size(last, last);
        }
        poses.extractBlock(b_i, Toi);
        poses.extractBlock(j, Toj);
        maxNodeDim = PoseGraphHelpers::costBetweenTwoNodes(
            Toi, Toj, Tij, OmegaIn, b, gradi_data, idx, gradj_data, last,
            hessii_data, hessii_size, hessij_data, hessij_size, hessji_data,
            hessji_size, hessjj_data, hessjj_size);
        cost += maxNodeDim;
        bi.insertGradientBlock(b_i, gradi_data, idx);
        bi.insertGradientBlock(j, gradj_data, last);
        maxNodeDim = bi.NodeDims[static_cast<int>(b_i) - 1];
        numEntries1 = maxNodeDim * maxNodeDim;
        if (std::isnan(numEntries1)) {
            y.set_size(1, 1);
            y[0] = rtNaN;
        } else if (numEntries1 < 1.0) {
            y.set_size(1, 0);
        } else {
            y.set_size(1, static_cast<int>(numEntries1 - 1.0) + 1);
            loop_ub = static_cast<int>(numEntries1 - 1.0);
            for (i1 = 0; i1 <= loop_ub; i1++) {
                y[i1] = static_cast<double>(i1) + 1.0;
            }
        }
        v1.set_size(1, y.size(1));
        loop_ub = y.size(1);
        for (i1 = 0; i1 < loop_ub; i1++) {
            v1[i1] = static_cast<int>(y[i1]) - 1;
        }
        vk.set_size(1, v1.size(1));
        loop_ub = v1.size(1);
        for (i1 = 0; i1 < loop_ub; i1++) {
            vk[i1] = div_s32(v1[i1], static_cast<int>(maxNodeDim));
        }
        v1.set_size(1, v1.size(1));
        loop_ub = v1.size(1) - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            v1[i1] = v1[i1] - vk[i1] * static_cast<int>(maxNodeDim);
        }
        b_m1.set_size(1, v1.size(1));
        loop_ub = v1.size(1);
        n1.set_size(1, vk.size(1));
        for (i1 = 0; i1 < loop_ub; i1++) {
            b_m1[i1] = v1[i1] + 1;
            n1[i1] = vk[i1] + 1;
        }
        numEntries2_tmp = bi.NodeDims[static_cast<int>(j) - 1];
        numEntries2 = maxNodeDim * numEntries2_tmp;
        b = std::isnan(numEntries2);
        if (b) {
            y.set_size(1, 1);
            y[0] = rtNaN;
        } else if (numEntries2 < 1.0) {
            y.set_size(1, 0);
        } else {
            y.set_size(1, static_cast<int>(numEntries2 - 1.0) + 1);
            loop_ub = static_cast<int>(numEntries2 - 1.0);
            for (i1 = 0; i1 <= loop_ub; i1++) {
                y[i1] = static_cast<double>(i1) + 1.0;
            }
        }
        v1.set_size(1, y.size(1));
        loop_ub = y.size(1);
        for (i1 = 0; i1 < loop_ub; i1++) {
            v1[i1] = static_cast<int>(y[i1]) - 1;
        }
        vk.set_size(1, v1.size(1));
        loop_ub = v1.size(1);
        for (i1 = 0; i1 < loop_ub; i1++) {
            vk[i1] = div_s32(v1[i1], static_cast<int>(maxNodeDim));
        }
        v1.set_size(1, v1.size(1));
        loop_ub = v1.size(1) - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            v1[i1] = v1[i1] - vk[i1] * static_cast<int>(maxNodeDim);
        }
        m2.set_size(1, v1.size(1));
        loop_ub = v1.size(1);
        n2.set_size(1, vk.size(1));
        for (i1 = 0; i1 < loop_ub; i1++) {
            m2[i1] = v1[i1] + 1;
            n2[i1] = vk[i1] + 1;
        }
        if (b) {
            y.set_size(1, 1);
            y[0] = rtNaN;
        } else if (numEntries2 < 1.0) {
            y.set_size(1, 0);
        } else {
            y.set_size(1, static_cast<int>(numEntries2 - 1.0) + 1);
            loop_ub = static_cast<int>(numEntries2 - 1.0);
            for (i1 = 0; i1 <= loop_ub; i1++) {
                y[i1] = static_cast<double>(i1) + 1.0;
            }
        }
        v1.set_size(1, y.size(1));
        loop_ub = y.size(1);
        for (i1 = 0; i1 < loop_ub; i1++) {
            v1[i1] = static_cast<int>(y[i1]) - 1;
        }
        vk.set_size(1, v1.size(1));
        loop_ub = v1.size(1);
        for (i1 = 0; i1 < loop_ub; i1++) {
            vk[i1] = div_s32(v1[i1], static_cast<int>(numEntries2_tmp));
        }
        v1.set_size(1, v1.size(1));
        loop_ub = v1.size(1) - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            v1[i1] = v1[i1] - vk[i1] * static_cast<int>(numEntries2_tmp);
        }
        m3.set_size(1, v1.size(1));
        loop_ub = v1.size(1);
        n3.set_size(1, vk.size(1));
        for (i1 = 0; i1 < loop_ub; i1++) {
            m3[i1] = v1[i1] + 1;
            n3[i1] = vk[i1] + 1;
        }
        maxNodeDim = numEntries2_tmp * numEntries2_tmp;
        if (std::isnan(maxNodeDim)) {
            y.set_size(1, 1);
            y[0] = rtNaN;
        } else if (maxNodeDim < 1.0) {
            y.set_size(1, 0);
        } else {
            y.set_size(1, static_cast<int>(maxNodeDim - 1.0) + 1);
            loop_ub = static_cast<int>(maxNodeDim - 1.0);
            for (i1 = 0; i1 <= loop_ub; i1++) {
                y[i1] = static_cast<double>(i1) + 1.0;
            }
        }
        v1.set_size(1, y.size(1));
        loop_ub = y.size(1);
        for (i1 = 0; i1 < loop_ub; i1++) {
            v1[i1] = static_cast<int>(y[i1]) - 1;
        }
        vk.set_size(1, v1.size(1));
        loop_ub = v1.size(1);
        for (i1 = 0; i1 < loop_ub; i1++) {
            vk[i1] = div_s32(v1[i1], static_cast<int>(numEntries2_tmp));
        }
        v1.set_size(1, v1.size(1));
        loop_ub = v1.size(1) - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            v1[i1] = v1[i1] - vk[i1] * static_cast<int>(numEntries2_tmp);
        }
        m4.set_size(1, v1.size(1));
        loop_ub = v1.size(1);
        y.set_size(1, vk.size(1));
        for (i1 = 0; i1 < loop_ub; i1++) {
            m4[i1] = v1[i1] + 1;
            y[i1] = vk[i1] + 1;
        }
        loop_ub = hessii_size[0] * hessii_size[1];
        idx = hessij_size[0] * hessij_size[1];
        last = hessji_size[0] * hessji_size[1];
        b_loop_ub = hessjj_size[0] * hessjj_size[1];
        if (loop_ub - 1 >= 0) {
            std::copy(&hessii_data[0], &hessii_data[loop_ub], &H_data[0]);
        }
        for (i1 = 0; i1 < idx; i1++) {
            H_data[i1 + loop_ub] = hessij_data[i1];
        }
        for (i1 = 0; i1 < last; i1++) {
            H_data[(i1 + loop_ub) + idx] = hessji_data[i1];
        }
        for (i1 = 0; i1 < b_loop_ub; i1++) {
            H_data[((i1 + loop_ub) + idx) + last] = hessjj_data[i1];
        }
        numEntries1 =
            bi.HessianCSCCount + ((numEntries1 + 2.0 * numEntries2) + maxNodeDim);
        if (bi.HessianCSCCount > numEntries1 - 1.0) {
            i1 = 0;
        } else {
            i1 = static_cast<int>(bi.HessianCSCCount) - 1;
        }
        m1 = b_m1.size(1);
        m2_idx_0 = m2.size(1);
        m3_idx_0 = m3.size(1);
        varargin_1.set_size(((b_m1.size(1) + m2.size(1)) + m3.size(1)) + m4.size(1),
                            2);
        c_loop_ub = b_m1.size(1);
        for (i2 = 0; i2 < c_loop_ub; i2++) {
            maxNodeDim = bi.NodeMap[static_cast<int>(b_i) - 1] - 1.0;
            varargin_1[i2] = maxNodeDim + static_cast<double>(b_m1[i2]);
            varargin_1[i2 + varargin_1.size(0)] =
                maxNodeDim + static_cast<double>(n1[i2]);
        }
        c_loop_ub = m2.size(1);
        for (i2 = 0; i2 < c_loop_ub; i2++) {
            varargin_1_tmp = i2 + m1;
            varargin_1[varargin_1_tmp] =
                (bi.NodeMap[static_cast<int>(b_i) - 1] - 1.0) +
                static_cast<double>(m2[i2]);
            varargin_1[varargin_1_tmp + varargin_1.size(0)] =
                (bi.NodeMap[static_cast<int>(j) - 1] - 1.0) +
                static_cast<double>(n2[i2]);
        }
        c_loop_ub = m3.size(1);
        for (i2 = 0; i2 < c_loop_ub; i2++) {
            varargin_1_tmp = (i2 + m1) + m2_idx_0;
            varargin_1[varargin_1_tmp] = (bi.NodeMap[static_cast<int>(j) - 1] - 1.0) +
                                         static_cast<double>(m3[i2]);
            varargin_1[varargin_1_tmp + varargin_1.size(0)] =
                (bi.NodeMap[static_cast<int>(b_i) - 1] - 1.0) +
                static_cast<double>(n3[i2]);
        }
        c_loop_ub = m4.size(1);
        for (i2 = 0; i2 < c_loop_ub; i2++) {
            maxNodeDim = bi.NodeMap[static_cast<int>(j) - 1] - 1.0;
            varargin_1_tmp = ((i2 + m1) + m2_idx_0) + m3_idx_0;
            varargin_1[varargin_1_tmp] = maxNodeDim + static_cast<double>(m4[i2]);
            varargin_1[varargin_1_tmp + varargin_1.size(0)] = maxNodeDim + y[i2];
        }
        if (varargin_1.size(0) != 0) {
            idx = varargin_1.size(0);
        } else {
            idx = ((loop_ub + idx) + last) + b_loop_ub;
        }
        if ((idx == 0) || (varargin_1.size(0) != 0)) {
            i3 = 2;
        } else {
            i3 = 0;
        }
        loop_ub = i3;
        for (i2 = 0; i2 < loop_ub; i2++) {
            for (last = 0; last < idx; last++) {
                bi.HessianCSC[(i1 + last) + bi.HessianCSC.size(0) * i2] =
                    varargin_1[last + idx * i2];
            }
        }
        for (i2 = 0; i2 < idx; i2++) {
            bi.HessianCSC[(i1 + i2) + bi.HessianCSC.size(0) * i3] = H_data[i2];
        }
        bi.HessianCSCCount = numEntries1;
    }
    if (args_poseDeltaLength < 0.0) {
        maxNodeDim = 0.0;
        idx = 0;
    } else {
        maxNodeDim = args_poseDeltaLength;
        idx = static_cast<int>(args_poseDeltaLength);
    }
    b_I.set_size(static_cast<int>(maxNodeDim), static_cast<int>(maxNodeDim));
    last = static_cast<int>(maxNodeDim) * static_cast<int>(maxNodeDim);
    for (i = 0; i < last; i++) {
        b_I[i] = 0;
    }
    if (static_cast<int>(maxNodeDim) > 0) {
        for (k = 0; k < idx; k++) {
            b_I[k + b_I.size(0) * k] = 1;
        }
    }
    if (last < 1) {
        y.set_size(1, 0);
    } else {
        y.set_size(1, last);
        loop_ub = last - 1;
        for (i = 0; i <= loop_ub; i++) {
            y[i] = static_cast<double>(i) + 1.0;
        }
    }
    idx = b_I.size(0);
    v1.set_size(1, y.size(1));
    loop_ub = y.size(1);
    for (i = 0; i < loop_ub; i++) {
        v1[i] = static_cast<int>(y[i]) - 1;
    }
    vk.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    for (i = 0; i < loop_ub; i++) {
        vk[i] = div_s32(v1[i], idx);
    }
    v1.set_size(1, v1.size(1));
    loop_ub = v1.size(1) - 1;
    for (i = 0; i <= loop_ub; i++) {
        v1[i] = v1[i] - vk[i] * idx;
    }
    y.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    m4.set_size(1, vk.size(1));
    for (i = 0; i < loop_ub; i++) {
        y[i] = v1[i] + 1;
        m4[i] = vk[i] + 1;
    }
    numEntries1 = bi.HessianCSCCount + static_cast<double>(last);
    if (bi.HessianCSCCount > numEntries1 - 1.0) {
        i = 0;
    } else {
        i = static_cast<int>(bi.HessianCSCCount) - 1;
    }
    loop_ub = y.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
        i2 = i + i1;
        bi.HessianCSC[i2] = (bi.NodeMap[0] + y[i1]) - 1.0;
        bi.HessianCSC[i2 + bi.HessianCSC.size(0)] =
            (bi.NodeMap[0] + static_cast<double>(m4[i1])) - 1.0;
    }
    for (i1 = 0; i1 < last; i1++) {
        bi.HessianCSC[(i + i1) + bi.HessianCSC.size(0) * 2] = b_I[i1];
    }
    maxNodeDim = (args_nodeMap[static_cast<int>(poses.NumRowBlocks) - 1] +
                  args_nodeDims[static_cast<int>(poses.NumRowBlocks) - 1]) -
                 1.0;
    if (maxNodeDim < 1.0) {
        loop_ub = 0;
    } else {
        loop_ub = static_cast<int>(maxNodeDim);
    }
    gradient.set_size(loop_ub);
    for (i = 0; i < loop_ub; i++) {
        gradient[i] = bi.Gradient[i];
    }
    if (numEntries1 - 1.0 < 1.0) {
        loop_ub = 0;
    } else {
        loop_ub = static_cast<int>(numEntries1 - 1.0);
    }
    b_bi.set_size(loop_ub);
    c_bi.set_size(loop_ub);
    d_bi.set_size(loop_ub);
    for (i = 0; i < loop_ub; i++) {
        b_bi[i] = bi.HessianCSC[i];
        c_bi[i] = bi.HessianCSC[i + bi.HessianCSC.size(0)];
        d_bi[i] = bi.HessianCSC[i + bi.HessianCSC.size(0) * 2];
    }
    b_sparse(b_bi, c_bi, d_bi, hessian);
    return cost;
}

}  // namespace internal
}  // namespace algs
}  // namespace nav
}  // namespace coder
}  // namespace SlamGraph2D

///
/// File trailer for PoseGraphHelpers.cpp
///
/// [EOF]
///
