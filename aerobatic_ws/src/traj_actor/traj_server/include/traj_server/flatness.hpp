/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef FLATNESS_HPP
#define FLATNESS_HPP

#include <Eigen/Eigen>

#include <cmath>

namespace flatness
{
    class FlatnessMap // See https://github.com/ZJU-FAST-Lab/GCOPTER/blob/main/misc/flatness.pdf
    {
    public:
        FlatnessMap() : mass(0.98), grav(9.81), dh(0), dv(0), cp(0), veps(0.0001) {}

        inline void reset(const double &vehicle_mass,
                          const double &gravitational_acceleration,
                          const double &horitonral_drag_coeff,
                          const double &vertical_drag_coeff,
                          const double &parasitic_drag_coeff,
                          const double &speed_smooth_factor)
        {
            mass = vehicle_mass;
            grav = gravitational_acceleration;
            dh = horitonral_drag_coeff;
            dv = vertical_drag_coeff;
            cp = parasitic_drag_coeff;
            veps = speed_smooth_factor;

            return;
        }

        inline bool forward(const Eigen::Vector3d &vel,
                            const Eigen::Vector3d &acc,
                            const Eigen::Vector3d &jer,
                            const double &psi,
                            const double &dpsi,
                            double &thr,
                            Eigen::Vector4d &quat,
                            Eigen::Vector3d &omg)
        {
            v0 = vel(0);
            v1 = vel(1);
            v2 = vel(2);
            a0 = acc(0);
            a1 = acc(1);
            a2 = acc(2);
            cp_term = sqrt(v0 * v0 + v1 * v1 + v2 * v2 + veps);
            w_term = 1.0 + cp * cp_term;
            w0 = w_term * v0;
            w1 = w_term * v1;
            w2 = w_term * v2;
            dh_over_m = dh / mass;
            zu0 = a0 + dh_over_m * w0;
            zu1 = a1 + dh_over_m * w1;
            zu2 = a2 + dh_over_m * w2 + grav;
            zu_sqr0 = zu0 * zu0;
            zu_sqr1 = zu1 * zu1;
            zu_sqr2 = zu2 * zu2;
            zu01 = zu0 * zu1;
            zu12 = zu1 * zu2;
            zu02 = zu0 * zu2;
            zu_sqr_norm = zu_sqr0 + zu_sqr1 + zu_sqr2;
            zu_norm = sqrt(zu_sqr_norm);
            z0 = zu0 / zu_norm;
            z1 = zu1 / zu_norm;
            z2 = zu2 / zu_norm;
            ng_den = zu_sqr_norm * zu_norm;
            ng00 = (zu_sqr1 + zu_sqr2) / ng_den;
            ng01 = -zu01 / ng_den;
            ng02 = -zu02 / ng_den;
            ng11 = (zu_sqr0 + zu_sqr2) / ng_den;
            ng12 = -zu12 / ng_den;
            ng22 = (zu_sqr0 + zu_sqr1) / ng_den;
            v_dot_a = v0 * a0 + v1 * a1 + v2 * a2;
            dw_term = cp * v_dot_a / cp_term;
            dw0 = w_term * a0 + dw_term * v0;
            dw1 = w_term * a1 + dw_term * v1;
            dw2 = w_term * a2 + dw_term * v2;
            dz_term0 = jer(0) + dh_over_m * dw0;
            dz_term1 = jer(1) + dh_over_m * dw1;
            dz_term2 = jer(2) + dh_over_m * dw2;
            dz0 = ng00 * dz_term0 + ng01 * dz_term1 + ng02 * dz_term2;
            dz1 = ng01 * dz_term0 + ng11 * dz_term1 + ng12 * dz_term2;
            dz2 = ng02 * dz_term0 + ng12 * dz_term1 + ng22 * dz_term2;
            f_term0 = mass * a0 + dv * w0;
            f_term1 = mass * a1 + dv * w1;
            f_term2 = mass * (a2 + grav) + dv * w2;
            thr = z0 * f_term0 + z1 * f_term1 + z2 * f_term2;
            tilt_den = sqrt(2.0 * (1.0 + z2));
            tilt0 = 0.5 * tilt_den;
            tilt1 = -z1 / tilt_den;
            tilt2 = z0 / tilt_den;
            c_half_psi = cos(0.5 * psi);
            s_half_psi = sin(0.5 * psi);
            quat(0) = tilt0 * c_half_psi;
            quat(1) = tilt1 * c_half_psi + tilt2 * s_half_psi;
            quat(2) = tilt2 * c_half_psi - tilt1 * s_half_psi;
            quat(3) = tilt0 * s_half_psi;
            c_psi = cos(psi);
            s_psi = sin(psi);
            omg_den = z2 + 1.0;
            omg_term = dz2 / omg_den;
            omg(0) = dz0 * s_psi - dz1 * c_psi -
                     (z0 * s_psi - z1 * c_psi) * omg_term;
            omg(1) = dz0 * c_psi + dz1 * s_psi -
                     (z0 * c_psi + z1 * s_psi) * omg_term;
            omg(2) = (z1 * dz0 - z0 * dz1) / omg_den + dpsi;
            return true;
        }

        inline bool forward(const Eigen::Vector3d &vel,
                            const Eigen::Vector3d &acc,
                            const Eigen::Vector3d &jer,
                            const Eigen::Vector3d &dir,
                            const Eigen::Vector3d &ddir,
                            double &psi,
                            double &dpsi,
                            double &thr,
                            Eigen::Vector4d &quat,
                            Eigen::Vector3d &omg)
        {
            v0 = vel(0);
            v1 = vel(1);
            v2 = vel(2);
            a0 = acc(0);
            a1 = acc(1);
            a2 = acc(2);
            cp_term = sqrt(v0 * v0 + v1 * v1 + v2 * v2 + veps);
            w_term = 1.0 + cp * cp_term;
            w0 = w_term * v0;
            w1 = w_term * v1;
            w2 = w_term * v2;
            dh_over_m = dh / mass;
            zu0 = a0 + dh_over_m * w0;
            zu1 = a1 + dh_over_m * w1;
            zu2 = a2 + dh_over_m * w2 + grav;
            zu_sqr0 = zu0 * zu0;
            zu_sqr1 = zu1 * zu1;
            zu_sqr2 = zu2 * zu2;
            zu01 = zu0 * zu1;
            zu12 = zu1 * zu2;
            zu02 = zu0 * zu2;
            zu_sqr_norm = zu_sqr0 + zu_sqr1 + zu_sqr2;
            zu_norm = sqrt(zu_sqr_norm);
            z0 = zu0 / zu_norm;
            z1 = zu1 / zu_norm;
            z2 = zu2 / zu_norm;
            ng_den = zu_sqr_norm * zu_norm;
            ng00 = (zu_sqr1 + zu_sqr2) / ng_den;
            ng01 = -zu01 / ng_den;
            ng02 = -zu02 / ng_den;
            ng11 = (zu_sqr0 + zu_sqr2) / ng_den;
            ng12 = -zu12 / ng_den;
            ng22 = (zu_sqr0 + zu_sqr1) / ng_den;
            v_dot_a = v0 * a0 + v1 * a1 + v2 * a2;
            dw_term = cp * v_dot_a / cp_term;
            dw0 = w_term * a0 + dw_term * v0;
            dw1 = w_term * a1 + dw_term * v1;
            dw2 = w_term * a2 + dw_term * v2;
            dz_term0 = jer(0) + dh_over_m * dw0;
            dz_term1 = jer(1) + dh_over_m * dw1;
            dz_term2 = jer(2) + dh_over_m * dw2;
            dz0 = ng00 * dz_term0 + ng01 * dz_term1 + ng02 * dz_term2;
            dz1 = ng01 * dz_term0 + ng11 * dz_term1 + ng12 * dz_term2;
            dz2 = ng02 * dz_term0 + ng12 * dz_term1 + ng22 * dz_term2;

            f_term0 = mass * a0 + dv * w0;
            f_term1 = mass * a1 + dv * w1;
            f_term2 = mass * (a2 + grav) + dv * w2;
            thr = z0 * f_term0 + z1 * f_term1 + z2 * f_term2;

            tilt_den = sqrt(2.0 * (1.0 + z2));
            tilt0 = 0.5 * tilt_den;
            tilt1 = -z1 / tilt_den;
            tilt2 = z0 / tilt_den;
            dtilt_den = dz2 / tilt_den;
            dtilt0 = 0.5 * dtilt_den;
            dtilt1 = -((dz1 - z1 * dtilt_den / tilt_den) / tilt_den);
            dtilt2 = (dz0 - z0 * dtilt_den / tilt_den) / tilt_den;

            qw = tilt0;
            qx = tilt1;
            qy = tilt2;
            dqw = dtilt0;
            dqx = dtilt1;
            dqy = dtilt2;

            x0 = 1 - 2 * qy * qy;
            x1 = 2 * qx * qy;
            x2 = -2 * qw * qy;
            dx0 = -4 * dqy * qy;
            dx1 = 2 * (dqx * qy + dqy * qx);
            dx2 = -2 * (dqw * qy + dqy * qw);

            dir_dot_zb = dir(0) * z0 + dir(1) * z1 + dir(2) * z2;
            dir_xb0 = dir(0) - dir_dot_zb * z0;
            dir_xb1 = dir(1) - dir_dot_zb * z1;
            dir_xb2 = dir(2) - dir_dot_zb * z2;
            dir_xb_sqr_norm = dir_xb0 * dir_xb0 + dir_xb1 * dir_xb1 + dir_xb2 * dir_xb2;
            dir_xb_norm = sqrt(dir_xb_sqr_norm);
            ddir_dot_zb = z0 * ddir(0) + dir(0) * dz0 + z1 * ddir(1) + dir(1) * dz1 + z2 * ddir(2) + dir(2) * dz2;
            ddir_xb0 = ddir(0) - z0 * ddir_dot_zb - dir_dot_zb * dz0;
            ddir_xb1 = ddir(1) - z1 * ddir_dot_zb - dir_dot_zb * dz1;
            ddir_xb2 = ddir(2) - z2 * ddir_dot_zb - dir_dot_zb * dz2;
            ddir_xb_sqr_norm = 2 * dir_xb0 * ddir_xb0 + 2 * dir_xb1 * ddir_xb1 + 2 * dir_xb2 * ddir_xb2;
            ddir_xb_norm = ddir_xb_sqr_norm / (2.0 * dir_xb_norm);
            xb_cross_dir_xb0 = x1 * dir_xb2 - x2 * dir_xb1;
            xb_cross_dir_xb1 = x2 * dir_xb0 - x0 * dir_xb2;
            xb_cross_dir_xb2 = x0 * dir_xb1 - x1 * dir_xb0;

            if ((z0 * xb_cross_dir_xb0 + z1 * xb_cross_dir_xb1 + z2 * xb_cross_dir_xb2) > 0)
                sign = 1.0;
            else
                sign = -1.0;
            xb_dot_dir_xb = (x0 * dir_xb0 + x1 * dir_xb1 + x2 * dir_xb2) / dir_xb_norm;
            psi = acos(xb_dot_dir_xb) * sign;
            dxb_dot_dir_xb = (dir_xb0 * dx0 + x0 * ddir_xb0 + dir_xb1 * dx1 + x1 * ddir_xb1 + dir_xb2 * dx2 + x2 * ddir_xb2 - xb_dot_dir_xb * ddir_xb_norm) / dir_xb_norm;
            dpsi = -sign * (dxb_dot_dir_xb / sqrt(1.0 - xb_dot_dir_xb * xb_dot_dir_xb));

            c_half_psi = cos(0.5 * psi);
            s_half_psi = sin(0.5 * psi);
            quat(0) = tilt0 * c_half_psi;
            quat(1) = tilt1 * c_half_psi + tilt2 * s_half_psi;
            quat(2) = tilt2 * c_half_psi - tilt1 * s_half_psi;
            quat(3) = tilt0 * s_half_psi;

            c_psi = cos(psi);
            s_psi = sin(psi);
            omg_den = z2 + 1.0;
            omg_term = dz2 / omg_den;
            omg(0) = dz0 * s_psi - dz1 * c_psi -
                     (z0 * s_psi - z1 * c_psi) * omg_term;
            omg(1) = dz0 * c_psi + dz1 * s_psi -
                     (z0 * c_psi + z1 * s_psi) * omg_term;
            omg(2) = (z1 * dz0 - z0 * dz1) / omg_den + dpsi;
            return true;
        }

    private:
        double mass, grav, dh, dv, cp, veps;

        double w0, w1, w2, dw0, dw1, dw2;
        double v0, v1, v2, a0, a1, a2, v_dot_a;
        double z0, z1, z2, dz0, dz1, dz2;
        double cp_term, w_term, dh_over_m;
        double zu_sqr_norm, zu_norm, zu0, zu1, zu2;
        double zu_sqr0, zu_sqr1, zu_sqr2, zu01, zu12, zu02;
        double ng00, ng01, ng02, ng11, ng12, ng22, ng_den;
        double dw_term, dz_term0, dz_term1, dz_term2, f_term0, f_term1, f_term2;
        double tilt_den, tilt0, tilt1, tilt2, c_half_psi, s_half_psi;
        double c_psi, s_psi, omg_den, omg_term;
        double qw, qx, qy, qz, dqw, dqx, dqy;
        double dtilt0, dtilt1, dtilt2, dtilt_den;
        double x0, x1, x2, dx0, dx1, dx2;
        double ddir_dot_zb, ddir_xb0, ddir_xb1, ddir_xb2, ddir_xb_norm;
        double dir_dot_zb, dir_xb0, dir_xb1, dir_xb2, dir_xb_norm;
        double xb_dot_dir_xb, dxb_dot_dir_xb;
        double dir_xb_sqr_norm, ddir_xb_sqr_norm;
        double sign;
        double xb_cross_dir_xb0, xb_cross_dir_xb1, xb_cross_dir_xb2;
    };
}

#endif
