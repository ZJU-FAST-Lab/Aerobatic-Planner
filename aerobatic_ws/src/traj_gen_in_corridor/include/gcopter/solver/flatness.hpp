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

            tilt_den_sqr = 2.0 * (1.0 + z2);
            tilt_den = sqrt(tilt_den_sqr);
            tilt0 = 0.5 * tilt_den;
            tilt1 = -z1 / tilt_den;
            tilt2 = z0 / tilt_den;
            // if (z2 + 1.0 < 1.0e-5)
            // {
            //     // return false;
            //     std::cout << "1";
            //     dtilt_den = 0.0;
            // }
            // else
            // {
            // std::cout << "0";
            dtilt_den = dz2 / tilt_den;
            // }
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

            // dir_xb = dir - dir.dot(zb) * zb;
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
            // if (dir_xb_sqr_norm < 1.0e-5)
            // {
            //     std::cout << "2";
            //     std::cout << "dir = " << dir.transpose() << std::endl;
            //     std::cout << "dir_xb = " << dir_xb0 << ", " << dir_xb1 << ", " << dir_xb2 << ", " << std::endl;;
            //     return false;
            //     ddir_xb_norm = 0.0;
            // }
            // else
            // {
            // std::cout << "0";
            ddir_xb_norm = ddir_xb_sqr_norm / (2.0 * dir_xb_norm);
            // }

            // sign = xb.cross(dir_xb).dot(zb) > 0 ? 1.0 : -1.0;
            xb_cross_dir_xb0 = x1 * dir_xb2 - x2 * dir_xb1;
            xb_cross_dir_xb1 = x2 * dir_xb0 - x0 * dir_xb2;
            xb_cross_dir_xb2 = x0 * dir_xb1 - x1 * dir_xb0;

            if ((z0 * xb_cross_dir_xb0 + z1 * xb_cross_dir_xb1 + z2 * xb_cross_dir_xb2) > 0)
                sign = 1.0;
            else
                sign = -1.0;
            // xvdot = xb.dot(dir_xb.normalized());
            xb_dot_dir_xb = (x0 * dir_xb0 + x1 * dir_xb1 + x2 * dir_xb2) / dir_xb_norm;
            // yaw = acos(xvdot) * sign;
            psi = acos(xb_dot_dir_xb) * sign;

            dxb_dot_dir_xb = (dir_xb0 * dx0 + x0 * ddir_xb0 + dir_xb1 * dx1 + x1 * ddir_xb1 + dir_xb2 * dx2 + x2 * ddir_xb2 - xb_dot_dir_xb * ddir_xb_norm) / dir_xb_norm;
            // yaw = acos(xvdot) * sign;
            // if (1.0 - xb_dot_dir_xb < 1.0e-5 || xb_dot_dir_xb + 1.0 < 1.0e-5)
            // {
            //     std::cout << "3";
            //     std::cout << "xb = " << x0 << "," << x1 << "," << x2 << " | " << dir.transpose() <<  std::endl;
            //     dpsi = -sign * (dxb_dot_dir_xb / sqrt(1.0 - xb_dot_dir_xb * xb_dot_dir_xb));
            //     std::cout << "dpsi = " << dpsi << std::endl;
            //     return false;
            //     dpsi = 0.0;
            // }
            // else
            // {
            // std::cout << "0";
            dpsi = -sign * (dxb_dot_dir_xb / sqrt(1.0 - xb_dot_dir_xb * xb_dot_dir_xb));
            // }

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

            // std::cout << " | ";
            return true;
        }

        inline void backward(const Eigen::Vector3d &pos_grad,
                             const Eigen::Vector3d &vel_grad,
                             const double &thr_grad,
                             const Eigen::Vector4d &quat_grad,
                             const Eigen::Vector3d &omg_grad,
                             Eigen::Vector3d &pos_total_grad,
                             Eigen::Vector3d &vel_total_grad,
                             Eigen::Vector3d &acc_total_grad,
                             Eigen::Vector3d &jer_total_grad,
                             double &psi_total_grad,
                             double &dpsi_total_grad)
        {
            tilt0b = s_half_psi * (quat_grad(3)) + c_half_psi * (quat_grad(0));
            head3b = tilt0 * (quat_grad(3)) + tilt2 * (quat_grad(1)) - tilt1 * (quat_grad(2));
            tilt2b = c_half_psi * (quat_grad(2)) + s_half_psi * (quat_grad(1));
            head0b = tilt2 * (quat_grad(2)) + tilt1 * (quat_grad(1)) + tilt0 * (quat_grad(0));
            tilt1b = c_half_psi * (quat_grad(1)) - s_half_psi * (quat_grad(2));
            tilt_den_sqr = tilt_den * tilt_den;
            tilt_denb = (z1 * tilt1b - z0 * tilt2b) / tilt_den_sqr + 0.5 * tilt0b;
            omg_termb = -((z0 * c_psi + z1 * s_psi) * (omg_grad(1))) -
                        (z0 * s_psi - z1 * c_psi) * (omg_grad(0));
            tempb = omg_grad(2) / omg_den;
            dpsi_total_grad = omg_grad(2);
            z1b = dz0 * tempb;
            dz0b = z1 * tempb + c_psi * (omg_grad(1)) + s_psi * (omg_grad(0));
            z0b = -(dz1 * tempb);
            dz1b = s_psi * (omg_grad(1)) - z0 * tempb - c_psi * (omg_grad(0));
            omg_denb = -((z1 * dz0 - z0 * dz1) * tempb / omg_den) -
                       dz2 * omg_termb / (omg_den * omg_den);
            tempb = -(omg_term * (omg_grad(1)));
            c_psib = dz0 * (omg_grad(1)) + z0 * tempb;
            s_psib = dz1 * (omg_grad(1)) + z1 * tempb;
            z0b += c_psi * tempb;
            z1b += s_psi * tempb;
            tempb = -(omg_term * (omg_grad(0)));
            s_psib += dz0 * (omg_grad(0)) + z0 * tempb;
            c_psib += -dz1 * (omg_grad(0)) - z1 * tempb;
            z0b += s_psi * tempb + tilt2b / tilt_den + f_term0 * (thr_grad);
            z1b += -c_psi * tempb - tilt1b / tilt_den + f_term1 * (thr_grad);
            dz2b = omg_termb / omg_den;
            z2b = omg_denb + tilt_denb / tilt_den + f_term2 * (thr_grad);
            psi_total_grad = c_psi * s_psib + 0.5 * c_half_psi * head3b -
                             s_psi * c_psib - 0.5 * s_half_psi * head0b;
            f_term0b = z0 * (thr_grad);
            f_term1b = z1 * (thr_grad);
            f_term2b = z2 * (thr_grad);
            ng02b = dz_term0 * dz2b + dz_term2 * dz0b;
            dz_term0b = ng02 * dz2b + ng01 * dz1b + ng00 * dz0b;
            ng12b = dz_term1 * dz2b + dz_term2 * dz1b;
            dz_term1b = ng12 * dz2b + ng11 * dz1b + ng01 * dz0b;
            ng22b = dz_term2 * dz2b;
            dz_term2b = ng22 * dz2b + ng12 * dz1b + ng02 * dz0b;
            ng01b = dz_term0 * dz1b + dz_term1 * dz0b;
            ng11b = dz_term1 * dz1b;
            ng00b = dz_term0 * dz0b;
            jer_total_grad(2) = dz_term2b;
            dw2b = dh_over_m * dz_term2b;
            jer_total_grad(1) = dz_term1b;
            dw1b = dh_over_m * dz_term1b;
            jer_total_grad(0) = dz_term0b;
            dw0b = dh_over_m * dz_term0b;
            tempb = cp * (v2 * dw2b + v1 * dw1b + v0 * dw0b) / cp_term;
            acc_total_grad(2) = mass * f_term2b + w_term * dw2b + v2 * tempb;
            acc_total_grad(1) = mass * f_term1b + w_term * dw1b + v1 * tempb;
            acc_total_grad(0) = mass * f_term0b + w_term * dw0b + v0 * tempb;
            vel_total_grad(2) = dw_term * dw2b + a2 * tempb;
            vel_total_grad(1) = dw_term * dw1b + a1 * tempb;
            vel_total_grad(0) = dw_term * dw0b + a0 * tempb;
            cp_termb = -(v_dot_a * tempb / cp_term);
            tempb = ng22b / ng_den;
            zu_sqr0b = tempb;
            zu_sqr1b = tempb;
            ng_denb = -((zu_sqr0 + zu_sqr1) * tempb / ng_den);
            zu12b = -(ng12b / ng_den);
            tempb = ng11b / ng_den;
            ng_denb += zu12 * ng12b / (ng_den * ng_den) -
                       (zu_sqr0 + zu_sqr2) * tempb / ng_den;
            zu_sqr0b += tempb;
            zu_sqr2b = tempb;
            zu02b = -(ng02b / ng_den);
            zu01b = -(ng01b / ng_den);
            tempb = ng00b / ng_den;
            ng_denb += zu02 * ng02b / (ng_den * ng_den) +
                       zu01 * ng01b / (ng_den * ng_den) -
                       (zu_sqr1 + zu_sqr2) * tempb / ng_den;
            zu_normb = zu_sqr_norm * ng_denb -
                       (zu2 * z2b + zu1 * z1b + zu0 * z0b) / zu_sqr_norm;
            zu_sqr_normb = zu_norm * ng_denb + zu_normb / (2.0 * zu_norm);
            tempb += zu_sqr_normb;
            zu_sqr1b += tempb;
            zu_sqr2b += tempb;
            zu2b = z2b / zu_norm + zu0 * zu02b + zu1 * zu12b + 2 * zu2 * zu_sqr2b;
            w2b = dv * f_term2b + dh_over_m * zu2b;
            zu1b = z1b / zu_norm + zu2 * zu12b + zu0 * zu01b + 2 * zu1 * zu_sqr1b;
            w1b = dv * f_term1b + dh_over_m * zu1b;
            zu_sqr0b += zu_sqr_normb;
            zu0b = z0b / zu_norm + zu2 * zu02b + zu1 * zu01b + 2 * zu0 * zu_sqr0b;
            w0b = dv * f_term0b + dh_over_m * zu0b;
            w_termb = a2 * dw2b + a1 * dw1b + a0 * dw0b +
                      v2 * w2b + v1 * w1b + v0 * w0b;
            acc_total_grad(2) += zu2b;
            acc_total_grad(1) += zu1b;
            acc_total_grad(0) += zu0b;
            cp_termb += cp * w_termb;
            v_sqr_normb = cp_termb / (2.0 * cp_term);
            vel_total_grad(2) += w_term * w2b + 2 * v2 * v_sqr_normb + vel_grad(2);
            vel_total_grad(1) += w_term * w1b + 2 * v1 * v_sqr_normb + vel_grad(1);
            vel_total_grad(0) += w_term * w0b + 2 * v0 * v_sqr_normb + vel_grad(0);
            pos_total_grad(2) = pos_grad(2);
            pos_total_grad(1) = pos_grad(1);
            pos_total_grad(0) = pos_grad(0);
        }

        inline void backward(const Eigen::Vector3d &dir,
                             const Eigen::Vector3d &ddir,
                             const Eigen::Vector3d &pos_grad,
                             const Eigen::Vector3d &vel_grad,
                             const double &thr_grad,
                             const Eigen::Vector4d &quat_grad,
                             const Eigen::Vector3d &omg_grad,
                             Eigen::Vector3d &pos_total_grad,
                             Eigen::Vector3d &vel_total_grad,
                             Eigen::Vector3d &acc_total_grad,
                             Eigen::Vector3d &jer_total_grad,
                             Eigen::Vector3d &dir_total_grad,
                             Eigen::Vector3d &ddir_total_grad)
        {
            s_half_psib = tilt0 * quat_grad(3) - tilt1 * quat_grad(2) + tilt2 * quat_grad(1);
            c_half_psib = tilt2 * quat_grad(2) + tilt1 * quat_grad(1) + tilt0 * quat_grad(0);
            omg_termb = (z1 * c_psi - z0 * s_psi) * omg_grad(0) - (z0 * c_psi + z1 * s_psi) * omg_grad(1);
            omg_temp0b = -(omg_term * omg_grad(0));
            omg_temp1b = -(omg_term * omg_grad(1));
            omg_temp2b = omg_grad(2) / omg_den;
            omg_denb = (z0 * dz1 - z1 * dz0) * omg_temp2b / omg_den - dz2 * omg_termb / (omg_den * omg_den);
            z0b = s_psi * omg_temp0b + c_psi * omg_temp1b - dz1 * omg_temp2b;
            z1b = s_psi * omg_temp1b - c_psi * omg_temp0b + dz0 * omg_temp2b;
            z2b = omg_denb;
            dz0b = s_psi * omg_grad(0) + c_psi * omg_grad(1) + z1 * omg_temp2b;
            dz1b = s_psi * omg_grad(1) - c_psi * omg_grad(0) - z0 * omg_temp2b;
            dz2b = omg_termb / omg_den;
            c_psib = dz0 * omg_grad(1) - dz1 * omg_grad(0) + z0 * omg_temp1b - z1 * omg_temp0b;
            s_psib = dz1 * omg_grad(1) + dz0 * omg_grad(0) + z1 * omg_temp1b + z0 * omg_temp0b;
            psib = c_psi * s_psib + 0.5 * c_half_psi * s_half_psib - s_psi * c_psib - 0.5 * s_half_psi * c_half_psib;
            temp = sqrt(1.0 - xb_dot_dir_xb * xb_dot_dir_xb);
            dxb_dot_dir_xbb = -(sign * omg_grad(2) / temp);
            dxb_dot_dir_xb_term = dxb_dot_dir_xbb / dir_xb_norm;
            xb_dot_dir_xbb = xb_dot_dir_xb * dxb_dot_dir_xb * dxb_dot_dir_xbb / (temp * temp) - ddir_xb_norm * dxb_dot_dir_xb_term - sign * psib / temp;
            xb_dot_dir_xb_term = xb_dot_dir_xbb / dir_xb_norm;
            ddir_xb_normb = -xb_dot_dir_xb * dxb_dot_dir_xb_term;
            ddir_xb_sqr_normb = ddir_xb_normb / (2.0 * dir_xb_norm);
            x0b = ddir_xb0 * dxb_dot_dir_xb_term + dir_xb0 * xb_dot_dir_xb_term;
            x1b = ddir_xb1 * dxb_dot_dir_xb_term + dir_xb1 * xb_dot_dir_xb_term;
            x2b = ddir_xb2 * dxb_dot_dir_xb_term + dir_xb2 * xb_dot_dir_xb_term;
            dx0b = dir_xb0 * dxb_dot_dir_xb_term;
            dx1b = dir_xb1 * dxb_dot_dir_xb_term;
            dx2b = dir_xb2 * dxb_dot_dir_xb_term;
            temp = dir_xb0 * dx0 + x0 * ddir_xb0 + dir_xb1 * dx1 + x1 * ddir_xb1 + dir_xb2 * dx2 + x2 * ddir_xb2;
            dir_xb_normb = -(temp - xb_dot_dir_xb * ddir_xb_norm) * dxb_dot_dir_xb_term / dir_xb_norm - (x0 * dir_xb0 + x1 * dir_xb1 + x2 * dir_xb2) * xb_dot_dir_xb_term / dir_xb_norm - ddir_xb_sqr_norm * ddir_xb_sqr_normb / dir_xb_norm;
            dir_xb_sqr_normb = dir_xb_normb / (2.0 * sqrt(dir_xb_sqr_norm));
            qwb = -dqy * 2 * dx2b - qy * 2 * x2b;
            qxb = dqy * 2 * dx1b + qy * 2 * x1b;
            qyb = -dqw * 2 * dx2b + dqx * 2 * dx1b + qx * 2 * x1b - dqy * 4 * dx0b - qw * 2 * x2b - 4 * qy * x0b;
            dqwb = -qy * 2 * dx2b;
            dqxb = qy * 2 * dx1b;
            dqyb = -qw * 2 * dx2b + qx * 2 * dx1b - qy * 4 * dx0b;
            dir_xb0b = dx0 * dxb_dot_dir_xb_term + x0 * xb_dot_dir_xb_term + 2 * ddir_xb0 * ddir_xb_sqr_normb + 2 * dir_xb0 * dir_xb_sqr_normb;
            dir_xb1b = dx1 * dxb_dot_dir_xb_term + x1 * xb_dot_dir_xb_term + 2 * ddir_xb1 * ddir_xb_sqr_normb + 2 * dir_xb1 * dir_xb_sqr_normb;
            dir_xb2b = dx2 * dxb_dot_dir_xb_term + x2 * xb_dot_dir_xb_term + 2 * ddir_xb2 * ddir_xb_sqr_normb + 2 * dir_xb2 * dir_xb_sqr_normb;
            ddir_xb0b = x0 * dxb_dot_dir_xb_term + 2 * dir_xb0 * ddir_xb_sqr_normb;
            ddir_xb1b = x1 * dxb_dot_dir_xb_term + 2 * dir_xb1 * ddir_xb_sqr_normb;
            ddir_xb2b = x2 * dxb_dot_dir_xb_term + 2 * dir_xb2 * ddir_xb_sqr_normb;
            ddir_dot_zbb = -(z2 * ddir_xb2b) - z1 * ddir_xb1b - z0 * ddir_xb0b;
            dir_dot_zbb = -(dz2 * ddir_xb2b) - dz1 * ddir_xb1b - dz0 * ddir_xb0b - z2 * dir_xb2b - z1 * dir_xb1b - z0 * dir_xb0b;
            ddir_total_grad(0) = ddir_xb0b + z0 * ddir_dot_zbb;
            ddir_total_grad(1) = ddir_xb1b + z1 * ddir_dot_zbb;
            ddir_total_grad(2) = ddir_xb2b + z2 * ddir_dot_zbb;
            dir_total_grad(0) = dz0 * ddir_dot_zbb + dir_xb0b + z0 * dir_dot_zbb;
            dir_total_grad(1) = dz1 * ddir_dot_zbb + dir_xb1b + z1 * dir_dot_zbb;
            dir_total_grad(2) = dz2 * ddir_dot_zbb + dir_xb2b + z2 * dir_dot_zbb;
            tilt0b = s_half_psi * quat_grad(3) + c_half_psi * quat_grad(0) + qwb;
            tilt1b = c_half_psi * quat_grad(1) - s_half_psi * quat_grad(2) + qxb;
            tilt2b = c_half_psi * quat_grad(2) + s_half_psi * quat_grad(1) + qyb;
            dtilt_denb = (z1 * dqxb - z0 * dqyb) / tilt_den_sqr + 0.5 * dqwb;
            tilt_denb = -(dz0 - z0 * dtilt_den / tilt_den) * dqyb / tilt_den_sqr + z0 * dtilt_den / tilt_den * dqyb / tilt_den_sqr + (dz1 - z1 * dtilt_den / tilt_den) * dqxb / tilt_den_sqr - z1 * dtilt_den / tilt_den * dqxb / tilt_den_sqr - (dz2 * dtilt_denb - z1 * tilt1b + z0 * tilt2b) / tilt_den_sqr + 0.5 * tilt0b;
            dz0b = dz0b + dir(0) * ddir_dot_zbb - dir_dot_zb * ddir_xb0b + dqyb / tilt_den;
            dz1b = dz1b + dir(1) * ddir_dot_zbb - dir_dot_zb * ddir_xb1b - dqxb / tilt_den;
            dz2b = dz2b + dir(2) * ddir_dot_zbb - dir_dot_zb * ddir_xb2b + dtilt_denb / tilt_den;
            z0b = z0b + ddir(0) * ddir_dot_zbb - ddir_dot_zb * ddir_xb0b + dir(0) * dir_dot_zbb - dir_dot_zb * dir_xb0b - dtilt_den * dqyb / tilt_den_sqr + tilt2b / tilt_den + f_term0 * thr_grad;
            z1b = z1b + ddir(1) * ddir_dot_zbb - ddir_dot_zb * ddir_xb1b + dir(1) * dir_dot_zbb - dir_dot_zb * dir_xb1b + dtilt_den * dqxb / tilt_den_sqr - tilt1b / tilt_den + f_term1 * thr_grad;
            z2b = z2b + ddir(2) * ddir_dot_zbb - ddir_dot_zb * ddir_xb2b + dir(2) * dir_dot_zbb - dir_dot_zb * dir_xb2b + tilt_denb / tilt_den + f_term2 * thr_grad;
            f_term0b = z0 * thr_grad;
            f_term1b = z1 * thr_grad;
            f_term2b = z2 * thr_grad;
            dz_term0b = ng02 * dz2b + ng01 * dz1b + ng00 * dz0b;
            dz_term1b = ng12 * dz2b + ng11 * dz1b + ng01 * dz0b;
            dz_term2b = ng22 * dz2b + ng12 * dz1b + ng02 * dz0b;
            ng00b = dz_term0 * dz0b;
            ng11b = dz_term1 * dz1b;
            ng22b = dz_term2 * dz2b;
            ng01b = dz_term0 * dz1b + dz_term1 * dz0b;
            ng02b = dz_term0 * dz2b + dz_term2 * dz0b;
            ng12b = dz_term1 * dz2b + dz_term2 * dz1b;
            dw0b = dh_over_m * dz_term0b;
            dw1b = dh_over_m * dz_term1b;
            dw2b = dh_over_m * dz_term2b;
            jer_total_grad(0) = dz_term0b;
            jer_total_grad(1) = dz_term1b;
            jer_total_grad(2) = dz_term2b;
            dw_termb = v2 * dw2b + v1 * dw1b + v0 * dw0b;
            v_dot_ab = cp * dw_termb / cp_term;
            zu00b = ng00b / ng_den;
            zu11b = ng11b / ng_den;
            zu22b = ng22b / ng_den;
            zu01b = -(ng01b / ng_den);
            zu02b = -(ng02b / ng_den);
            zu12b = -(ng12b / ng_den);
            ng_denb = zu12 * ng12b / (ng_den * ng_den) - (zu_sqr0 + zu_sqr1) * zu22b / ng_den
                    + zu02 * ng02b / (ng_den * ng_den) - (zu_sqr0 + zu_sqr2) * zu11b / ng_den
                    + zu01 * ng01b / (ng_den * ng_den) - (zu_sqr1 + zu_sqr2) * zu00b / ng_den;
            zu_normb = zu_sqr_norm * ng_denb - (zu2 * z2b + zu1 * z1b + zu0 * z0b) / zu_sqr_norm;
            zu_sqr_normb = zu_norm * ng_denb + zu_normb / (2.0 * zu_norm);
            zu_sqr0b = zu11b + zu22b + zu_sqr_normb;
            zu_sqr1b = zu22b + zu00b + zu_sqr_normb;
            zu_sqr2b = zu00b + zu11b + zu_sqr_normb;
            zu0b = z0b / zu_norm + zu2 * zu02b + zu1 * zu01b + 2 * zu0 * zu_sqr0b;
            zu1b = z1b / zu_norm + zu2 * zu12b + zu0 * zu01b + 2 * zu1 * zu_sqr1b;
            zu2b = z2b / zu_norm + zu0 * zu02b + zu1 * zu12b + 2 * zu2 * zu_sqr2b;
            w0b = dv * f_term0b + dh_over_m * zu0b;
            w1b = dv * f_term1b + dh_over_m * zu1b;
            w2b = dv * f_term2b + dh_over_m * zu2b;
            w_termb = a2 * dw2b + a1 * dw1b + a0 * dw0b + v2 * w2b + v1 * w1b + v0 * w0b;
            cp_termb = cp * w_termb - v_dot_a * v_dot_ab / cp_term;
            temp = cp_termb / (2.0 * sqrt(veps + v0 * v0 + v1 * v1 + v2 * v2));
            acc_total_grad(0) = mass * f_term0b + w_term * dw0b + v0 * v_dot_ab + zu0b;
            acc_total_grad(1) = mass * f_term1b + w_term * dw1b + v1 * v_dot_ab + zu1b;
            acc_total_grad(2) = mass * f_term2b + w_term * dw2b + v2 * v_dot_ab + zu2b;
            vel_total_grad(0) = vel_grad(0) + dw_term * dw0b + a0 * v_dot_ab + w_term * w0b + 2 * v0 * temp;
            vel_total_grad(1) = vel_grad(1) + dw_term * dw1b + a1 * v_dot_ab + w_term * w1b + 2 * v1 * temp;
            vel_total_grad(2) = vel_grad(2) + dw_term * dw2b + a2 * v_dot_ab + w_term * w2b + 2 * v2 * temp;
            pos_total_grad(0) = pos_grad(0);
            pos_total_grad(1) = pos_grad(1);
            pos_total_grad(2) = pos_grad(2);
            return;
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
        double tilt_den_sqr, tilt_den, tilt0, tilt1, tilt2, c_half_psi, s_half_psi;
        double c_psi, s_psi, omg_den, omg_term;
        double qw, qx, qy, qz, dqw, dqx, dqy;
        double dtilt0, dtilt1, dtilt2, dtilt_den;
        double x0, x1, x2, dx0, dx1, dx2;
        double ddir_dot_zb, ddir_xb0, ddir_xb1, ddir_xb2, ddir_xb_norm;
        double dir_dot_zb, dir_xb0, dir_xb1, dir_xb2, dir_xb_norm;
        double xb_dot_dir_xb, dxb_dot_dir_xb;
        double dir_xb_sqr_norm, ddir_xb_sqr_norm, sign;
        double xb_cross_dir_xb0, xb_cross_dir_xb1, xb_cross_dir_xb2;

        double tilt0b, tilt1b, tilt2b, head0b, head3b;
        double s_half_psib, c_half_psib, s_psib, c_psib, psib;
        double omg_denb, omg_termb, omg_temp0b, omg_temp1b, omg_temp2b, temp, tempb;
        double z0b, z1b, z2b, dz0b, dz1b, dz2b;
        double x0b, x1b, x2b, dx0b, dx1b, dx2b;
        double dir_xb0b, dir_xb1b, dir_xb2b, ddir_xb0b, ddir_xb1b, ddir_xb2b;
        double xb_dot_dir_xbb, dxb_dot_dir_xbb;
        double xb_dot_dir_xb_term, dxb_dot_dir_xb_term;
        double dir_xb_normb, ddir_xb_normb;
        double dir_xb_sqr_normb, ddir_xb_sqr_normb;
        double dir_dot_zbb, ddir_dot_zbb;
        double qxb, dqxb, qyb, dqyb, qwb, dqwb;
        double tilt_denb, dtilt_denb;
        double f_term0b, f_term1b, f_term2b;
        double dz_term0b, dz_term1b, dz_term2b;
        double ng00b, ng01b, ng02b, ng11b, ng12b, ng22b, ng_denb;
        double w0b, w1b, w2b, dw0b, dw1b, dw2b, w_termb, dw_termb;
        double v_dot_ab, cp_termb, v_sqr_normb;
        double zu0b, zu1b, zu2b, zu01b, zu02b, zu12b, zu00b, zu11b, zu22b;
        double zu_sqr0b, zu_sqr1b, zu_sqr2b, zu_normb, zu_sqr_normb;
    };
}

#endif
