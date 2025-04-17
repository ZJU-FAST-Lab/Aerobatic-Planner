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

#ifndef GCOPTER_HPP
#define GCOPTER_HPP

#include "gcopter/minco.hpp"
#include "gcopter/solver/lbfgs.hpp"
#include "traj_gen_in_corridor/visualizer.hpp"
#include "ros/ros.h"
#include "quadrotor_msgs/OptCostDebug.h"
#include "gcopter/trajectory.hpp"
#include "gcopter/solver/flatness.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>

namespace gcopter
{

    class GCOPTER
    {
    public:
        typedef Eigen::Matrix3Xd PolyhedronV;
        typedef Eigen::Matrix<double, 6, -1> PolyhedronH;
        typedef std::vector<PolyhedronV> PolyhedraV;
        typedef std::vector<PolyhedronH> PolyhedraH;

    private:
#if TRAJ_ORDER == 3
        minco::MINCO_S2NU minco;
#elif TRAJ_ORDER == 5
        minco::MINCO_S3NU minco;
#elif TRAJ_ORDER == 7
        minco::MINCO_S4NU minco;
#endif

        bool isdebug;
        bool isStartOptFlip;
        bool isStartOptOmgZ;
        bool isoptsettime;
        bool isStartOptAtt;
        flatness::FlatnessMap flatmap;
        Visualizer *visualTraj;
        ros::Publisher *optPub;
        Trajectory<TRAJ_ORDER> Traj;
        quadrotor_msgs::OptCostDebug optCost;

        double rho;                           
        Eigen::Matrix<double, 3, 4> headPVAJ; 
        Eigen::Matrix<double, 3, 4> tailPVAJ; 

        PolyhedraV vPolytopes;      
        PolyhedraH hPolytopes;      
        Eigen::Matrix3Xd shortPath; 

        Eigen::VectorXd usekeypos; 
        Eigen::Matrix3Xd set_att;  
        Eigen::Matrix3Xd set_pos;  
        Eigen::VectorXd set_time;  
        Eigen::Matrix3Xd flip_pos; 
        Eigen::VectorXd flip_err;  

        Eigen::VectorXi pieceIdx;
        Eigen::VectorXi vPolyIdx;
        Eigen::VectorXi hPolyIdx;

        int polyN;       
        int pieceN;      
        int optCount;    
        int spatialDim;  
        int temporalDim; 
        int flipDim;     

        double smoothEps;            
        int integralRes;             
        Eigen::VectorXi resVector;   
        int flipRes;                 
        Eigen::VectorXd magnitudeBd; 
        Eigen::VectorXd penaltyWt;   
        double allocSpeed;           
        double deltaT;               

        double total_cost, minctrl_cost, time_cost;
        double pos_cost, vel_cost, acc_cost, omg_cost;
        Eigen::VectorXd flip_cost, flip_pos_cost;

        lbfgs::lbfgs_parameter_t lbfgs_params; 

        Eigen::Vector3i totalOptCount;          
        Eigen::Matrix3Xd points;                
        Eigen::VectorXd times;                  
        Eigen::VectorXd flipTimes;              
        Eigen::Matrix3Xd gradByPoints;          
        Eigen::VectorXd gradByTimes;            
        Eigen::VectorXd gradByFliptimes;        
        Eigen::MatrixX3d partialGradByCoeffs;   
        Eigen::VectorXd partialGradByTimes;     
        Eigen::VectorXd partialGradByFliptimes; 

    private:
        static inline void forwardFlipTime(const Eigen::VectorXd &flipTime,
                                           const double &deltaT,
                                           Eigen::VectorXd &set_alpha)
        {
            Eigen::VectorXd norm_T = flipTime / (flipTime.sum() + 1);
            const int N = norm_T.size();
            set_alpha.resize(N);
            for (int i = 0; i < N; i++)
            {
                set_alpha(i) = (1 - 2 * deltaT) * norm_T.head(i + 1).sum() + deltaT;
            }
        }

        static inline void backwardGradFlipTime(const Eigen::VectorXd &set_alpha,
                                                const Eigen::VectorXd &flipTime,
                                                const Eigen::VectorXd &gradAlpha,
                                                const double &deltaT,
                                                Eigen::VectorXd &gradFlipTime)
        {
            const int N = set_alpha.size();
            double normsumpT;
            const double sumT = flipTime.sum() + 1, sumT2 = sumT * sumT;
            gradFlipTime.resize(N);
            gradFlipTime.setZero();

            for (int i = 0; i < N; i++)
            { // alpha
                for (int j = 0; j < N; j++)
                { // fliptime
                    for (int m = 0; m < N; m++)
                    {
                        normsumpT = 0;
                        if (m <= i)
                        {
                            normsumpT = -flipTime(m) / sumT2;
                            if (m == j)
                            {
                                normsumpT += 1 / sumT;
                            }
                            else
                            {
                            }
                            gradFlipTime(j) += gradAlpha(i) * (1 - 2 * deltaT) * normsumpT;
                        }
                    }
                }
            }
        }

        static inline void forwardT(const Eigen::VectorXd &tau,
                                    Eigen::VectorXd &T)
        {
            const int sizeTau = tau.size();
            T.resize(sizeTau);
            for (int i = 0; i < sizeTau; i++)
            {
                T(i) = tau(i) > 0.0
                           ? ((0.5 * tau(i) + 1.0) * tau(i) + 1.0)
                           : 1.0 / ((0.5 * tau(i) - 1.0) * tau(i) + 1.0);
            }
            return;
        }

        template <typename EIGENVEC>
        static inline void backwardT(const Eigen::VectorXd &T,
                                     EIGENVEC &tau)
        {
            const int sizeT = T.size();
            tau.resize(sizeT);
            for (int i = 0; i < sizeT; i++)
            {
                tau(i) = T(i) > 1.0
                             ? (sqrt(2.0 * T(i) - 1.0) - 1.0)
                             : (1.0 - sqrt(2.0 / T(i) - 1.0));
            }

            return;
        }

        template <typename EIGENVEC>
        static inline void backwardGradT(const Eigen::VectorXd &tau,
                                         const Eigen::VectorXd &gradT,
                                         EIGENVEC &gradTau)
        {
            const int sizeTau = tau.size();
            gradTau.resize(sizeTau);
            double denSqrt;
            for (int i = 0; i < sizeTau; i++)
            {
                if (tau(i) > 0)
                {
                    gradTau(i) = gradT(i) * (tau(i) + 1.0);
                }
                else
                {
                    denSqrt = (0.5 * tau(i) - 1.0) * tau(i) + 1.0;
                    gradTau(i) = gradT(i) * (1.0 - tau(i)) / (denSqrt * denSqrt);
                }
            }

            return;
        }

        static inline void forwardP(const Eigen::VectorXd &xi,
                                    const Eigen::VectorXi &vIdx,
                                    const PolyhedraV &vPolys,
                                    Eigen::Matrix3Xd &P)
        {
            const int sizeP = vIdx.size();
            P.resize(3, sizeP);
            Eigen::VectorXd q;
            for (int i = 0, j = 0, k, l; i < sizeP; i++, j += k)
            {
                l = vIdx(i);
                k = vPolys[l].cols();
                q = xi.segment(j, k).normalized().head(k - 1);
                P.col(i) = vPolys[l].rightCols(k - 1) * q.cwiseProduct(q) +
                           vPolys[l].col(0);
            }
            return;
        }

        static inline double costTinyNLS(void *ptr,
                                         const Eigen::VectorXd &xi,
                                         Eigen::VectorXd &gradXi)
        {
            const int n = xi.size();
            const Eigen::Matrix3Xd &ovPoly = *(Eigen::Matrix3Xd *)ptr;

            const double sqrNormXi = xi.squaredNorm();
            const double invNormXi = 1.0 / sqrt(sqrNormXi);
            const Eigen::VectorXd unitXi = xi * invNormXi;
            const Eigen::VectorXd r = unitXi.head(n - 1);
            const Eigen::Vector3d delta = ovPoly.rightCols(n - 1) * r.cwiseProduct(r) +
                                          ovPoly.col(1) - ovPoly.col(0);

            double cost = delta.squaredNorm();
            gradXi.head(n - 1) = (ovPoly.rightCols(n - 1).transpose() * (2 * delta)).array() *
                                 r.array() * 2.0;
            gradXi(n - 1) = 0.0;
            gradXi = (gradXi - unitXi.dot(gradXi) * unitXi).eval() * invNormXi;

            const double sqrNormViolation = sqrNormXi - 1.0;
            if (sqrNormViolation > 0.0)
            {
                double c = sqrNormViolation * sqrNormViolation;
                const double dc = 3.0 * c;
                c *= sqrNormViolation;
                cost += c;
                gradXi += dc * 2.0 * xi;
            }

            return cost;
        }

        template <typename EIGENVEC>
        static inline void backwardP(const Eigen::Matrix3Xd &P,
                                     const Eigen::VectorXi &vIdx,
                                     const PolyhedraV &vPolys,
                                     EIGENVEC &xi)
        {
            const int sizeP = P.cols();

            double minSqrD;
            lbfgs::lbfgs_parameter_t tiny_nls_params;
            tiny_nls_params.g_epsilon = FLT_EPSILON;
            tiny_nls_params.max_iterations = 128;

            Eigen::Matrix3Xd ovPoly;
            for (int i = 0, j = 0, k, l; i < sizeP; i++, j += k)
            {
                l = vIdx(i);
                k = vPolys[l].cols();
                xi.segment(j, k).setConstant(sqrt(1.0 / k));
                ovPoly.resize(3, k + 1);
                ovPoly.col(0) = P.col(i);
                ovPoly.rightCols(k) = vPolys[l];
                Eigen::VectorXd xii = xi.segment(j, k);
                lbfgs::lbfgs_optimize(xii,
                                      minSqrD,
                                      &GCOPTER::costTinyNLS,
                                      nullptr,
                                      nullptr,
                                      &ovPoly,
                                      tiny_nls_params);
                xi.segment(j, k) = xii;
            }

            return;
        }

        template <typename EIGENVEC>
        static inline void backwardGradP(const Eigen::VectorXd &xi,
                                         const Eigen::VectorXi &vIdx,
                                         const PolyhedraV &vPolys,
                                         const Eigen::Matrix3Xd &gradP,
                                         EIGENVEC &gradXi)
        {
            const int sizeP = vIdx.size();
            gradXi.resize(xi.size());

            double normInv;
            Eigen::VectorXd q, gradQ, unitQ;
            for (int i = 0, j = 0, k, l; i < sizeP; i++, j += k)
            {
                l = vIdx(i);
                k = vPolys[l].cols();
                q = xi.segment(j, k);
                normInv = 1.0 / q.norm();
                unitQ = q * normInv;
                gradQ.resize(k);
                gradQ.head(k - 1) = (vPolys[l].rightCols(k - 1).transpose() * gradP.col(i)).array() *
                                    unitQ.head(k - 1).array() * 2.0;
                gradQ(k - 1) = 0.0;
                gradXi.segment(j, k) = (gradQ - unitQ * unitQ.dot(gradQ)) * normInv;
            }

            return;
        }

        template <typename EIGENVEC>
        static inline void normRetrictionLayer(const Eigen::VectorXd &xi,
                                               const Eigen::VectorXi &vIdx,
                                               const PolyhedraV &vPolys,
                                               double &cost,
                                               EIGENVEC &gradXi)
        {
            const int sizeP = vIdx.size();
            gradXi.resize(xi.size());

            double sqrNormQ, sqrNormViolation, c, dc;
            Eigen::VectorXd q;
            for (int i = 0, j = 0, k; i < sizeP; i++, j += k)
            {
                k = vPolys[vIdx(i)].cols();

                q = xi.segment(j, k);
                sqrNormQ = q.squaredNorm();
                sqrNormViolation = sqrNormQ - 1.0;
                if (sqrNormViolation > 0.0)
                {
                    c = sqrNormViolation * sqrNormViolation;
                    dc = 3.0 * c;
                    c *= sqrNormViolation;
                    cost += c;
                    gradXi.segment(j, k) += dc * 2.0 * q;
                }
            }

            return;
        }

        static inline void setResVector(const Eigen::Vector3d &start,
                                        const Eigen::Vector3d &end,
                                        const Eigen::Matrix3Xd &points,
                                        const Eigen::VectorXd &times,
                                        const double integralRes,
                                        Eigen::VectorXi &resVector)
        {
            int N = times.size();
            double dis = sqrt(pow(points(0, 0) - start(0), 2) + pow(points(1, 0) - start(1), 2) + pow(points(2, 0) - start(2), 2));
            resVector(0) = (int)(dis * integralRes + 1);
            for (int i = 1; i < N - 1; i++)
            {
                dis = sqrt(pow(points(0, i) - points(0, i - 1), 2) + pow(points(1, i) - points(1, i - 1), 2) + pow(points(2, i) - points(2, i - 1), 2));
                resVector(i) = (int)(dis * integralRes + 1);
            }
            dis = sqrt(pow(points(0, N - 2) - end(0), 2) + pow(points(1, N - 2) - end(1), 2) + pow(points(2, N - 2) - end(2), 2));
            resVector(N - 1) = (int)(dis * integralRes + 1);
        }

        static inline bool smoothedL1(const double &x,
                                      const double &mu,
                                      double &f,
                                      double &df)
        {
            if (x < 0.0)
            {
                return false;
            }
            else if (x > mu)
            {
                f = x - 0.5 * mu;
                df = 1.0;
                return true;
            }
            else
            {
                const double xdmu = x / mu;
                const double sqrxdmu = xdmu * xdmu;
                const double mumxd2 = mu - 0.5 * x;
                f = mumxd2 * sqrxdmu * xdmu;
                df = sqrxdmu * ((-0.5) * xdmu + 3.0 * mumxd2 / mu);
                return true;
            }
        }

        static inline void cal_yaw(const Eigen::Vector3d &vel,
                                   const Eigen::Vector4d &quat,
                                   double &yaw_cal,
                                   double &yaw_dot_cal)
        {
            // calculate yaw
            Eigen::Vector3d zb, zb_norm, dir, dir_xb, xb, g;
            Eigen::Quaterniond q, q_, ori;
            Eigen::Matrix3d R;

            ori = Eigen::Quaterniond(quat(0), quat(1), quat(2), quat(3));
            R = ori.normalized().toRotationMatrix();
            xb = Eigen::Vector3d(R(0, 0), R(1, 0), R(2, 0));
            zb = Eigen::Vector3d(R(0, 2), R(1, 2), R(2, 2));
            zb_norm = zb.normalized();

            dir = vel;
            dir_xb = dir - dir.dot(zb_norm) * zb_norm;

            double theta = acos(xb.dot(dir_xb.normalized()));
            if (zb.dot(xb.cross(dir_xb)) > 0)
                yaw_cal = theta;
            else
                yaw_cal = -theta;
        }

        static inline void attachPenaltyFunctional(const Eigen::VectorXd &T,
                                                   const Eigen::MatrixX3d &coeffs,
                                                   const Eigen::VectorXi &hIdx,
                                                   const PolyhedraH &hPolys,
                                                   const double &smoothFactor,
                                                   const Eigen::VectorXi &resolutionVector,
                                                   const int &integralResolution,
                                                   const Eigen::VectorXd &magnitudeBounds,
                                                   const Eigen::VectorXd &penaltyWeights,
                                                   const bool &isStartOptOmgZ,
                                                   const bool &isStartOptFlip,
                                                   flatness::FlatnessMap &flatmap,
                                                   double &cost,
                                                   Eigen::VectorXd &gradT,
                                                   Eigen::MatrixX3d &gradC,
                                                   double &pos_cost,
                                                   double &vel_cost,
                                                   double &acc_cost,
                                                   double &omg_cost)
        {
            pos_cost = 0;
            vel_cost = 0;
            acc_cost = 0;
            omg_cost = 0;
            const double velSqrMax = magnitudeBounds(0) * magnitudeBounds(0);
            const double omgSqrMax = magnitudeBounds(1) * magnitudeBounds(1);
            const double thrustMean = 0.5 * (magnitudeBounds(2) + magnitudeBounds(3));
            const double thrustRadi = 0.5 * fabs(magnitudeBounds(2) - magnitudeBounds(3));
            const double thrustSqrRadi = thrustRadi * thrustRadi;

            const double weightPos = penaltyWeights(0);
            const double weightVel = penaltyWeights(1);
            const double weightOmg = penaltyWeights(2);
            const double weightThr = penaltyWeights(3);
            Eigen::Vector3d g(0, 0, 10);

            Eigen::Vector3d pos, vel, acc, jer, sna, dir, ddir;
            Eigen::Vector3d gradPos, gradVel, gradOmg;
            Eigen::Vector3d totalgradPos, totalgradVel, totalgradAcc, totalgradJer;
            Eigen::Vector3d totalgradDir, totalgradDirD;
            double totalgradpsi, totalgraddpsi;
            double gradThr;

            double step, alpha, thr, yaw, yaw_dot, omg_z;
            double s1, s2, s3, s4, s5, s6, s7;
            Eigen::Matrix<double, TRAJ_ORDER + 1, 1> beta0, beta1, beta2, beta3, beta4;
            Eigen::Vector3d outerNormal, omg;
            Eigen::Vector4d quat;

            int K, L;
            double violaPos, violaVel, violaThr, violaOmg, violaOmgZ;
            double violaPosPenaD, violaVelPenaD, violaThrPenaD, violaOmgPenaD, violaOmgZPenaD;
            double violaPosPena, violaVelPena, violaThrPena, violaOmgPena, violaOmgZPena;
            double node, pena;

            const int pieceNum = T.size();
            double integralFrac;
            for (int i = 0; i < pieceNum; i++)
            {
                const Eigen::Matrix<double, TRAJ_ORDER + 1, 3> &c = coeffs.block<TRAJ_ORDER + 1, 3>(i * (TRAJ_ORDER + 1), 0);
                integralFrac = 1.0 / resolutionVector(i);
                step = T(i) * integralFrac;
                for (int j = 0; j <= resolutionVector(i); j++)
                {
                    s1 = j * step;
                    s2 = s1 * s1;
                    s3 = s2 * s1;
                    s4 = s2 * s2;
                    s5 = s4 * s1;
                    s6 = s3 * s3;
                    s7 = s6 * s1;
#if TRAJ_ORDER == 3
                    beta0(0) = 1.0, beta0(1) = s1, beta0(2) = s2, beta0(3) = s3;
                    beta1(0) = 0.0, beta1(1) = 1.0, beta1(2) = 2.0 * s1, beta1(3) = 3.0 * s2;
                    beta2(0) = 0.0, beta2(1) = 0.0, beta2(2) = 2.0, beta2(3) = 6.0 * s1;
                    beta3(0) = 0.0, beta3(1) = 0.0, beta3(2) = 0.0, beta3(3) = 6.0;
                    beta4(0) = 0.0, beta4(1) = 0.0, beta4(2) = 0.0, beta4(3) = 0.0;
#elif TRAJ_ORDER == 5
                    beta0(0) = 1.0, beta0(1) = s1, beta0(2) = s2, beta0(3) = s3, beta0(4) = s4, beta0(5) = s5;
                    beta1(0) = 0.0, beta1(1) = 1.0, beta1(2) = 2.0 * s1, beta1(3) = 3.0 * s2, beta1(4) = 4.0 * s3, beta1(5) = 5.0 * s4;
                    beta2(0) = 0.0, beta2(1) = 0.0, beta2(2) = 2.0, beta2(3) = 6.0 * s1, beta2(4) = 12.0 * s2, beta2(5) = 20.0 * s3;
                    beta3(0) = 0.0, beta3(1) = 0.0, beta3(2) = 0.0, beta3(3) = 6.0, beta3(4) = 24.0 * s1, beta3(5) = 60.0 * s2;
                    beta4(0) = 0.0, beta4(1) = 0.0, beta4(2) = 0.0, beta4(3) = 0.0, beta4(4) = 24.0, beta4(5) = 120.0 * s1;
#elif TRAJ_ORDER == 7
                    beta0(0) = 1.0, beta0(1) = s1, beta0(2) = s2, beta0(3) = s3, beta0(4) = s4, beta0(5) = s5, beta0(6) = s6, beta0(7) = s7;
                    beta1(0) = 0.0, beta1(1) = 1.0, beta1(2) = 2.0 * s1, beta1(3) = 3.0 * s2, beta1(4) = 4.0 * s3, beta1(5) = 5.0 * s4, beta1(6) = 6.0 * s5, beta1(7) = 7.0 * s6;
                    beta2(0) = 0.0, beta2(1) = 0.0, beta2(2) = 2.0, beta2(3) = 6.0 * s1, beta2(4) = 12.0 * s2, beta2(5) = 20.0 * s3, beta2(6) = 30.0 * s4, beta2(7) = 42.0 * s5;
                    beta3(0) = 0.0, beta3(1) = 0.0, beta3(2) = 0.0, beta3(3) = 6.0, beta3(4) = 24.0 * s1, beta3(5) = 60.0 * s2, beta3(6) = 120.0 * s3, beta3(7) = 210.0 * s4;
                    beta4(0) = 0.0, beta4(1) = 0.0, beta4(2) = 0.0, beta4(3) = 0.0, beta4(4) = 24.0, beta4(5) = 120.0 * s1, beta4(6) = 360.0 * s2, beta4(7) = 840.0 * s3;
#endif
                    pos = c.transpose() * beta0;
                    vel = c.transpose() * beta1;
                    acc = c.transpose() * beta2;
                    jer = c.transpose() * beta3;
                    sna = c.transpose() * beta4;
                    dir = vel;
                    ddir = acc;

                    if (isStartOptFlip)
                    {
                        if (vel.norm() < 1.0e-3)
                        {
                            continue;
                        }
                        flatmap.forward(vel, acc, jer, dir, ddir, yaw, yaw_dot, thr, quat, omg);
                        if (isnan(omg(0)) || isnan(omg(1)) || isnan(omg(2)) || isinf(omg(0)) || isinf(omg(1)) || isinf(omg(2)))
                        {
                            continue;
                        }
                    }
                    else
                    {
                        flatmap.forward(vel, acc, jer, 0, 0, thr, quat, omg);
                    }
                    omg_z = omg(2);
                    omg(2) = 0;
                    violaOmg = omg.squaredNorm() - omgSqrMax;
                    violaOmgZ = omg_z * omg_z - M_PI * M_PI;
                    violaVel = vel.squaredNorm() - velSqrMax;
                    violaThr = (thr - thrustMean) * (thr - thrustMean) - thrustSqrRadi;

                    gradPos.setZero(), gradVel.setZero(), gradOmg.setZero();
                    pena = 0.0, gradThr = 0.0;
                    node = (j == 0 || j == resolutionVector(i)) ? 0.5 : 1.0;

                    L = hIdx(i);
                    K = hPolys[L].cols();
                    for (int k = 0; k < K; k++)
                    {
                        outerNormal = hPolys[L].col(k).head<3>();
                        violaPos = outerNormal.dot(pos - hPolys[L].col(k).tail<3>());
                        if (smoothedL1(violaPos, smoothFactor, violaPosPena, violaPosPenaD))
                        {
                            gradPos += weightPos * violaPosPenaD * outerNormal;
                            pos_cost += node * step * weightPos * violaPosPena;
                            pena += weightPos * violaPosPena;
                        }
                    }

                    if (smoothedL1(violaVel, smoothFactor, violaVelPena, violaVelPenaD))
                    {
                        gradVel += weightVel * violaVelPenaD * 2.0 * vel;
                        vel_cost += node * step * weightVel * violaVelPena;
                        pena += weightVel * violaVelPena;
                    }

                    if (smoothedL1(violaOmg, smoothFactor, violaOmgPena, violaOmgPenaD))
                    {
                        gradOmg += weightOmg * violaOmgPenaD * 2.0 * omg;
                        omg_cost += node * step * weightOmg * violaOmgPena;
                        pena += weightOmg * violaOmgPena;
                    }

                    if (isStartOptFlip && isStartOptOmgZ)
                    {
                        if (smoothedL1(violaOmgZ, smoothFactor, violaOmgZPena, violaOmgZPenaD))
                        {
                            gradOmg(2) += weightOmg * violaOmgZPenaD * 2.0 * omg_z;
                            omg_cost += node * step * weightOmg * violaOmgZPena;
                            pena += weightOmg * violaOmgZPena;
                        }
                    }

                    if (smoothedL1(violaThr, smoothFactor, violaThrPena, violaThrPenaD))
                    {
                        gradThr += weightThr * violaThrPenaD * 2.0 * (thr - thrustMean);
                        acc_cost += node * step * weightThr * violaThrPena;
                        pena += weightThr * violaThrPena;
                    }

                    if (isStartOptFlip)
                    {
                        flatmap.backward(dir, ddir, gradPos, gradVel, gradThr, Eigen::Vector4d::Zero(), gradOmg,
                                         totalgradPos, totalgradVel, totalgradAcc, totalgradJer, totalgradDir, totalgradDirD);
                        totalgradVel += totalgradDir;
                        totalgradAcc += totalgradDirD;
                    }
                    else
                    {
                        flatmap.backward(gradPos, gradVel, gradThr, Eigen::Vector4d::Zero(), gradOmg,
                                         totalgradPos, totalgradVel, totalgradAcc, totalgradJer, totalgradpsi, totalgraddpsi);
                    }

                    alpha = j * integralFrac;
                    gradC.block(i * (TRAJ_ORDER + 1), 0, TRAJ_ORDER + 1, 3) += (beta0 * totalgradPos.transpose() +
                                                                                beta1 * totalgradVel.transpose() +
                                                                                beta2 * totalgradAcc.transpose() +
                                                                                beta3 * totalgradJer.transpose()) *
                                                                               node * step;
                    gradT(i) += (totalgradPos.dot(vel) +
                                 totalgradVel.dot(acc) +
                                 totalgradAcc.dot(jer) +
                                 totalgradJer.dot(sna)) *
                                    alpha * node * step +
                                node * integralFrac * pena;
                    cost += node * step * pena;
                }
            }
            // std::cout << std::endl;
            return;
        }

        static inline void attachFlippenalty(const Eigen::VectorXd &T,
                                             const Eigen::MatrixX3d &coeffs,
                                             const Eigen::VectorXi &hIdx,
                                             const PolyhedraH &hPolys,
                                             const double &smoothFactor,
                                             const int &flipResolution,
                                             const Eigen::VectorXd &magnitudeBounds,
                                             const Eigen::VectorXd &penaltyWeights,
                                             const Eigen::VectorXd &usekeypos,
                                             const Eigen::Matrix3Xd &set_pos,
                                             const Eigen::Matrix3Xd &set_att,
                                             const Eigen::VectorXd &set_time,
                                             const bool &isOptSetTime,
                                             double &cost,
                                             Eigen::VectorXd &gradT,
                                             Eigen::MatrixX3d &gradC,
                                             Eigen::VectorXd &gradA,
                                             Eigen::VectorXd &flip_cost,
                                             Eigen::VectorXd &flip_pos_cost,
                                             Eigen::Matrix3Xd &flip_pos,
                                             Eigen::VectorXd &flip_err)
        {
            const int attN = set_time.size();
            if (attN == 0)
                return;

            flip_cost.setZero();
            flip_pos_cost.setZero();
            gradA.setZero();
            const double deltaT = magnitudeBounds(4);
            const double maxrange = magnitudeBounds(5);
            const double weightFlip = penaltyWeights(4);
            const double weightFlipPose = penaltyWeights(5);
            const double sampleN = flipResolution;
            const double Tsum = T.sum(), step = Tsum / sampleN * 2 * deltaT;
            const double integralSampleN = 1.0 / sampleN;
            const double rangeMax = maxrange * maxrange;
            const Eigen::Vector3d g(0, 0, 10);
            int trajIdx = 0;
            double sumHeadT = T[0];

            double violaFlip, violaFlipPena, violaFlipPenaD, pena, theta;
            double violaFlipPose, violaFlipPosePena, violaFlipPosePenaD;
            Eigen::Vector3d pos, vel, acc, jer, gradFlipC, gradFlipPose, att, tau, pxptau, setPose;
            double s1, s2, s3, s4, s5, s6, s7, nowTime, att_time, dot, gamma, graddotvel, graddotjer;
            Eigen::Matrix<double, TRAJ_ORDER + 1, 1> beta0, beta1, beta2, beta3;

            for (int j = 0; j < attN; j++)
            {
                att_time = Tsum * set_time(j);
                att = set_att.col(j).normalized();
                setPose = set_pos.col(j);
                sumHeadT = T[0];
                trajIdx = 0;
                for (int i = -sampleN / 2; i < sampleN / 2; i++)
                {
                    nowTime = att_time + i * step;
                    while (nowTime > sumHeadT)
                    {
                        trajIdx++;
                        if (trajIdx >= int(T.size()))
                        {
                            std::cout << "DeltaT = " << deltaT << " | T = " << T.transpose() << std::endl;
                            std::cout << "step = " << step << " | start_t = " << att_time + -sampleN / 2 * step << " | end_t = " << att_time + sampleN / 2 * step << std::endl;
                            std::cout << "trajIdx = " << trajIdx << "| time=" << nowTime << " | sumHeadT = " << sumHeadT << " | i = " << i << std::endl;
                        }
                        sumHeadT = T.head(trajIdx + 1).sum();
                    }

                    const Eigen::Matrix<double, TRAJ_ORDER + 1, 3> &c = coeffs.block<TRAJ_ORDER + 1, 3>(trajIdx * (TRAJ_ORDER + 1), 0);
                    s1 = nowTime - sumHeadT + T[trajIdx];
                    s2 = s1 * s1;
                    s3 = s2 * s1;
                    s4 = s2 * s2;
                    s5 = s4 * s1;
                    s6 = s3 * s3;
                    s7 = s6 * s1;
#if TRAJ_ORDER == 3
                    beta0(0) = 1.0, beta0(1) = s1, beta0(2) = s2, beta0(3) = s3;
                    beta1(0) = 0.0, beta1(1) = 1.0, beta1(2) = 2.0 * s1, beta1(3) = 3.0 * s2;
                    beta2(0) = 0.0, beta2(1) = 0.0, beta2(2) = 2.0, beta2(3) = 6.0 * s1;
                    beta3(0) = 0.0, beta3(1) = 0.0, beta3(2) = 0.0, beta3(3) = 6.0;
#elif TRAJ_ORDER == 5
                    beta0(0) = 1.0, beta0(1) = s1, beta0(2) = s2, beta0(3) = s3, beta0(4) = s4, beta0(5) = s5;
                    beta1(0) = 0.0, beta1(1) = 1.0, beta1(2) = 2.0 * s1, beta1(3) = 3.0 * s2, beta1(4) = 4.0 * s3, beta1(5) = 5.0 * s4;
                    beta2(0) = 0.0, beta2(1) = 0.0, beta2(2) = 2.0, beta2(3) = 6.0 * s1, beta2(4) = 12.0 * s2, beta2(5) = 20.0 * s3;
                    beta3(0) = 0.0, beta3(1) = 0.0, beta3(2) = 0.0, beta3(3) = 6.0, beta3(4) = 24.0 * s1, beta3(5) = 60.0 * s2;
#elif TRAJ_ORDER == 7
                    beta0(0) = 1.0, beta0(1) = s1, beta0(2) = s2, beta0(3) = s3, beta0(4) = s4, beta0(5) = s5, beta0(6) = s6, beta0(7) = s7;
                    beta1(0) = 0.0, beta1(1) = 1.0, beta1(2) = 2.0 * s1, beta1(3) = 3.0 * s2, beta1(4) = 4.0 * s3, beta1(5) = 5.0 * s4, beta1(6) = 6.0 * s5, beta1(7) = 7.0 * s6;
                    beta2(0) = 0.0, beta2(1) = 0.0, beta2(2) = 2.0, beta2(3) = 6.0 * s1, beta2(4) = 12.0 * s2, beta2(5) = 20.0 * s3, beta2(6) = 30.0 * s4, beta2(7) = 42.0 * s5;
                    beta3(0) = 0.0, beta3(1) = 0.0, beta3(2) = 0.0, beta3(3) = 6.0, beta3(4) = 24.0 * s1, beta3(5) = 60.0 * s2, beta3(6) = 120.0 * s3, beta3(7) = 210.0 * s4;
#endif
                    pos = c.transpose() * beta0;
                    vel = c.transpose() * beta1;
                    acc = c.transpose() * beta2;
                    jer = c.transpose() * beta3;

                    gradFlipC.setZero();
                    pena = 0.0;
                    tau = acc + g;

                    gamma = set_time(j) + 2.0 * deltaT * i / sampleN;
                    theta = 2 * i / sampleN * M_PI;
                    dot = tau.dot(att) / tau.norm();
                    violaFlip = std::cos(theta) - dot;
                    if (smoothedL1(violaFlip, smoothFactor, violaFlipPena, violaFlipPenaD))
                    {
                        pxptau = att / tau.norm() - att.dot(tau) * tau / pow(tau.norm(), 3);
                        gradFlipC -= weightFlip * violaFlipPenaD * pxptau;
                        flip_cost(j) += weightFlip * violaFlipPena * step;
                        pena = weightFlip * violaFlipPena;
                    }

                    if (abs(theta / M_PI * 180) < 0.01)
                    {
                        // for visualize
                        flip_pos(0, j) = pos(0);
                        flip_pos(1, j) = pos(1);
                        flip_pos(2, j) = pos(2);
                        violaFlipPose = (pos - setPose).squaredNorm() - rangeMax;
                        flip_err(j) = acos(dot) / M_PI * 180;

                        // set pose cost
                        if (usekeypos(j) && smoothedL1(violaFlipPose, smoothFactor, violaFlipPosePena, violaFlipPosePenaD))
                        {
                            gradFlipPose = 2 * weightFlipPose * violaFlipPosePenaD * (pos - setPose);
                            graddotvel = gradFlipPose.dot(vel);
                            // cal grad and cost
                            gradC.block(trajIdx * (TRAJ_ORDER + 1), 0, TRAJ_ORDER + 1, 3) += beta0 * gradFlipPose.transpose();
                            if (isOptSetTime)
                            {
                                gradA(j) += graddotvel * Tsum;
                            }
                            for (int k = 0; k < trajIdx; k++)
                            {
                                gradT(k) += graddotvel * (gamma - 1);
                            }
                            for (int k = trajIdx; k < T.size(); k++)
                            {
                                gradT(k) += graddotvel * gamma;
                            }
                            flip_pos_cost(j) += weightFlipPose * violaFlipPosePena;
                            cost += weightFlipPose * violaFlipPosePena;
                        }
                    }

                    graddotjer = gradFlipC.dot(jer);
                    if (isOptSetTime)
                    {
                        gradA(j) += graddotjer * step * Tsum;
                    }

                    gradC.block(trajIdx * (TRAJ_ORDER + 1), 0, TRAJ_ORDER + 1, 3) += beta2 * gradFlipC.transpose() * step;

                    for (int k = 0; k < trajIdx; k++)
                    {
                        gradT(k) += graddotjer * (gamma - 1) * step + pena * integralSampleN * 2 * deltaT;
                    }
                    for (int k = trajIdx; k < T.size(); k++)
                    {
                        gradT(k) += graddotjer * gamma * step + pena * integralSampleN * 2 * deltaT;
                    }

                    cost += pena * step;
                }
            }
            return;
        }

        static inline double costFunctional(void *ptr,
                                            const Eigen::VectorXd &x,
                                            Eigen::VectorXd &grad)
        {
            GCOPTER &obj = *(GCOPTER *)ptr;
            const int dimTau = obj.temporalDim;
            const int dimXi = obj.spatialDim;
            const int dimFlip = obj.flipDim;
            const double weightT = obj.rho;

            forwardT(x.segment(dimTau, dimFlip), obj.flipTimes);
            forwardFlipTime(obj.flipTimes, obj.deltaT, obj.set_time);
            obj.partialGradByFliptimes.setZero();
            forwardT(x.head(dimTau), obj.times);
            if (obj.times.maxCoeff() > 100 || obj.times.minCoeff() < 5.0e-4)
            {
                return DBL_MAX;
            }
            forwardP(x.tail(dimXi), obj.vPolyIdx, obj.vPolytopes, obj.points);

            double cost = 0;
            obj.minco.setParameters(obj.points, obj.times);
            obj.minco.getEnergy(cost);
            obj.minctrl_cost = cost;
            obj.minco.getEnergyPartialGradByCoeffs(obj.partialGradByCoeffs);
            obj.minco.getEnergyPartialGradByTimes(obj.partialGradByTimes);
            cost *= obj.penaltyWt(6);
            obj.minctrl_cost = cost;
            obj.partialGradByCoeffs *= obj.penaltyWt(6);
            obj.partialGradByTimes *= obj.penaltyWt(6);

            attachPenaltyFunctional(obj.times, obj.minco.getCoeffs(),
                                    obj.hPolyIdx, obj.hPolytopes,
                                    obj.smoothEps, obj.resVector, obj.integralRes,
                                    obj.magnitudeBd, obj.penaltyWt, obj.isStartOptOmgZ,
                                    obj.isStartOptFlip, obj.flatmap, cost,
                                    obj.partialGradByTimes, obj.partialGradByCoeffs,
                                    obj.pos_cost, obj.vel_cost, obj.acc_cost, obj.omg_cost);

            if (obj.isStartOptFlip)
            {
                attachFlippenalty(obj.times, obj.minco.getCoeffs(),
                                  obj.hPolyIdx, obj.hPolytopes, obj.smoothEps, obj.flipRes,
                                  obj.magnitudeBd, obj.penaltyWt, obj.usekeypos, obj.set_pos,
                                  obj.set_att, obj.set_time, obj.isoptsettime, cost,
                                  obj.partialGradByTimes, obj.partialGradByCoeffs,
                                  obj.partialGradByFliptimes, obj.flip_cost,
                                  obj.flip_pos_cost, obj.flip_pos, obj.flip_err);
            }

            obj.minco.propogateGrad(obj.partialGradByCoeffs, obj.partialGradByTimes,
                                    obj.gradByPoints, obj.gradByTimes);

            obj.backwardGradFlipTime(obj.set_time, obj.flipTimes, obj.partialGradByFliptimes,
                                     obj.deltaT, obj.gradByFliptimes);
                                     
            obj.time_cost = weightT * obj.times.sum();
            cost += obj.time_cost;
            obj.gradByTimes.array() += weightT;

            if (!obj.isStartOptAtt && obj.pos_cost + obj.vel_cost + obj.acc_cost + obj.omg_cost > obj.flip_cost.sum() + obj.flip_pos_cost.sum())
                obj.gradByFliptimes.setZero();
            else
                obj.isStartOptAtt = true;

            Eigen::VectorXd gradTau = grad.head(dimTau);
            Eigen::VectorXd gradFliptau = grad.segment(dimTau, dimFlip);
            Eigen::VectorXd gradXi = grad.tail(dimXi);
            backwardGradT(x.head(dimTau), obj.gradByTimes, gradTau);
            backwardGradT(x.segment(dimTau, dimFlip), obj.gradByFliptimes, gradFliptau);
            
            backwardGradP(x.tail(dimXi), obj.vPolyIdx, obj.vPolytopes, obj.gradByPoints, gradXi);
            normRetrictionLayer(x.tail(dimXi), obj.vPolyIdx, obj.vPolytopes, cost, gradXi);
            grad.head(dimTau) = gradTau;
            grad.segment(dimTau, dimFlip) = gradFliptau;
            grad.tail(dimXi) = gradXi;
            obj.total_cost = cost;
            return cost;
        }

        static inline int costProgress(void *ptr,
                                       const Eigen::VectorXd &x,
                                       const Eigen::VectorXd &grad,
                                       const double fx,
                                       const double step,
                                       const int k,
                                       const int ls)
        {
            GCOPTER &obj = *(GCOPTER *)ptr;
            obj.optCount = k;

            if (obj.isdebug)
            {
                // publish msg
                obj.optCost.iter = k;
                obj.optCost.line_search_count = ls;
                obj.optCost.gnorm = grad.norm();
                obj.optCost.total_cost = fx;
                obj.optCost.ctrl_cost = obj.minctrl_cost;
                obj.optCost.time_cost = obj.time_cost;
                obj.optCost.pos_cost = obj.pos_cost;
                obj.optCost.vel_cost = obj.vel_cost;
                obj.optCost.acc_cost = obj.acc_cost;
                obj.optCost.omg_cost = obj.omg_cost;
                obj.optCost.flip_cost_sum = obj.flip_cost.sum();
                obj.optCost.flip_pos_cost_sum = obj.flip_pos_cost.sum();
                for (int i = 0; i < obj.flipDim; i++)
                {
                    obj.optCost.flip_cost[i] = obj.flip_cost(i);
                    obj.optCost.flip_pos_cost[i] = obj.flip_pos_cost(i);
                }
                obj.optPub->publish(obj.optCost);
            }

            if (obj.isdebug && (k % 50 == 0 || k == 1))
            {
                const int dimTau = obj.temporalDim, dimXi = obj.spatialDim; //, dimFlip = obj.flipDim;
                forwardT(x.head(dimTau), obj.times);
                forwardP(x.tail(dimXi), obj.vPolyIdx, obj.vPolytopes, obj.points);

                obj.minco.getTrajectory(obj.Traj);
                obj.visualTraj->visualize(obj.Traj, 0);
                if (obj.isStartOptFlip)
                {
                    obj.visualTraj->visualizeEllipsoid(obj.Traj, int(obj.Traj.getTotalDuration() / 0.03));
                }
                obj.visualTraj->visualizeAttArrow(obj.set_att, obj.flip_pos);
                // usleep(200);
                // getchar(); // pause
            }

            // if (obj.isdebug && k % 50 == 0)
            // {
            //     std::cout << "===========================================================" << std::endl
            //               //   << "optT = " << x.head(dimTau).transpose() << std::endl
            //               //   << "times = " << obj.times.transpose() << std::endl
            //               //   << "points = " << std::endl << obj.points.transpose() << std::endl
            //               << "In the " << k << " iteration | startOptAtt = " << obj.isStartOptAtt << std::endl
            //               << "fx = " << fx << "\t gnorm = " << grad.norm() << std::endl
            //               << "ls = " << ls << "\t step = " << step << std::endl
            //               << "cost:\t" << obj.total_cost - obj.time_cost - obj.minctrl_cost << std::endl
            //               << "time:\t" << obj.time_cost << "\t ctrl:\t" << obj.minctrl_cost << std::endl
            //               << "pos:\t" << obj.pos_cost << "\t vel:\t" << obj.vel_cost << std::endl
            //               << "acc:\t" << obj.acc_cost << "\t omg:\t" << obj.omg_cost << std::endl
            //               << "flip:\t" << obj.flip_cost.sum() << " | deltaT = " << obj.deltaT << std::endl
            //               << "flip_cost:" << obj.flip_cost.transpose() << std::endl
            //               << "set_alpha = " << obj.set_time.transpose() * 100 << std::endl
            //               << "flip_err = " << obj.flip_err.transpose() << std::endl
            //               << "flip_pos_cost:\t" << obj.flip_pos_cost.sum() << std::endl
            //               << "flip_pos_cost:" << obj.flip_pos_cost.transpose() << std::endl;
            // }
            return 0;
        }

        static inline double costDistance(void *ptr,
                                          const Eigen::VectorXd &x,
                                          Eigen::VectorXd &grad)
        {
            void **dataPtrs = (void **)ptr;
            const double &dEps = *((const double *)(dataPtrs[0]));
            const Eigen::Vector3d &ini = *((const Eigen::Vector3d *)(dataPtrs[1]));
            const Eigen::Vector3d &fin = *((const Eigen::Vector3d *)(dataPtrs[2]));
            const PolyhedraV &vPolys = *((PolyhedraV *)(dataPtrs[3]));

            double cost = 0.0;
            const int overlaps = vPolys.size() / 2;

            Eigen::Matrix3Xd gradP = Eigen::Matrix3Xd::Zero(3, overlaps);
            Eigen::Vector3d a, b, d;
            Eigen::VectorXd r;
            double smoothedDistance;
            for (int i = 0, j = 0, k = 0; i <= overlaps; i++, j += k)
            {
                a = i == 0 ? ini : b;
                if (i < overlaps)
                {
                    k = vPolys[2 * i + 1].cols();
                    Eigen::Map<const Eigen::VectorXd> q(x.data() + j, k);
                    r = q.normalized().head(k - 1);
                    b = vPolys[2 * i + 1].rightCols(k - 1) * r.cwiseProduct(r) +
                        vPolys[2 * i + 1].col(0);
                }
                else
                {
                    b = fin;
                }

                d = b - a;
                smoothedDistance = sqrt(d.squaredNorm() + dEps);
                cost += smoothedDistance;

                if (i < overlaps)
                {
                    gradP.col(i) += d / smoothedDistance;
                }
                if (i > 0)
                {
                    gradP.col(i - 1) -= d / smoothedDistance;
                }
            }

            Eigen::VectorXd unitQ;
            double sqrNormQ, invNormQ, sqrNormViolation, c, dc;
            for (int i = 0, j = 0, k; i < overlaps; i++, j += k)
            {
                k = vPolys[2 * i + 1].cols();
                Eigen::Map<const Eigen::VectorXd> q(x.data() + j, k);
                Eigen::Map<Eigen::VectorXd> gradQ(grad.data() + j, k);
                sqrNormQ = q.squaredNorm();
                invNormQ = 1.0 / sqrt(sqrNormQ);
                unitQ = q * invNormQ;
                gradQ.head(k - 1) = (vPolys[2 * i + 1].rightCols(k - 1).transpose() * gradP.col(i)).array() *
                                    unitQ.head(k - 1).array() * 2.0;
                gradQ(k - 1) = 0.0;
                gradQ = (gradQ - unitQ * unitQ.dot(gradQ)).eval() * invNormQ;

                sqrNormViolation = sqrNormQ - 1.0;
                if (sqrNormViolation > 0.0)
                {
                    c = sqrNormViolation * sqrNormViolation;
                    dc = 3.0 * c;
                    c *= sqrNormViolation;
                    cost += c;
                    gradQ += dc * 2.0 * q;
                }
            }

            return cost;
        }

        static inline void getShortestPath(const Eigen::Vector3d &ini,
                                           const Eigen::Vector3d &fin,
                                           const PolyhedraV &vPolys,
                                           const double &smoothD,
                                           Eigen::Matrix3Xd &path)
        {
            const int overlaps = vPolys.size() / 2;

            Eigen::VectorXi vSizes(overlaps);
            for (int i = 0; i < overlaps; i++)
            {
                vSizes(i) = vPolys[2 * i + 1].cols();
            }

            Eigen::VectorXd xi(vSizes.sum());
            for (int i = 0, j = 0; i < overlaps; i++)
            {
                xi.segment(j, vSizes(i)).setConstant(sqrt(1.0 / vSizes(i)));
                j += vSizes(i);
            }

            double minDistance;
            void *dataPtrs[4];
            dataPtrs[0] = (void *)(&smoothD);
            dataPtrs[1] = (void *)(&ini);
            dataPtrs[2] = (void *)(&fin);
            dataPtrs[3] = (void *)(&vPolys);
            lbfgs::lbfgs_parameter_t shortest_path_params;
            shortest_path_params.past = 3;
            shortest_path_params.delta = 1.0e-3;

            lbfgs::lbfgs_optimize(xi,
                                  minDistance,
                                  &GCOPTER::costDistance,
                                  nullptr,
                                  nullptr,
                                  dataPtrs,
                                  shortest_path_params);

            path.resize(3, overlaps + 2);
            path.leftCols<1>() = ini;
            path.rightCols<1>() = fin;
            Eigen::VectorXd r;
            for (int i = 0, j = 0, k; i < overlaps; i++, j += k)
            {
                k = vPolys[2 * i + 1].cols();
                Eigen::Map<const Eigen::VectorXd> q(xi.data() + j, k);
                r = q.normalized().head(k - 1);

                path.col(i + 1) = vPolys[2 * i + 1].rightCols(k - 1) * r.cwiseProduct(r) +
                                  vPolys[2 * i + 1].col(0);
            }

            return;
        }

        static inline void sortCorridor(PolyhedronV &vpoly)
        {
            int minidx, n = vpoly.cols();
            double minpos, temppos;
            for (int i = 0; i < n; i++)
            {
                minidx = i;
                minpos = vpoly.col(i).sum();
                for (int j = i + 1; j < n; j++)
                {
                    temppos = vpoly.col(j).sum();
                    if (minpos > temppos)
                    {
                        minidx = j;
                        minpos = temppos;
                    }
                }
                vpoly.col(i).swap(vpoly.col(minidx));
            }
        }

        static inline bool processCorridor(const PolyhedraH &hPs,
                                           PolyhedraV &vPs)
        {
            const int sizeCorridor = hPs.size() - 1;

            vPs.clear();
            vPs.reserve(2 * sizeCorridor + 1);

            int nv;
            PolyhedronH curIH;
            PolyhedronV curIV, curIOB;
            for (int i = 0; i < sizeCorridor; i++)
            {

                if (!geoutils::enumerateVs(hPs[i], curIV))
                {
                    return false;
                }
                nv = curIV.cols();
                sortCorridor(curIV);
                curIOB.resize(3, nv);
                curIOB.col(0) = curIV.col(0);
                curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
                vPs.push_back(curIOB);

                curIH.resize(6, hPs[i].cols() + hPs[i + 1].cols());
                curIH.leftCols(hPs[i].cols()) = hPs[i];
                curIH.rightCols(hPs[i + 1].cols()) = hPs[i + 1];
                if (!geoutils::enumerateVs(curIH, curIV))
                {
                    return false;
                }
                nv = curIV.cols();
                sortCorridor(curIV);
                curIOB.resize(3, nv);
                curIOB.col(0) = curIV.col(0);
                curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
                vPs.push_back(curIOB);
            }

            if (!geoutils::enumerateVs(hPs.back(), curIV))
            {
                return false;
            }
            nv = curIV.cols();
            sortCorridor(curIV);
            curIOB.resize(3, nv);
            curIOB.col(0) = curIV.col(0);
            curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
            vPs.push_back(curIOB);

            return true;
        }

        static inline void setInitial(const Eigen::Matrix3Xd &path,
                                      const double &speed,
                                      const Eigen::VectorXi &intervalNs,
                                      Eigen::Matrix3Xd &innerPoints,
                                      Eigen::VectorXd &timeAlloc)
        {
            const int sizeM = intervalNs.size();
            const int sizeN = intervalNs.sum();
            innerPoints.resize(3, sizeN - 1);
            timeAlloc.resize(sizeN);

            Eigen::Vector3d a, b, c;
            for (int i = 0, j = 0, k = 0, l; i < sizeM; i++)
            {
                l = intervalNs(i);
                a = path.col(i);
                b = path.col(i + 1);
                c = (b - a) / l;
                timeAlloc.segment(j, l).setConstant(c.norm() / speed);
                j += l;
                for (int m = 0; m < l; m++)
                {
                    if (i > 0 || m > 0)
                    {
                        innerPoints.col(k++) = a + c * m;
                    }
                }
            }
            // timeAlloc(0) *= 2.0;
            // timeAlloc(timeAlloc.size() - 1) *= 2.0;
        }

    public:
        inline bool setup(const double &timeWeight,
                          const Eigen::Matrix<double, 3, 4> &initialPVAJ,
                          const Eigen::Matrix<double, 3, 4> &terminalPVAJ,
                          const PolyhedraH &safeCorridor,
                          const double &lengthPerPiece,
                          const double &smoothingFactor,
                          const int &integralResolution,
                          const int &flipResolution,
                          const Eigen::VectorXd &magnitudeBounds,
                          const Eigen::VectorXd &penaltyWeights,
                          const Eigen::VectorXd &useKeyPos,
                          const Eigen::Matrix3Xd &attSequence,
                          const Eigen::Matrix3Xd &attPose,
                          const Eigen::VectorXd &attTimeProportion,
                          const bool &isDebug,
                          const bool &isOptSetTime,
                          Visualizer &visual,
                          ros::Publisher &optPuber)
        {
            isStartOptFlip = false;
            isStartOptOmgZ = false;
            isStartOptAtt = false;
            isdebug = isDebug;
            isoptsettime = isOptSetTime;
            visualTraj = &visual;
            optPub = &optPuber;
            rho = timeWeight;
            headPVAJ = initialPVAJ;
            tailPVAJ = terminalPVAJ;

            usekeypos = useKeyPos;
            set_att = attSequence;
            set_pos = attPose;
            set_time = attTimeProportion;

            hPolytopes = safeCorridor;
            for (size_t i = 0; i < hPolytopes.size(); i++)
            { 
                hPolytopes[i].topRows<3>().colwise().normalize();
            }
            
            bool corridor_failed;
            corridor_failed = !processCorridor(hPolytopes, vPolytopes);
            if (corridor_failed)
            {

                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "                Corridor process failed!" << std::endl
                          << "      Press ENTER to check which intention is wrong!" << std::endl
                          << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                // for corridor debug, visualize all corridor both H and V reprensentation in order
                gcopter::GCOPTER::PolyhedraH tempH;
                for (int i = 0; i < int(vPolytopes.size()); i++)
                {
                    if (i % 2)
                    {
                        tempH.push_back(hPolytopes[i / 2]);
                        tempH.push_back(hPolytopes[i / 2 + 1]);
                    }
                    else
                    {
                        tempH.push_back(hPolytopes[i / 2]);
                    }
                    std::cout << std::endl
                              << vPolytopes[i].transpose() << std::endl;
                    visualTraj->visualizeVCorridorPoints(vPolytopes[i]);
                    visualTraj->visualizePolytope(tempH);
                    getchar();
                    tempH.clear();
                }
                visualTraj->visualizePolytope(hPolytopes);

                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "              Corridor visualization finished!" << std::endl
                          << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                return false;
            }

            polyN = hPolytopes.size();
            smoothEps = smoothingFactor;
            integralRes = integralResolution;
            flipRes = flipResolution;
            magnitudeBd = magnitudeBounds;
            penaltyWt = penaltyWeights;
            allocSpeed = magnitudeBd(0);
            deltaT = magnitudeBd(4);

            getShortestPath(headPVAJ.col(0), tailPVAJ.col(0),
                            vPolytopes, smoothEps, shortPath);

            const Eigen::Matrix3Xd deltas = shortPath.rightCols(polyN) - shortPath.leftCols(polyN);
            pieceIdx = (deltas.colwise().norm() / lengthPerPiece).cast<int>().transpose();
            pieceIdx.array() += 1;
            pieceN = pieceIdx.sum();

            temporalDim = pieceN;
            flipDim = set_time.size();
            optCost.flip_cost.resize(flipDim);
            optCost.flip_pos_cost.resize(flipDim);

            spatialDim = 0;
            vPolyIdx.resize(pieceN - 1);
            hPolyIdx.resize(pieceN);
            resVector.resize(pieceN);

            for (int i = 0, j = 0, k; i < polyN; i++)
            {
                k = pieceIdx(i);
                for (int l = 0; l < k; l++, j++)
                {                 
                    if (l < k - 1) 
                    {              
                        vPolyIdx(j) = 2 * i;
                        spatialDim += vPolytopes[2 * i].cols();
                    }
                    else if (i < polyN - 1) 
                    {
                        
                        vPolyIdx(j) = 2 * i + 1;
                        spatialDim += vPolytopes[2 * i + 1].cols();
                    }
                    hPolyIdx(j) = i;
                }
            }

#if TRAJ_ORDER == 3
            minco.setConditions(headPVAJ.leftCols(2), tailPVAJ.leftCols(2), pieceN);
#elif TRAJ_ORDER == 5
            minco.setConditions(headPVAJ.leftCols(3), tailPVAJ.leftCols(3), pieceN);
#elif TRAJ_ORDER == 7
            minco.setConditions(headPVAJ, tailPVAJ, pieceN);
#else
            return false;
#endif

            // Allocate temp variables
            points.resize(3, pieceN - 1);
            times.resize(pieceN);
            flipTimes.resize(flipDim);
            gradByPoints.resize(3, pieceN - 1);
            gradByTimes.resize(pieceN);
            gradByFliptimes.resize(flipDim);
            partialGradByCoeffs.resize((TRAJ_ORDER + 1) * pieceN, 3);
            partialGradByTimes.resize(pieceN);
            partialGradByFliptimes.resize(flipDim);
            flip_cost.resize(flipDim);
            flip_pos_cost.resize(flipDim);
            flip_pos.resize(3, flipDim);
            flip_err.resize(flipDim);
            flip_cost.setZero();
            flip_pos_cost.setZero();
            flip_err.setZero();

            return true;
        }

        inline double optimize(Trajectory<TRAJ_ORDER> &traj,
                               const double &relCostTol)
        {
            Traj = traj;
            Eigen::VectorXd x;
            x.resize(temporalDim + spatialDim + flipDim);
            x.setZero();

            setInitial(shortPath, allocSpeed, pieceIdx, points, times);

            // cal a new deltaT
            deltaT = deltaT / times.sum();
            magnitudeBd(4) = deltaT;

            // set sample point number
            setResVector(headPVAJ.col(0), tailPVAJ.col(0), points, times, integralRes, resVector);
            Eigen::VectorXd tau = x.head(temporalDim);
            Eigen::VectorXd fliptau = x.segment(temporalDim, flipDim);
            Eigen::VectorXd xi = x.tail(spatialDim);
            flipTimes = set_time;

            backwardT(times, tau);
            backwardT(flipTimes, fliptau);
            backwardP(points, vPolyIdx, vPolytopes, xi);
            x.head(temporalDim) = tau;
            x.segment(temporalDim, flipDim) = fliptau;
            x.tail(spatialDim) = xi;

            double minCostFunctional;
            lbfgs_params.mem_size = 32;
            lbfgs_params.min_step = 1.0e-32;
            lbfgs_params.g_epsilon = 1.0e-6;
            lbfgs_params.delta = relCostTol;
            lbfgs_params.max_linesearch = 256;

            int ret = lbfgs::lbfgs_optimize(x,
                                            minCostFunctional,
                                            &GCOPTER::costFunctional,
                                            nullptr,
                                            &GCOPTER::costProgress,
                                            this,
                                            lbfgs_params);

            totalOptCount[0] = optCount;
            forwardT(x.head(temporalDim), times);
            forwardP(x.tail(spatialDim), vPolyIdx, vPolytopes, points);
            setResVector(headPVAJ.col(0), tailPVAJ.col(0), points, times, integralRes, resVector);
            isStartOptFlip = true;

            ret = lbfgs::lbfgs_optimize(x,
                                        minCostFunctional,
                                        &GCOPTER::costFunctional,
                                        nullptr,
                                        &GCOPTER::costProgress,
                                        this,
                                        lbfgs_params);
            totalOptCount[1] = optCount;
            isStartOptOmgZ = true;
            forwardT(x.head(temporalDim), times);
            forwardP(x.tail(spatialDim), vPolyIdx, vPolytopes, points);
            setResVector(headPVAJ.col(0), tailPVAJ.col(0), points, times, integralRes, resVector);

            ret = lbfgs::lbfgs_optimize(x,
                                        minCostFunctional,
                                        &GCOPTER::costFunctional,
                                        nullptr,
                                        &GCOPTER::costProgress,
                                        this,
                                        lbfgs_params);

            totalOptCount[2] = optCount;
            forwardT(x.head(temporalDim), times);
            forwardP(x.tail(spatialDim), vPolyIdx, vPolytopes, points);
            minco.setParameters(points, times);
            minco.getTrajectory(traj);
            visualTraj->visualizeAttArrow(set_att, flip_pos);

            // check the opt result
            traj.resetMaxValues();
            Eigen::Vector3d pos, vel, acc, jer, sna, dir, ddir;
            Eigen::Vector3d omg, outerNormal;
            Eigen::Vector4d quat;
            double yaw, yaw_dot, thr, omg_xy, omg_z, dis, seg_t;
            int K, L;

            Eigen::VectorXd Time = traj.getDurations();
            double totalT = traj.getTotalDuration();
            const double checkResolution = 0.01;
            for (double t = 0; t < totalT; t += checkResolution)
            {
                seg_t = t;
                int idx = traj.locatePieceIdx(seg_t);
                pos = traj.getPos(t);
                vel = traj.getVel(t);
                acc = traj.getAcc(t);
                jer = traj.getJer(t);

                dir = vel;
                ddir = acc;

                if (vel.norm() < 1.0e-3)
                    continue;
                flatmap.forward(vel, acc, jer, dir, ddir, yaw, yaw_dot, thr, quat, omg);
                omg_z = abs(omg(2));
                omg(2) = 0;
                omg_xy = omg.norm();

                if (vel.norm() > traj.maxVel.norm())
                    traj.maxVel = vel;
                if (thr > traj.maxTau)
                    traj.maxTau = thr;
                if (thr < traj.minTau)
                    traj.minTau = thr;
                if (omg_z > traj.maxOmgz)
                    traj.maxOmgz = omg_z;
                if (omg_xy > traj.maxOmgxy)
                    traj.maxOmgxy = omg_xy;

                L = hPolyIdx(idx);;
                K = hPolytopes[L].cols();
                for (int k = 0; k < K; k++)
                {
                    outerNormal = hPolytopes[L].col(k).head<3>();
                    dis = outerNormal.dot(pos - hPolytopes[L].col(k).tail<3>());
                    dis = dis / outerNormal.norm();
                    if (dis > traj.maxDis)
                        traj.maxDis = dis;
                }
            }

            if (traj.maxVel.norm() > magnitudeBd(0) * 1.05){
                minCostFunctional = -101;
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "Max velocity: " << traj.maxVel.norm()
                          << "m/s, bound value: " << magnitudeBd(0) << "m/s"
                          << std::endl;
            } else if (traj.maxOmgxy > magnitudeBd(1) * 1.05){
                minCostFunctional = -102;
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "Max omega xy: " << traj.maxOmgxy
                          << "rad/s, bound value: " << magnitudeBd(1) << "rad/s."
                          << std::endl;
            } else if (traj.maxOmgz > magnitudeBd(1) * 1.05){
                minCostFunctional = -103;
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "Max omega z: " << traj.maxOmgz
                          << "rad/s, bound value: " << magnitudeBd(1) << "rad/s."
                          << std::endl;
            } else if (traj.maxTau > magnitudeBd(2) * 1.05){
                minCostFunctional = -104;
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "Max net thrust: " << traj.maxTau
                          << "m/s^2, bound value: " << magnitudeBd(2) << "m/s^2."
                          << std::endl;
            } else if (traj.minTau < magnitudeBd(3) * 0.95){
                minCostFunctional = -105;
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "Min net thrust: " << traj.minTau
                          << "m/s^2, bound value: " << magnitudeBd(3) << "m/s^2."
                          << std::endl;
            } else if (traj.maxDis > 0.1){
                minCostFunctional = -106;
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "Max distance out of corridor: " << traj.maxDis
                          << "m, bound value: 0.1m."
                          << std::endl;
            }

            if (ret < 0)
            {
                minCostFunctional = -100;
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "Optimization Failed. :(" << std::endl
                          << lbfgs::lbfgs_strerror(ret) << std::endl;
                return minCostFunctional;
            }
            else if (minCostFunctional < 0)
            {
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl
                          << "Safety check failed!" << std::endl
                          << "maxDis = " << traj.maxDis << " " << std::endl
                          << "maxVel = " << traj.maxVel.norm() << " " << std::endl
                          << "maxTau = " << traj.maxTau << " " << std::endl
                          << "minTau = " << traj.minTau << " " << std::endl
                          << "maxOxy = " << traj.maxOmgxy << " " << std::endl
                          << "maxOmz = " << traj.maxOmgz << std::endl
                          << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                return minCostFunctional;
            }
            else
            {
                std::cout << "Optimization Success!" << std::endl
                          << lbfgs::lbfgs_strerror(ret) << std::endl;
            }

            std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
                      << "total optimization count is " << totalOptCount.sum() << ", each is " << totalOptCount.transpose() << std::endl
                      << "cost:\t" << total_cost - time_cost - minctrl_cost << std::endl
                      << "time:\t" << time_cost << "\t ctrl:\t" << minctrl_cost << std::endl
                      << "pos:\t" << pos_cost << "\t vel:\t" << vel_cost << std::endl
                      << "acc:\t" << acc_cost << "\t omg:\t" << omg_cost << std::endl
                      << "flip:\t" << flip_cost.sum() << std::endl
                      << "flip_cost:" << flip_cost.transpose() << std::endl
                      << "set_alpha = " << set_time.transpose() * 100 << std::endl
                      << "flip_err = " << flip_err.transpose() << std::endl
                      << "flip_pos_cost:\t" << flip_pos_cost.sum() << std::endl
                      << "flip_pos_cost:" << flip_pos_cost.transpose() << std::endl;
                    //   << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;


            return minCostFunctional;
        }
    };

} // namespace gcopter

#endif
