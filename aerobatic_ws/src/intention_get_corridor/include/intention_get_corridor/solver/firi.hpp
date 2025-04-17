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

#ifndef FIRI_HPP
#define FIRI_HPP

#include "lbfgs.hpp"
#include "sdlp.hpp"

#include <Eigen/Eigen>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <vector>

namespace firi
{

    inline Eigen::Vector3d origin2Triangle(const Eigen::Vector3d &a,
                                           const Eigen::Vector3d &b,
                                           const Eigen::Vector3d &c)
    {
        const Eigen::Vector3d ab = b - a;
        const Eigen::Vector3d ac = c - a;
        const double d1 = ab.dot(a);
        const double d2 = ac.dot(a);
        if (!(d1 < 0.0 || d2 < 0.0))
        {
            return a;
        }

        const double d3 = ab.dot(b);
        const double d4 = ac.dot(b);
        if (!(d3 > 0.0 || d4 < d3))
        {
            return b;
        }

        const double vc = d1 * d4 - d3 * d2;
        if (!(vc > 0.0 || d1 > 0.0 || d3 < 0.0))
        {
            return a + ab * d1 / (d1 - d3);
        }

        const double d5 = ab.dot(c);
        const double d6 = ac.dot(c);
        if (!(d6 > 0.0 || d5 < d6))
        {
            return c;
        }

        const double vb = d5 * d2 - d1 * d6;
        if (!(vb > 0.0 || d2 > 0.0 || d6 < 0.0))
        {
            return a + ac * d2 / (d2 - d6);
        }

        const double va = d3 * d6 - d5 * d4;
        if (!(va > 0.0 || (d4 - d3) > 0.0 || (d5 - d6) > 0.0))
        {
            return b + (c - b) * (d4 - d3) / ((d4 - d3) + (d5 - d6));
        }

        return a + (ab * vb + ac * vc) / (va + vb + vc);
    }

    inline void chol3d(const Eigen::Matrix3d &A,
                       Eigen::Matrix3d &L)
    {
        L(0, 0) = sqrt(A(0, 0));
        L(0, 1) = 0.0;
        L(0, 2) = 0.0;
        L(1, 0) = 0.5 * (A(0, 1) + A(1, 0)) / L(0, 0);
        L(1, 1) = sqrt(A(1, 1) - L(1, 0) * L(1, 0));
        L(1, 2) = 0.0;
        L(2, 0) = 0.5 * (A(0, 2) + A(2, 0)) / L(0, 0);
        L(2, 1) = (0.5 * (A(1, 2) + A(2, 1)) - L(2, 0) * L(1, 0)) / L(1, 1);
        L(2, 2) = sqrt(A(2, 2) - L(2, 0) * L(2, 0) - L(2, 1) * L(2, 1));
        return;
    }

    inline bool smoothedL1(const double &mu,
                           const double &x,
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

    inline double costMVIE(void *data,
                           const double *x,
                           double *grad,
                           const int n)
    {
        const int *pM = (int *)data;
        const double *pSmoothEps = (double *)(pM + 1);
        const double *pPenaltyWt = pSmoothEps + 1;
        const double *pA = pPenaltyWt + 1;

        const int M = *pM;
        const double smoothEps = *pSmoothEps;
        const double penaltyWt = *pPenaltyWt;
        Eigen::Map<const Eigen::MatrixX3d> A(pA, M, 3);
        Eigen::Map<const Eigen::Vector3d> p(x);
        Eigen::Map<const Eigen::Vector3d> rtd(x + 3);
        Eigen::Map<const Eigen::Vector3d> cde(x + 6);
        Eigen::Map<Eigen::Vector3d> gdp(grad);
        Eigen::Map<Eigen::Vector3d> gdrtd(grad + 3);
        Eigen::Map<Eigen::Vector3d> gdcde(grad + 6);

        double cost = 0;
        gdp.setZero();
        gdrtd.setZero();
        gdcde.setZero();

        Eigen::Matrix3d L;
        L(0, 0) = rtd(0) * rtd(0) + DBL_EPSILON;
        L(0, 1) = 0.0;
        L(0, 2) = 0.0;
        L(1, 0) = cde(0);
        L(1, 1) = rtd(1) * rtd(1) + DBL_EPSILON;
        L(1, 2) = 0.0;
        L(2, 0) = cde(2);
        L(2, 1) = cde(1);
        L(2, 2) = rtd(2) * rtd(2) + DBL_EPSILON;

        const Eigen::MatrixX3d AL = A * L;
        const Eigen::VectorXd normAL = AL.rowwise().norm();
        const Eigen::Matrix3Xd adjNormAL = (AL.array().colwise() / normAL.array()).transpose();
        const Eigen::VectorXd consViola = (normAL + A * p).array() - 1.0;

        double c, dc;
        Eigen::Vector3d vec;
        for (int i = 0; i < M; i++)
        {
            if (smoothedL1(smoothEps, consViola(i), c, dc))
            {
                cost += c;
                vec = dc * A.row(i).transpose();
                gdp += vec;
                gdrtd += adjNormAL.col(i).cwiseProduct(vec);
                gdcde(0) += adjNormAL(0, i) * vec(1);
                gdcde(1) += adjNormAL(1, i) * vec(2);
                gdcde(2) += adjNormAL(0, i) * vec(2);
            }
        }
        cost *= penaltyWt;
        gdp *= penaltyWt;
        gdrtd *= penaltyWt;
        gdcde *= penaltyWt;

        cost -= log(L(0, 0)) + log(L(1, 1)) + log(L(2, 2));
        gdrtd(0) -= 1.0 / L(0, 0);
        gdrtd(1) -= 1.0 / L(1, 1);
        gdrtd(2) -= 1.0 / L(2, 2);

        gdrtd(0) *= 2.0 * rtd(0);
        gdrtd(1) *= 2.0 * rtd(1);
        gdrtd(2) *= 2.0 * rtd(2);

        return cost;
    }

    // Each col of hPoly denotes a facet (outter_normal^T,point^T)^T
    // The outter_normal is assumed to be NORMALIZED
    // R, p, r are ALWAYS taken as the initial guess
    // R is also assumed to be a rotation matrix
    inline bool maxVolInsEllipsoid(const Eigen::Matrix<double, 6, -1> &hPoly,
                                   Eigen::Matrix3d &R,
                                   Eigen::Vector3d &p,
                                   Eigen::Vector3d &r)
    {
        // Find the deepest interior point
        int M = hPoly.cols();
        Eigen::MatrixXd Alp(M, 4);
        Eigen::VectorXd blp(M), clp(4), xlp(4);
        Alp.leftCols<3>() = hPoly.topRows<3>().transpose();
        Alp.rightCols<1>().setConstant(1.0);
        blp = hPoly.topRows<3>().cwiseProduct(hPoly.bottomRows<3>()).colwise().sum().transpose();
        clp.setZero();
        clp(3) = -1.0;
        double maxdepth = -sdlp::linprog(clp, Alp, blp, xlp);
        if (!(maxdepth > 0.0) || std::isinf(maxdepth))
        {
            return false;
        }
        Eigen::Vector3d interior = xlp.head<3>();

        // Prepare the data for MVIE optimization
        uint8_t *optData = new uint8_t[sizeof(int) + (2 + 3 * M) * sizeof(double)];
        int *pM = (int *)optData;
        double *pSmoothEps = (double *)(pM + 1);
        double *pPenaltyWt = pSmoothEps + 1;
        double *pA = pPenaltyWt + 1;

        *pM = M;
        Eigen::Map<Eigen::MatrixX3d> A(pA, M, 3);
        A = Alp.leftCols<3>().array().colwise() /
            (blp - Alp.leftCols<3>() * interior).array();

        Eigen::Matrix<double, 9, 1> x;
        const Eigen::Matrix3d Q = R * (r.cwiseProduct(r)).asDiagonal() * R.transpose();
        Eigen::Matrix3d L;
        chol3d(Q, L);

        x.head<3>() = p - interior;
        x(3) = sqrt(L(0, 0));
        x(4) = sqrt(L(1, 1));
        x(5) = sqrt(L(2, 2));
        x(6) = L(1, 0);
        x(7) = L(2, 1);
        x(8) = L(2, 0);

        double minCost;
        lbfgs::lbfgs_parameter_t paramsMVIE;
        lbfgs::lbfgs_load_default_parameters(&paramsMVIE);
        paramsMVIE.mem_size = 18;
        paramsMVIE.g_epsilon = 0.0;
        paramsMVIE.min_step = 1.0e-32;
        paramsMVIE.past = 3;
        paramsMVIE.delta = 1.0e-7;
        paramsMVIE.line_search_type = 0;
        *pSmoothEps = 1.0e-2;
        *pPenaltyWt = 1.0e+3;

        int ret = lbfgs::lbfgs_optimize(9,
                                        x.data(),
                                        &minCost,
                                        &costMVIE,
                                        nullptr,
                                        nullptr,
                                        optData,
                                        &paramsMVIE);

        if (ret < 0)
        {
            printf("WARNING: %s\n", lbfgs::lbfgs_strerror(ret));
        }

        p = x.head<3>() + interior;
        L(0, 0) = x(3) * x(3);
        L(0, 1) = 0.0;
        L(0, 2) = 0.0;
        L(1, 0) = x(6);
        L(1, 1) = x(4) * x(4);
        L(1, 2) = 0.0;
        L(2, 0) = x(8);
        L(2, 1) = x(7);
        L(2, 2) = x(5) * x(5);
        Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::FullPivHouseholderQRPreconditioner> svd(L, Eigen::ComputeFullU);
        const Eigen::Matrix3d U = svd.matrixU();
        const Eigen::Vector3d S = svd.singularValues();
        if (U.determinant() < 0.0)
        {
            R.col(0) = U.col(1);
            R.col(1) = U.col(0);
            R.col(2) = U.col(2);
            r(0) = S(1);
            r(1) = S(0);
            r(2) = S(2);
        }
        else
        {
            R = U;
            r = S;
        }

        delete[] optData;

        return ret >= 0;
    }

    inline void maximalVolInsPolytope(const Eigen::Matrix3Xd &mesh,
                                      const Eigen::Matrix3Xd &pc,
                                      const Eigen::Vector3d &seed,
                                      Eigen::Matrix<double, 6, -1> &hPolytope,
                                      const int iterations = 4,
                                      const double epsilon = 1.0e-6)
    {
        const int triangleNum = mesh.cols() / 3;
        const int pointNum = pc.cols();
        const int totalNum = triangleNum + pointNum;

        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
        Eigen::Vector3d p = seed;
        Eigen::Vector3d r = Eigen::Vector3d::Ones();
        std::vector<Eigen::Vector3d> normals;
        normals.reserve(totalNum);
        std::vector<Eigen::Vector3d> points;
        points.reserve(totalNum);
        for (int loop = 0; loop < iterations; loop++)
        {
            const Eigen::Matrix3d backward = R * r.asDiagonal();
            const Eigen::Matrix3d forward = r.cwiseInverse().asDiagonal() * R.transpose();
            const Eigen::Matrix3Xd forwardMesh = forward * (mesh.colwise() - p);
            const Eigen::Matrix3Xd forwardPC = forward * (pc.colwise() - p);

            Eigen::Matrix3Xd cPts(3, triangleNum);
            Eigen::VectorXd sqrRs(totalNum);
            std::vector<uint8_t> flags(totalNum, 1);
            for (int i = 0; i < triangleNum; i++)
            {
                cPts.col(i) = origin2Triangle(forwardMesh.col(3 * i),
                                              forwardMesh.col(3 * i + 1),
                                              forwardMesh.col(3 * i + 2));
            }
            sqrRs.head(triangleNum) = cPts.colwise().squaredNorm();
            sqrRs.tail(pointNum) = forwardPC.colwise().squaredNorm();

            normals.clear();
            points.clear();

            bool completed = false;
            int minId, cmpId;
            double minSqrR = sqrRs.minCoeff(&cmpId);
            for (int i = 0; !completed && i < totalNum; i++)
            {
                minId = cmpId;
                if (minId < triangleNum)
                {
                    points.push_back(cPts.col(minId));
                }
                else
                {
                    points.push_back(forwardPC.col(minId - triangleNum));
                }
                normals.push_back(points.back().normalized());
                const double b = normals.back().dot(points.back()) - epsilon;

                minSqrR = INFINITY;
                completed = true;
                flags[minId] = 0;
                for (int j = 0; j < triangleNum; j++)
                {
                    if (flags[j])
                    {
                        if (!(normals.back().dot(forwardMesh.col(3 * j)) < b ||
                              normals.back().dot(forwardMesh.col(3 * j + 1)) < b ||
                              normals.back().dot(forwardMesh.col(3 * j + 2)) < b))
                        {
                            flags[j] = 0;
                        }
                        else
                        {
                            completed = false;
                            if (minSqrR > sqrRs(j))
                            {
                                cmpId = j;
                                minSqrR = sqrRs(j);
                            }
                        }
                    }
                }
                for (int j = triangleNum; j < totalNum; j++)
                {
                    if (flags[j])
                    {
                        if (!(normals.back().dot(forwardPC.col(j - triangleNum)) < b))
                        {
                            flags[j] = 0;
                        }
                        else
                        {
                            completed = false;
                            if (minSqrR > sqrRs(j))
                            {
                                cmpId = j;
                                minSqrR = sqrRs(j);
                            }
                        }
                    }
                }
            }

            hPolytope.resize(6, normals.size());
            for (int i = 0; i < hPolytope.cols(); i++)
            {
                hPolytope.col(i).head<3>() = (forward.transpose() * normals[i]).normalized();
                hPolytope.col(i).tail<3>() = backward * points[i] + p;
            }

            if (loop == iterations - 1)
            {
                break;
            }

            maxVolInsEllipsoid(hPolytope, R, p, r);
        }

        return;
    }

}

#endif