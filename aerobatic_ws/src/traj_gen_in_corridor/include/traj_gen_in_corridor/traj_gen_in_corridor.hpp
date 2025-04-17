#ifndef TRAJ_GEN_IN_CORRIDOR_HPP
#define TRAJ_GEN_IN_CORRIDOR_HPP

#include "traj_gen_in_corridor/config.hpp"
#include "traj_gen_in_corridor/visualizer.hpp"
#include "traj_gen_in_corridor/tictoc.hpp"
#include "gcopter/solver/geoutils.hpp"
#include "gcopter/trajectory.hpp"
#include "gcopter/gcopter.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <chrono>
#include <random>

#include <ros/ros.h>
#include <ros/console.h>
#include <ros/package.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PointStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/PoseArray.h>
#include <sensor_msgs/PointCloud2.h>
#include <nav_msgs/Odometry.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include "quadrotor_msgs/PolynomialTrajectory.h"
#include "quadrotor_msgs/Corridor.h"
#include "quadrotor_msgs/CorridorList.h"
#include "quadrotor_msgs/OptCostDebug.h"
#include "quadrotor_msgs/PositionCommand.h"
#include "quadrotor_msgs/TrajServerDebug.h"
#include "quadrotor_msgs/PointSeq.h"

class GlobalPlanner
{
private:
    Config config;

    ros::NodeHandle nh;
    ros::Subscriber targetSub, triggerSub;
    ros::Subscriber odomSub;
    ros::Subscriber corridorSub;
    ros::Subscriber keyposSub;
    ros::Subscriber intentionSub;
    ros::Publisher mapPub, optDebugPub;
    ros::Publisher trajPub, visTrajPub;

    bool corridorInitialized;
    bool odomInitialized;
    bool targetInitialized;
    bool postarInitialized;
    Eigen::Vector3d initialPos;
    Eigen::Vector3d finalPos;
    Visualizer visualizer;
    Eigen::Matrix3Xd key_pos;
    Eigen::Matrix3Xd key_att;
    Eigen::VectorXd key_time;
    Eigen::VectorXd useKeyPos;
    gcopter::GCOPTER gcopter;
    quadrotor_msgs::PolynomialTrajectory trajMsg;

    std::vector<Eigen::Matrix<double, 6, -1>> corridor;

public:
    GlobalPlanner(Config &conf, ros::NodeHandle &nh_)
        : config(conf), nh(nh_),
          corridorInitialized(false),
          odomInitialized(false),
          targetInitialized(false),
          postarInitialized(false),
          visualizer(config, nh)
    {
        trajPub = nh.advertise<quadrotor_msgs::PolynomialTrajectory>(config.trajTopic, 1);
        visTrajPub = nh.advertise<quadrotor_msgs::PolynomialTrajectory>(config.visTrajTopic, 1);
        optDebugPub = nh.advertise<quadrotor_msgs::OptCostDebug>("optCostChange", 1);
        odomSub = nh.subscribe(config.odomTopic, 1,
                               &GlobalPlanner::odomCallBack, this,
                               ros::TransportHints().tcpNoDelay());
        targetSub = nh.subscribe(config.targetTopic, 1,
                                 &GlobalPlanner::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
        corridorSub = nh.subscribe("corridor_list", 1,
                                   &GlobalPlanner::corridorCallback, this,
                                   ros::TransportHints().tcpNoDelay());
        triggerSub = nh.subscribe(config.triggerTopic, 1,
                                  &GlobalPlanner::triggerCallback, this,
                                  ros::TransportHints().tcpNoDelay());
        keyposSub = nh.subscribe("key_pos", 1,
                                 &GlobalPlanner::keyposCallback, this,
                                 ros::TransportHints().tcpNoDelay());
        intentionSub = nh.subscribe("/MyPointSeq", 1,
                                    &GlobalPlanner::IntentionCallback, this,
                                    ros::TransportHints().tcpNoDelay());
    }

    inline void IntentionCallback(const quadrotor_msgs::PointSeq::ConstPtr &msg)
    {
        if (!odomInitialized)
        {
            ROS_WARN("No Odom!");
            return;
        }
        int idx = msg->length - 1;
        finalPos(0) = msg->pos[idx * 3 + 0];
        finalPos(1) = msg->pos[idx * 3 + 1];
        finalPos(2) = std::max(msg->pos[idx * 3 + 2], 0.5);
    }

    inline void keyposCallback(const geometry_msgs::PoseArray::ConstPtr &msg)
    {
        int cnt_keypos = msg->poses.size();
        int cnt_keyatt = 0; // real cnt for acro opt
        Eigen::Matrix3Xd vis_pos;
        Eigen::Vector3d lth;
        vis_pos.resize(3, cnt_keypos);
        vis_pos.setZero();

        // save A* langth and topology pos
        auto iter = msg->poses.begin();
        for (int i = 0; i < cnt_keypos; i++)
        {
            lth(0) = iter->orientation.x;
            lth(1) = iter->orientation.y;
            lth(2) = iter->orientation.z;
            vis_pos(0, i) = iter->position.x;
            vis_pos(1, i) = iter->position.y;
            vis_pos(2, i) = iter->position.z;
            lth.normalize();
            if (lth(2) != 1.0)
                cnt_keyatt++;
            iter++;
        }

        // get key att, key pos and time
        double length = 0.0;
        useKeyPos.resize(cnt_keyatt);
        key_att.resize(3, cnt_keyatt);
        key_pos.resize(3, cnt_keyatt);
        key_time.resize(cnt_keyatt);
        iter = msg->poses.begin(); // idx for vis_pos
        int idx = 0;               // idx for key pos key att key time
        for (int i = 0; i < cnt_keypos; i++)
        {
            lth(0) = iter->orientation.x;
            lth(1) = iter->orientation.y;
            lth(2) = iter->orientation.z;
            length += lth.norm(); // A* path length between each vis_pos
            lth.normalize();
            if (lth(2) != 1.0) // all att except (0, 0, 1)
            {
                // key pos, may be not for opt
                key_pos(0, idx) = iter->position.x;
                key_pos(1, idx) = iter->position.y;
                key_pos(2, idx) = iter->position.z;
                // key att expect (0, 0, 1)
                key_att(0, idx) = lth(0);
                key_att(1, idx) = lth(1);
                key_att(2, idx) = lth(2);
                // weather use key pos as the key att constrain
                useKeyPos(idx) = iter->orientation.w;
                key_time(idx) = length;
                length = 0.0;
                idx++;
            }
            iter++;
        }

        key_time = key_time / length;
        std::cout << "key_time = " << key_time.transpose() << std::endl;
        std::cout << "key_pose = " << std::endl
                  << key_pos << std::endl;

        visualizer.visualizeSetPoints(vis_pos);
        postarInitialized = true;
    }

    inline void triggerCallback(const geometry_msgs::PointStamped::ConstPtr &msg)
    {
        std::cout << "Exec the trajectory !!" << std::endl;
        // wait for 3 seconds.
        int secs = 3;
        ros::Duration dr(1.0);
        while (ros::ok() && secs-- >= 0)
        {
            dr.sleep();
        }
        trajPub.publish(trajMsg);
    }

    inline void corridorCallback(const quadrotor_msgs::CorridorList::ConstPtr &msg)
    {
        corridor.clear();
        int totalNum = msg->corridor_cnt;
        for (int i = 0; i < totalNum; i++)
        {
            int size = msg->corridor_list[i].size;
            Eigen::Matrix<double, 6, -1> hP(6, size);
            if (msg->corridor_type == quadrotor_msgs::CorridorList::CORRIDOR_TYPE_H)
            {
                for (int j = 0; j < size; j++)
                {
                    hP.col(j)[0] = msg->corridor_list[i].nom_vec_list[j].x;
                    hP.col(j)[1] = msg->corridor_list[i].nom_vec_list[j].y;
                    hP.col(j)[2] = msg->corridor_list[i].nom_vec_list[j].z;

                    hP.col(j)[3] = msg->corridor_list[i].point_list[j].x;
                    hP.col(j)[4] = msg->corridor_list[i].point_list[j].y;
                    hP.col(j)[5] = msg->corridor_list[i].point_list[j].z;
                }
            }
            else
            {
                ROS_ERROR("Not ready for V reprensetion \n");
            }
            corridor.push_back(hP);
        }
        corridorInitialized = true;

        traj_generator();
        // exit(0);
    }

    inline void initializeCorridor(std::string path, int id)
    {
        std::ifstream fin(path);
        int totalNum = 0;
        fin >> totalNum;
        // std::cout << "totalNum : " << totalNum << std::endl;
        for (int i = 0; i < totalNum; i++)
        {
            int rowsize = 0, colsize = 0;
            fin >> rowsize;
            fin >> colsize;
            // std::cout << "curSize : " << rowsize << " " << colsize << std::endl;
            Eigen::Matrix<double, 6, -1> hP(6, colsize);
            double curVal;
            for (int j = 0; j < 6; j++)
            {
                for (int k = 0; k < colsize; k++)
                {
                    fin >> curVal;
                    hP(j, k) = curVal;
                }
            }
            if (id == 1)
            {
                corridor.push_back(hP);
            }
        }
        fin.close();
        corridorInitialized = true;
    }

    inline void odomCallBack(const nav_msgs::Odometry::ConstPtr &msg)
    {
        initialPos(0) = msg->pose.pose.position.x;
        initialPos(1) = msg->pose.pose.position.y;
        initialPos(2) = msg->pose.pose.position.z;

        odomInitialized = true;
    }

    inline void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        finalPos(0) = msg->pose.position.x;
        finalPos(1) = msg->pose.position.y;
        finalPos(2) = std::max(msg->pose.position.z, 0.5);
        targetInitialized = true;
    }

    inline void genPolyTrajMsg(const Trajectory<TRAJ_ORDER> &traj,
                               const Eigen::Isometry3d &tfR2L,
                               const ros::Time &iniStamp,
                               quadrotor_msgs::PolynomialTrajectory &trajMsg) const
    {
        trajMsg.header.stamp = iniStamp;
        static uint32_t traj_id = 0;
        traj_id++;
        trajMsg.trajectory_id = traj_id;
        trajMsg.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;
        trajMsg.num_order = traj[0].getDegree();
        trajMsg.num_segment = traj.getPieceNum();
        Eigen::Vector3d initialVel, finalVel;
        initialVel = tfR2L * traj.getVel(0.0);
        finalVel = tfR2L * traj.getVel(traj.getTotalDuration());
        trajMsg.start_yaw = 0.0;
        trajMsg.final_yaw = 0.0;

        for (size_t p = 0; p < (size_t)traj.getPieceNum(); p++)
        {
            trajMsg.time.push_back(traj[p].getDuration());
            trajMsg.order.push_back(traj[p].getCoeffMat().cols() - 1);

            Eigen::VectorXd linearTr(2);
            linearTr << 0.0, trajMsg.time[p];
            std::vector<Eigen::VectorXd> linearTrCoeffs;
            linearTrCoeffs.emplace_back(1);
            linearTrCoeffs[0] << 1;
            for (size_t k = 0; k < trajMsg.order[p]; k++)
            {
                linearTrCoeffs.push_back(RootFinder::polyConv(linearTrCoeffs[k], linearTr));
            }

            Eigen::MatrixXd coefMat(3, traj[p].getCoeffMat().cols());
            for (int i = 0; i < coefMat.cols(); i++)
            {
                coefMat.col(i) = tfR2L.rotation() * traj[p].getCoeffMat().col(coefMat.cols() - i - 1).head<3>();
            }
            coefMat.col(0) = (coefMat.col(0) + tfR2L.translation()).eval();

            for (int i = 0; i < coefMat.cols(); i++)
            {
                double coefx(0.0), coefy(0.0), coefz(0.0);
                for (int j = i; j < coefMat.cols(); j++)
                {
                    coefx += coefMat(0, j) * linearTrCoeffs[j](i);
                    coefy += coefMat(1, j) * linearTrCoeffs[j](i);
                    coefz += coefMat(2, j) * linearTrCoeffs[j](i);
                }
                trajMsg.coef_x.push_back(coefx);
                trajMsg.coef_y.push_back(coefy);
                trajMsg.coef_z.push_back(coefz);
            }
        }

        trajMsg.mag_coeff = 1.0;
        trajMsg.debug_info = "";
    }

    inline void traj_generator()
    {
        if (corridorInitialized && targetInitialized)
        {
            visualizer.visualizePolytope(corridor, 0);

            Eigen::Matrix<double, 3, 4> iniState;
            Eigen::Matrix<double, 3, 4> finState;

            iniState.setZero();
            finState.setZero();
            iniState(2, 0) = 1.5;
            finState(2, 0) = 1.5;

            if (odomInitialized)
            {
                iniState.col(0) = initialPos;
                finState.col(0) = finalPos;
            }

            double res = INFINITY;
            int itg = config.quadratureResolution;
            int smp = config.flipResolution;

            // initialize the trajectory planner
            Eigen::VectorXd magnitudeBounds(6);
            Eigen::VectorXd penaltyWeights(7);
            magnitudeBounds(0) = config.maxVelRate;
            magnitudeBounds(1) = config.maxOmgRate;
            magnitudeBounds(2) = config.maxTau;
            magnitudeBounds(3) = config.minTau;
            magnitudeBounds(4) = config.deltaT;
            magnitudeBounds(5) = config.maxRange;
            penaltyWeights(0) = config.chiVec[0];
            penaltyWeights(1) = config.chiVec[1];
            penaltyWeights(2) = config.chiVec[2];
            penaltyWeights(3) = config.chiVec[3];
            penaltyWeights(4) = config.chiVec[4];
            penaltyWeights(5) = config.chiVec[5];
            penaltyWeights(6) = config.chiVec[6];

            TicToc timer;
            Trajectory<TRAJ_ORDER> traj;

            // setup the trajectory planner with specified configuration
            if (!gcopter.setup(config.weightT, iniState, finState,
                               corridor, res, config.smoothingEps, itg, smp,
                               magnitudeBounds, penaltyWeights, useKeyPos,
                               key_att, key_pos, key_time, config.isDebug,
                               config.isOpt_settime, visualizer, optDebugPub))
            {
                std::cout << "WARNING: trajectory optimizer setup failed!" << std::endl;
                return;
            }

            // call the optimization considering collision avoidance with the other trajectory in gts with stamps in ginits
            int optFailedCount = 0;
            printf("\nOptimizing trajectory, please wait ...\n");
            while (gcopter.optimize(traj, config.relCostTol) < 0 && ros::ok() && optFailedCount < 10)
            {   // reset again and optimize again until success
                optFailedCount += 1;
                printf("Planning attempts is %d.\n\n", optFailedCount);
                printf("Attempint to replan, please be patient ...\n");
                gcopter.setup(config.weightT, iniState, finState,
                               corridor, res, config.smoothingEps, itg, smp,
                               magnitudeBounds, penaltyWeights, useKeyPos,
                               key_att, key_pos, key_time, config.isDebug,
                               config.isOpt_settime, visualizer, optDebugPub);
            }

            if (optFailedCount == 10){
                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                printf("          Unable to find a reasonable trajectory!\n");
                printf("  Please change intentions or constraints then try again!\n");
                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                return ;
            }

            double dtmsec = timer.toc();

            // Compute the trajectory length
            int segnum = 100000;
            double dt = traj.getTotalDuration() / segnum;
            Eigen::Vector3d xk, xk1 = traj.getPos(0.0);
            double length = 0.0;
            for (int i = 0; i < segnum; i++)
            {
                xk = xk1;
                xk1 = traj.getPos(dt * (i + 1.0));
                length += (xk1 - xk).norm();
            }

            // count the computation time
            printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
            printf("Spatial-temporal trajectory optimization called.\nTrajectory profile:\n");
            printf("       Total   computation  time: %7.2lf msecs.\n", dtmsec);
            printf("       Total   flight   duration: %7.2lf secs.\n", traj.getTotalDuration());
            printf("       Total  trajectory  length: %7.2lf m.\n", length);
            printf("       Average          velocity: %7.2lf m/s.\n", length / traj.getTotalDuration());
            printf("       Max  velocity   magnitude: %7.2lf m/s.\n", traj.getMaxVelRate());
            printf("       Max accleration magnitude: %7.2lf m/s^2.\n", traj.getMaxAccRate());
            printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");

            if (traj.getPieceNum() > 0)
            {
                visualizer.visualize(traj, 0, true);
                int samples = int(length / 0.1);
                visualizer.visualizeEllipsoid(traj, samples);
                quadrotor_msgs::PolynomialTrajectory vistrajMsg;
                genPolyTrajMsg(traj, Eigen::Isometry3d::Identity(), ros::Time::now(), vistrajMsg);
                int secs = 3;
                ros::Duration dr(1.0);
                while (ros::ok() && secs-- >= 0)
                {
                    dr.sleep();
                }
                trajMsg = vistrajMsg;
                visTrajPub.publish(vistrajMsg);
            }
        }
        return;
    }
};
#endif