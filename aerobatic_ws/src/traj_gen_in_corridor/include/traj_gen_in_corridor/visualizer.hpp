#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include "traj_gen_in_corridor/config.hpp"
#include "gcopter/trajectory.hpp"
#include "gcopter/solver/quickhull.hpp"
#include "gcopter/solver/geoutils.hpp"

#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include <ros/ros.h>
#include <std_msgs/Float64.h>
#include <quadrotor_msgs/MINCOTraj.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

// Visualizer for the planner
class Visualizer
{
private:
    // config contains the scale for some markers
    Config config;
    ros::NodeHandle nh;

    // These are publishers for path, waypoints on the trajectory,
    // the entire trajectory, the mesh of free-space polytopes,
    // the edge of free-space polytopes, and spheres for safety radius
    ros::Publisher routePub;
    ros::Publisher wayPointsPub;
    ros::Publisher setPointsPub;
    ros::Publisher trajectoryPub;
    ros::Publisher meshPub;
    ros::Publisher edgePub;
    ros::Publisher spherePub;
    ros::Publisher timetextPub;
    ros::Publisher attPosPub;
    ros::Publisher ellipsoidPub;
    ros::Publisher unityvisPub;

public:
    ros::Publisher speedPub;
    ros::Publisher thrPub;
    ros::Publisher tiltPub;
    ros::Publisher bdrPub;
    ros::Publisher vCorridorPub;

public:
    Visualizer(Config &conf, ros::NodeHandle &nh_)
        : config(conf), nh(nh_)
    {
        routePub = nh.advertise<visualization_msgs::Marker>("/visualizer/route", 10);
        wayPointsPub = nh.advertise<visualization_msgs::Marker>("/visualizer/waypoints", 10);
        setPointsPub = nh.advertise<visualization_msgs::Marker>("/visualizer/setpoints", 10);
        trajectoryPub = nh.advertise<visualization_msgs::Marker>("/visualizer/trajectory", 10);
        meshPub = nh.advertise<visualization_msgs::Marker>("/visualizer/mesh", 1000);
        edgePub = nh.advertise<visualization_msgs::Marker>("/visualizer/edge", 1000);
        spherePub = nh.advertise<visualization_msgs::Marker>("/visualizer/spheres", 1000);
        timetextPub = nh.advertise<visualization_msgs::Marker>("/visualizer/timetext", 1000);
        speedPub = nh.advertise<std_msgs::Float64>("/visualizer/speed", 1000);
        thrPub = nh.advertise<std_msgs::Float64>("/visualizer/thr", 1000);
        tiltPub = nh.advertise<std_msgs::Float64>("/visualizer/tilt", 1000);
        bdrPub = nh.advertise<std_msgs::Float64>("/visualizer/bdr", 1000);
        attPosPub = nh.advertise<visualization_msgs::Marker>("/visualizer/attArrow", 100);
        vCorridorPub = nh.advertise<visualization_msgs::Marker>("/visualizer/VcorridorPoints", 1000);
        ellipsoidPub = nh.advertise<visualization_msgs::MarkerArray>("/visualizer/ellipsoid", 1000);
        unityvisPub = nh.advertise<quadrotor_msgs::MINCOTraj>("/visualizer/unitytrajvis", 1000);
    }

    // Visualize the trajectory and its front-end path
    inline void visualize(const Trajectory<TRAJ_ORDER> &traj,
                          int trajId = 0, 
                          bool unity_pub = false)
    {
        visualization_msgs::Marker routeMarker, wayPointsMarker, trajMarker, textMarker;
        quadrotor_msgs::MINCOTraj unityTraj;
        unityTraj.drone_id = 0;
        unityTraj.traj_id = trajId;
        unityTraj.order = traj.getPieceNum();
        unityTraj.start_v[0] = 0;
        unityTraj.start_v[1] = 0;
        unityTraj.start_v[2] = 0;
        unityTraj.start_a[0] = 0;
        unityTraj.start_a[1] = 0;
        unityTraj.start_a[2] = 0;
        unityTraj.end_v[0] = 0;
        unityTraj.end_v[1] = 0;
        unityTraj.end_v[2] = 0;
        unityTraj.end_a[0] = 0;
        unityTraj.end_a[1] = 0;
        unityTraj.end_a[2] = 0;

        routeMarker.id = trajId;
        routeMarker.type = visualization_msgs::Marker::LINE_LIST;
        routeMarker.header.stamp = ros::Time::now();
        routeMarker.header.frame_id = "world";
        routeMarker.pose.orientation.w = 1.00;
        routeMarker.action = visualization_msgs::Marker::ADD;
        routeMarker.ns = "route";
        routeMarker.color.r = 1.00;
        routeMarker.color.g = 0.00;
        routeMarker.color.b = 0.00;
        routeMarker.color.a = 1.00;
        routeMarker.scale.x = 0.1;

        wayPointsMarker = routeMarker;
        wayPointsMarker.id = -wayPointsMarker.id - 1;
        wayPointsMarker.type = visualization_msgs::Marker::SPHERE_LIST;
        wayPointsMarker.ns = "waypoints";
        wayPointsMarker.color.r = 1.00;
        wayPointsMarker.color.g = 0.00;
        wayPointsMarker.color.b = 0.00;
        wayPointsMarker.scale.x = 0.25;
        wayPointsMarker.scale.y = 0.25;
        wayPointsMarker.scale.z = 0.5;

        trajMarker = routeMarker;
        trajMarker.header.frame_id = "world";
        trajMarker.id = trajId;
        trajMarker.ns = "trajectory";
        trajMarker.color.r = 78 / 255.0;
        trajMarker.color.g = 154 / 255.0;
        trajMarker.color.b = 6 / 255.0;
        trajMarker.scale.x = 0.10;

        textMarker = routeMarker;
        textMarker.id = 0;
        textMarker.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
        textMarker.header.frame_id = "world";
        textMarker.ns = "text";
        textMarker.color.r = 0.00;
        textMarker.color.g = 0.50;
        textMarker.color.b = 1.00;
        textMarker.scale.z = 0.50;

        if (traj.getPieceNum() > 0)
        {
            Eigen::VectorXd times = traj.getDurations();
            Eigen::MatrixXd wps = traj.getPositions();
            for (int i = 0; i < wps.cols(); i++)
            {
                std::ostringstream out;
                geometry_msgs::Point point;
                point.x = wps.col(i)(0);
                point.y = wps.col(i)(1);
                point.z = wps.col(i)(2);
                wayPointsMarker.points.push_back(point);
                if (i < times.size())
                {
                    unityTraj.duration.push_back(times(i));
                    out << times(i);
                    textMarker.text = out.str();
                    textMarker.id = i;
                    textMarker.pose.position.x = point.x;
                    textMarker.pose.position.y = point.y;
                    textMarker.pose.position.z = point.z;
                    timetextPub.publish(textMarker);
                }

                if (i == 0)
                {
                    unityTraj.start_p[0] = wps.col(i)(0);
                    unityTraj.start_p[1] = wps.col(i)(1);
                    unityTraj.start_p[2] = wps.col(i)(2);
                }
                else if (i == wps.cols() - 1)
                {
                    unityTraj.end_p[0] = wps.col(i)(0);
                    unityTraj.end_p[1] = wps.col(i)(1);
                    unityTraj.end_p[2] = wps.col(i)(2);
                } else {
                    unityTraj.inner_x.push_back(wps.col(i)(0));
                    unityTraj.inner_y.push_back(wps.col(i)(1));
                    unityTraj.inner_z.push_back(wps.col(i)(2));
                }
            }
            if (unity_pub)
            {
                unityvisPub.publish(unityTraj);
            }
            wayPointsPub.publish(wayPointsMarker);
        }

        if (traj.getPieceNum() > 0)
        {
            double T = 0.01;
            Eigen::Vector3d lastX = traj.getPos(0.0);
            for (double t = T; t < traj.getTotalDuration(); t += T)
            {
                geometry_msgs::Point point;
                Eigen::Vector3d X = traj.getPos(t);
                point.x = lastX(0);
                point.y = lastX(1);
                point.z = lastX(2);
                trajMarker.points.push_back(point);
                point.x = X(0);
                point.y = X(1);
                point.z = X(2);
                trajMarker.points.push_back(point);
                lastX = X;
            }
            trajectoryPub.publish(trajMarker);
        }
    }

    // Visualize some polytopes in H-representation
    inline void visualizePolytope(const std::vector<Eigen::Matrix<double, 6, -1>> &hPolys,
                                  const int id = 0)
    {

        // Due to the fact that H-representation cannot be directly visualized
        // We first conduct vertex enumeration of them, then apply quickhull
        // to obtain triangle meshs of polyhedra
        Eigen::Matrix3Xd mesh(3, 0), curTris(3, 0), oldTris(3, 0);
        for (size_t id = 0; id < hPolys.size(); id++)
        {
            oldTris = mesh;
            Eigen::Matrix<double, 3, -1, Eigen::ColMajor> vPoly;
            geoutils::enumerateVs(hPolys[id], vPoly);

            quickhull::QuickHull<double> tinyQH;
            const auto polyHull = tinyQH.getConvexHull(vPoly.data(), vPoly.cols(), false, true);
            const auto &idxBuffer = polyHull.getIndexBuffer();
            int hNum = idxBuffer.size() / 3;

            curTris.resize(3, hNum * 3);
            for (int i = 0; i < hNum * 3; i++)
            {
                curTris.col(i) = vPoly.col(idxBuffer[i]);
            }
            mesh.resize(3, oldTris.cols() + curTris.cols());
            mesh.leftCols(oldTris.cols()) = oldTris;
            mesh.rightCols(curTris.cols()) = curTris;
        }

        // RVIZ support tris for visualization
        visualization_msgs::Marker meshMarker, edgeMarker;

        meshMarker.id = id;
        meshMarker.header.stamp = ros::Time::now();
        meshMarker.header.frame_id = "world";
        meshMarker.pose.orientation.w = 1.00;
        meshMarker.action = visualization_msgs::Marker::ADD;
        meshMarker.type = visualization_msgs::Marker::TRIANGLE_LIST;
        meshMarker.ns = "mesh";
        meshMarker.color.r = 0.00;
        meshMarker.color.g = 0.00;
        meshMarker.color.b = 1.00;
        meshMarker.color.a = 0.15;
        meshMarker.scale.x = 1.0;
        meshMarker.scale.y = 1.0;
        meshMarker.scale.z = 1.0;

        edgeMarker = meshMarker;
        edgeMarker.type = visualization_msgs::Marker::LINE_LIST;
        edgeMarker.ns = "edge";
        edgeMarker.color.r = 0.00;
        edgeMarker.color.g = 1.00;
        edgeMarker.color.b = 0.00;
        edgeMarker.color.a = 1.00;
        edgeMarker.scale.x = 0.02;

        geometry_msgs::Point point;

        int ptnum = mesh.cols();

        for (int i = 0; i < ptnum; i++)
        {
            point.x = mesh(0, i);
            point.y = mesh(1, i);
            point.z = mesh(2, i);
            meshMarker.points.push_back(point);
        }

        for (int i = 0; i < ptnum / 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                point.x = mesh(0, 3 * i + j);
                point.y = mesh(1, 3 * i + j);
                point.z = mesh(2, 3 * i + j);
                edgeMarker.points.push_back(point);
                point.x = mesh(0, 3 * i + (j + 1) % 3);
                point.y = mesh(1, 3 * i + (j + 1) % 3);
                point.z = mesh(2, 3 * i + (j + 1) % 3);
                edgeMarker.points.push_back(point);
            }
        }

        meshPub.publish(meshMarker);
        edgePub.publish(edgeMarker);

        return;
    }

    // Visualize all spheres with centers sphs and the same radius
    inline void visualizeSpheres(const std::vector<Eigen::Vector3d> &sphs,
                                 const double &radius)
    {
        visualization_msgs::Marker sphereMarkers, sphereDeleter;

        sphereMarkers.id = 0;
        sphereMarkers.type = visualization_msgs::Marker::SPHERE_LIST;
        sphereMarkers.header.stamp = ros::Time::now();
        sphereMarkers.header.frame_id = "world";
        sphereMarkers.pose.orientation.w = 1.00;
        sphereMarkers.action = visualization_msgs::Marker::ADD;
        sphereMarkers.ns = "spheres";
        sphereMarkers.color.r = 0.00;
        sphereMarkers.color.g = 0.00;
        sphereMarkers.color.b = 1.00;
        sphereMarkers.color.a = 1.00;
        sphereMarkers.scale.x = radius * 2.0;
        sphereMarkers.scale.y = radius * 2.0;
        sphereMarkers.scale.z = radius * 2.0;

        sphereDeleter = sphereMarkers;
        sphereDeleter.action = visualization_msgs::Marker::DELETE;

        geometry_msgs::Point point;
        for (const Eigen::Vector3d &ele : sphs)
        {
            point.x = ele(0);
            point.y = ele(1);
            point.z = ele(2);
            sphereMarkers.points.push_back(point);
        }

        spherePub.publish(sphereDeleter);
        spherePub.publish(sphereMarkers);
    }

    // Visualize all the key attitude, visualize their pos and direction
    inline void visualizeAttArrow(const Eigen::Matrix3Xd &set_att,
                                  const Eigen::Matrix3Xd &flip_pos)
    {
        visualization_msgs::Marker attArrow;
        geometry_msgs::Point pt;
        int attN = set_att.cols();

        attArrow.header.stamp = ros::Time::now();
        attArrow.header.frame_id = "world";
        attArrow.type = visualization_msgs::Marker::ARROW;
        attArrow.action = visualization_msgs::Marker::ADD;
        attArrow.color.a = 1.0;
        attArrow.color.r = 0.627450980;
        attArrow.color.g = 0.125490196;
        attArrow.color.b = 0.941176471;
        attArrow.scale.x = 0.1;
        attArrow.scale.y = 0.2;
        attArrow.scale.z = 0.2;
        attArrow.pose.orientation.w = 1.0;

        for (int i = 0; i < attN; i++)
        {
            attArrow.points.clear();
            attArrow.id = i;
            pt.x = flip_pos(0, i);
            pt.y = flip_pos(1, i);
            pt.z = flip_pos(2, i);
            attArrow.points.push_back(pt);
            pt.x += set_att(0, i) / 1.1;
            pt.y += set_att(1, i) / 1.1;
            pt.z += set_att(2, i) / 1.1;
            attArrow.points.push_back(pt);
            attPosPub.publish(attArrow);
        }
    }

    // Visualize all the setpoints
    inline void visualizeSetPoints(const Eigen::Matrix3Xd &setPoints)
    {
        visualization_msgs::Marker setPointsMarker, setPointDeleter;
        setPointsMarker.header.stamp = ros::Time::now();
        setPointsMarker.header.frame_id = "world";
        setPointsMarker.pose.orientation.w = 1.00;
        setPointsMarker.action = visualization_msgs::Marker::ADD;
        setPointsMarker.type = visualization_msgs::Marker::SPHERE_LIST;
        setPointsMarker.ns = "waypoints";
        setPointsMarker.color.r = 0.627450980;
        setPointsMarker.color.g = 0.125490196;
        setPointsMarker.color.b = 0.941176471;
        setPointsMarker.color.a = 1.00;
        setPointsMarker.scale.x = 0.25;
        setPointsMarker.scale.y = 0.25;
        setPointsMarker.scale.z = 0.5;

        setPointDeleter = setPointsMarker;
        setPointDeleter.action = visualization_msgs::Marker::DELETEALL;

        for (int i = 0; i < setPoints.cols(); i++)
        {
            setPointsMarker.id = i;
            geometry_msgs::Point point;
            point.x = setPoints(0, i);
            point.y = setPoints(1, i);
            point.z = setPoints(2, i);
            setPointsMarker.points.push_back(point);
        }
        setPointsPub.publish(setPointDeleter);
        setPointsPub.publish(setPointsMarker);
    }

    inline void visualizeVCorridorPoints(const Eigen::Matrix3Xd &vCorridor, const bool isvcorridor = true)
    {
        visualization_msgs::Marker setPointsMarker, setPointDeleter;
        setPointsMarker.header.stamp = ros::Time::now();
        setPointsMarker.header.frame_id = "world";
        setPointsMarker.pose.orientation.w = 1.00;
        setPointsMarker.action = visualization_msgs::Marker::ADD;
        setPointsMarker.type = visualization_msgs::Marker::SPHERE_LIST;
        setPointsMarker.ns = "waypoints";
        setPointsMarker.color.r = 0.5;
        setPointsMarker.color.g = 0.9;
        setPointsMarker.color.b = 0.5;
        setPointsMarker.color.a = 1.00;
        setPointsMarker.scale.x = 0.25;
        setPointsMarker.scale.y = 0.25;
        setPointsMarker.scale.z = 0.5;

        setPointDeleter = setPointsMarker;
        setPointDeleter.action = visualization_msgs::Marker::DELETEALL;
        Eigen::Vector3d start = vCorridor.col(0);
        for (int i = 0; i < vCorridor.cols(); i++)
        {
            setPointsMarker.id = i;
            geometry_msgs::Point point;
            point.x = vCorridor(0, i);
            point.y = vCorridor(1, i);
            point.z = vCorridor(2, i);
            if (isvcorridor && i)
            {
                point.x += start(0);
                point.y += start(1);
                point.z += start(2);
            }
            setPointsMarker.points.push_back(point);
        }
        vCorridorPub.publish(setPointDeleter);
        vCorridorPub.publish(setPointsMarker);
    }

    void visualizeEllipsoid(const Trajectory<TRAJ_ORDER> &traj, const int samples)
    {
        if (ellipsoidPub.getNumSubscribers() > 0)
        {
            visualization_msgs::Marker ellipsoidMarker;
            visualization_msgs::MarkerArray ellipsoidMarkers;

            ellipsoidMarker.id = 0;
            ellipsoidMarker.type = visualization_msgs::Marker::MESH_RESOURCE;
            ellipsoidMarker.mesh_resource = "package://traj_gen_in_corridor/data/sphere.stl";
            ellipsoidMarker.header.stamp = ros::Time::now();
            ellipsoidMarker.header.frame_id = "world";
            ellipsoidMarker.pose.orientation.w = 1.00;
            ellipsoidMarker.action = visualization_msgs::Marker::ADD;
            ellipsoidMarker.ns = "ellipsoids";
            ellipsoidMarker.color.r = 0.0;
            ellipsoidMarker.color.g = 0.0;
            ellipsoidMarker.color.b = 1.0;
            ellipsoidMarker.color.a = 0.2;
            ellipsoidMarker.scale.x = 0.60;
            ellipsoidMarker.scale.y = 0.60;
            ellipsoidMarker.scale.z = 0.10;
            ellipsoidMarker.pose.orientation.x = 0.0;
            ellipsoidMarker.pose.orientation.y = 0.0;
            ellipsoidMarker.pose.orientation.z = 0.0;
            ellipsoidMarker.pose.orientation.w = 1.0;
            ellipsoidMarker.mesh_use_embedded_materials = true;

            ellipsoidMarker.action = visualization_msgs::Marker::DELETEALL;
            ellipsoidMarkers.markers.push_back(ellipsoidMarker);
            ellipsoidPub.publish(ellipsoidMarkers);
            ellipsoidMarker.action = visualization_msgs::Marker::ADD;
            ellipsoidMarkers.markers.clear();

            double dt = traj.getTotalDuration() / samples;
            geometry_msgs::Point point;
            Eigen::Vector3d pos, thr;
            Eigen::Matrix3d rotM;
            Eigen::Quaterniond quat;
            for (int i = 0; i <= samples; i++)
            {
                pos = traj.getPos(dt * i);
                ellipsoidMarker.pose.position.x = pos(0);
                ellipsoidMarker.pose.position.y = pos(1);
                ellipsoidMarker.pose.position.z = pos(2);
                rotM.col(2) = traj.getAcc(dt * i);
                rotM(2, 2) += 9.81;
                rotM.col(2).normalize();
                rotM.col(1) = rotM.col(2).cross(Eigen::Vector3d(cos(0), sin(0), 0.0));
                rotM.col(1).normalize();
                rotM.col(0) = rotM.col(1).cross(rotM.col(2));
                quat = Eigen::Quaterniond(rotM);
                thr = rotM.col(2).normalized();
                ellipsoidMarker.color.r = 0.5 * (1.0 - thr(2)) * 1.0;
                ellipsoidMarker.color.g = 0.0;
                ellipsoidMarker.color.b = 0.5 * (1.0 + thr(2)) * 1.0;
                ellipsoidMarker.pose.orientation.w = quat.w();
                ellipsoidMarker.pose.orientation.x = quat.x();
                ellipsoidMarker.pose.orientation.y = quat.y();
                ellipsoidMarker.pose.orientation.z = quat.z();
                ellipsoidMarkers.markers.push_back(ellipsoidMarker);
                ellipsoidMarker.id++;
            }

            ellipsoidPub.publish(ellipsoidMarkers);
        }
    }
};

#endif