#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include "intention_get_corridor/config.hpp"
#include "intention_get_corridor/solver/quickhull.hpp"

#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include <ros/ros.h>
#include <std_msgs/Float64.h>
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

    ros::Publisher meshPub;
    ros::Publisher edgePub;

public:
    Visualizer(Config &conf, ros::NodeHandle &nh_)
        : config(conf), nh(nh_)
    {
        meshPub = nh.advertise<visualization_msgs::Marker>("/visualizer/mesh", 1000);
        edgePub = nh.advertise<visualization_msgs::Marker>("/visualizer/edge", 1000);
    }

    // Visualize some polytopes in H-representation
    inline void visualizePolytope(const std::vector<Eigen::Matrix<double, 6, -1>> &hPolys,
                                  const int id = 0, const bool deleteAll = false)
    {
        // RVIZ support tris for visualization
        visualization_msgs::Marker meshMarker, edgeMarker;
        if (deleteAll){
            meshMarker.id = id;
            meshMarker.header.stamp = ros::Time::now();
            meshMarker.header.frame_id = "world";
            meshMarker.pose.orientation.w = 1.00;
            meshMarker.action = visualization_msgs::Marker::DELETE;
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
            meshPub.publish(meshMarker);
            edgePub.publish(edgeMarker);
        }
        if (hPolys.size() == 0){
            // ROS_WARN("visualizer corridor size is zero!!!\n");
            return;
        }
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
};

#endif