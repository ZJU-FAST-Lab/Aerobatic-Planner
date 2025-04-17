#ifndef _ASTAR_HPP_
#define _ASTAR_HPP_

#include <vector>
#include <list>
#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/PoseArray.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/PointCloud2.h>

#include <pcl/io/pcd_io.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/kdtree.h>
#include <pcl/search/impl/kdtree.hpp>

#include "config.hpp"
#include "intention_get_corridor/glbmap.hpp"

class astar_ros;

class Node
{
public:
    Eigen::Vector3i idx;
    Node *parent;
    double F, G, H;
    // static int node_cnt;

    Node(Eigen::Vector3i index) : idx(index), parent(NULL), F(0),
                                  G(0), H(0){/*node_cnt++;*/};
    Node(int x, int y, int z) : idx(Eigen::Vector3i(x, y, z)), parent(NULL),
                                F(0), G(0), H(0){/*node_cnt++;*/};
    ~Node(){/*node_cnt--;*/};
    inline bool setParent(Node *parents)
    {
        if (parents != NULL)
        {
            this->parent = parents;
            return true;
        }
        else
        {
            std::cout << "empty parents" << std::endl;
            return false;
        }
    }
    inline void coutVel()
    {
        std::cout << "G = " << this->G << ", H = " << this->H << ", F = " << this->F << std::endl;
    }
};

typedef Node *NodePtr;

class Astar
{
public:
    std::vector<Eigen::Vector3i> pathList;
    std::vector<NodePtr> closeList, openList;

    Astar(){};
    ~Astar(){};
private:
    double h_wet;
    Eigen::Vector3i mapsize;
    bool (*mapCollisionCheck)(astar_ros *, Eigen::Vector3i);
    astar_ros *thisptr;
    Eigen::Vector3i startIdx, goalIdx;
    std::vector<NodePtr> neighborList;

public:
    inline void setHweight(double val){
        h_wet = val;
    }

    inline void setMapSize(const Eigen::Vector3i &mapSize)
    {
        this->mapsize = mapSize;
    }

    inline void InitAstar(const Eigen::Vector3i &mapSize, bool (*mapCheckCollision)(astar_ros *, Eigen::Vector3i), astar_ros *ptr)
    {
        this->mapsize = mapSize;
        if (ptr == NULL)
        {
            std::cout << "this ptr is null" << std::endl;
        }
        else
        {
            this->thisptr = ptr;
        }
        this->mapCollisionCheck = mapCheckCollision;
    }

    inline void pathfinding(const Eigen::Vector3i &start, const Eigen::Vector3i &goal)
    {
        AstarClean();
        // std::cout << "start = " << start.transpose() << std::endl;
        // std::cout << "goal = " << goal.transpose() << std::endl;
        if (checkIllegal(start) || checkIllegal(goal))
        {
            std::cout << "illegal index, please check given point" << std::endl;
            return;
        }

        if (checkCollision(start) || checkCollision(goal))
        {
            std::cout << "start or goal is in collision" << std::endl;
            return;
        }

        if (isEqual(start, goal))
        {
            std::cout << "start == goal" << std::endl;
            return;
        }

        this->startIdx = start;
        this->goalIdx = goal;

        NodePtr curPt = new Node(startIdx);
        NodePtr nebPt = NULL;
        std::vector<NodePtr>::iterator nodeiter = openList.end();
        curPt->G = 0;
        calH(curPt);
        calF(curPt);

        openList.push_back(curPt);
        int cnt = 0;
        while (openList.size() != 0)
        {
            cnt++;
            if (cnt > 20000)
            {
                std::cout << "A* stop because explore node is too much!" << std::endl;
                return;
            }
            nodeiter = findMinF();
            curPt = (*nodeiter);

            openList.erase(nodeiter);
            closeList.push_back(curPt);
            if (isEqual(curPt->idx, goalIdx))
            {
                // std::cout << "find path, total explore node cnt = " << cnt << std::endl;
                sourcePath(curPt);
                return;
            }
            getNeighbors(curPt);
            for (auto iter = neighborList.begin(); iter != neighborList.end(); iter++)
            {
                nebPt = (*iter);

                if ((nodeiter = isInList(nebPt, closeList)) != closeList.end())
                {
                    delete nebPt;
                }
                else if ((nodeiter = isInList(nebPt, openList)) != openList.end())
                {
                    delete nebPt;
                    nebPt = *nodeiter;
                    if (calG(nebPt, nebPt->parent) > calG(nebPt, curPt))
                    {
                        nebPt->setParent(curPt);
                        calG(nebPt);
                        calF(nebPt);
                    }
                }
                else
                {
                    nebPt->setParent(curPt);
                    calG(nebPt);
                    calH(nebPt);
                    calF(nebPt);
                    openList.push_back(nebPt);
                }
            }
        }
        std::cout << "no path" << std::endl;
    }

private:
    inline void AstarClean()
    {
        for (auto iter = openList.begin(); iter != openList.end(); iter++)
        {
            delete *iter;
        }
        // std::cout << "delete open" << std::endl;
        for (auto iter = closeList.begin(); iter != closeList.end(); iter++)
        {
            delete *iter;
        }
        // std::cout << "delete close" << std::endl;
        openList.clear();
        closeList.clear();
        neighborList.clear();
        pathList.clear();
    }

    inline void getNeighbors(NodePtr cur)
    {
        this->neighborList.clear();
        NodePtr nebPt = NULL;
        Eigen::Vector3i idx;
        // std::cout << "idx = " << cur->idx(0) << " " << cur->idx(1) << " " << cur->idx(2) << std::endl;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    if (i == 1 && j == 1 && k == 1)
                    {
                        continue;
                    }
                    else
                    {
                        idx(0) = cur->idx(0) + i - 1;
                        idx(1) = cur->idx(1) + j - 1;
                        idx(2) = cur->idx(2) + k - 1;
                        if (!checkCollision(idx))
                        {
                            nebPt = new Node(idx);
                            this->neighborList.push_back(nebPt);
                        }
                    }
                }
            }
        }
        // std::cout << neighborList.size() << std::endl;
    }

    inline std::vector<NodePtr>::iterator findMinF()
    {
        std::vector<NodePtr>::iterator minF = openList.begin();
        for (auto iter = openList.begin(); iter != openList.end(); iter++)
        {
            if ((*minF)->F > (*iter)->F)
            {
                minF = iter;
            }
        }
        return minF;
    }

    inline std::vector<NodePtr>::iterator isInList(NodePtr curPt, std::vector<NodePtr> &list)
    {
        std::vector<NodePtr>::iterator retIter = list.end();
        for (auto iter = list.begin(); iter != list.end(); iter++)
        {
            if (curPt->idx == (*iter)->idx)
            {
                retIter = iter;
                break;
            }
        }
        return retIter;
    }

    inline void sourcePath(NodePtr goal)
    {
        while (goal->parent != NULL)
        {
            pathList.push_back(goal->idx);
            goal = goal->parent;
        }
    }

    inline void calF(NodePtr cur)
    {
        cur->F = cur->G + (1 + h_wet) * cur->H;
    }

    inline void calG(NodePtr cur)
    {
        Eigen::Vector3d det;
        det[0] = double(abs(cur->idx[0] - cur->parent->idx[0]));
        det[1] = double(abs(cur->idx[1] - cur->parent->idx[1]));
        det[2] = double(abs(cur->idx[2] - cur->parent->idx[2]));

        cur->G = cur->parent->G + det.norm();
    }

    inline double calG(NodePtr cur, NodePtr par)
    {
        Eigen::Vector3d det;
        det[0] = double(abs(cur->idx[0] - par->idx[0]));
        det[1] = double(abs(cur->idx[1] - par->idx[1]));
        det[2] = double(abs(cur->idx[2] - par->idx[2]));

        return par->G + det.norm();
    }

    inline void calH(NodePtr cur)
    {
        Eigen::Vector3d det1 = Eigen::Vector3d((cur->idx - this->goalIdx)[0], (cur->idx - this->goalIdx)[1], (cur->idx - this->goalIdx)[2]);
        Eigen::Vector3d det2 = Eigen::Vector3d((cur->idx - this->startIdx)[0], (cur->idx - this->startIdx)[1], (cur->idx - this->startIdx)[2]);
        double cross = (det1.cross(det2)).norm();
        if (det1.norm() < 1e-5 || det2.norm() < 1e-5)
        {
            cross = 0;
        }
        else
        {
            cross /= (det1.norm() * det2.norm());
        }

        double temp;
        Eigen::Vector3d det = det1.cwiseAbs();
        if (det(0) < det(1))
        {
            temp = det(0);
            det(0) = det(1);
            det(1) = temp;
        }

        if (det(0) < det(2))
        {
            temp = det(0);
            det(0) = det(2);
            det(2) = temp;
        }

        if (det(1) < det(2))
        {
            temp = det(1);
            det(1) = det(2);
            det(2) = temp;
        }

        cur->H = (sqrt(3) - sqrt(2)) * det(2) + (sqrt(2) - 1) * det(1) + det(0);
        cur->H += cross;
    }

    inline bool checkCollision(int x, int y, int z)
    {
        Eigen::Vector3i idx(x, y, z);
        return this->checkCollision(idx);
    }

    inline bool checkCollision(Eigen::Vector3i idx)
    {
        if (checkIllegal(idx))
        {
            return true;
        }
        return (*(this->mapCollisionCheck))(thisptr, idx);
    }

    inline bool checkIllegal(Eigen::Vector3i idx)
    {
        for (int i = 0; i < 3; i++)
        {
            if (idx(i) < 0 || idx(i) > this->mapsize(i))
            {
                return true;
            }
        }
        return false;
    }

    inline bool isEqual(Eigen::Vector3i idx_1, Eigen::Vector3i idx_2)
    {
        for (int i = 0; i < 3; i++)
        {
            if (idx_1(i) != idx_2(i))
            {
                return false;
            }
        }
        return true;
    }
};

typedef Astar *astarPtr;

class astar_ros
{
private:
    ros::NodeHandle nh;
    ros::Publisher path_pub, list_pub;
    ros::Subscriber start_sub, goal_sub;
    astarPtr astar_finder;
    Eigen::Vector3i goal, start;
    sensor_msgs::PointCloud2 pathcloud, listcloud;
    pcl::PointCloud<pcl::PointXYZ> cloudMap;
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtreeMap;
    double map_size_x, map_size_y, map_size_z;
    double _resolution;
    double _x_l, _x_h, _y_l, _y_h, _z_l, _z_h;
    double _x_size, _y_size, _z_size;
    double roofheight, floorheight;
    bool have_map, hav_start, hav_init;

    inline void set_param(const Config &cfg)
    {
        astar_finder->setHweight(cfg.astar_weight);
        roofheight = cfg.expectedHeight[1] - 0.1;
        floorheight = cfg.expectedHeight[0] + 0.1;
        _resolution = cfg.gridResolution;
        _x_l = cfg.r3Bound[0];
        _x_h = cfg.r3Bound[1];
        _y_l = cfg.r3Bound[2];
        _y_h = cfg.r3Bound[3];
        _z_l = std::max(cfg.r3Bound[4], floorheight);
        _z_h = std::min(cfg.r3Bound[5], roofheight);

        _x_size = _x_h - _x_l;
        _y_size = _y_h - _y_l;
        _z_size = _z_h - _z_l;

        if (_x_size > 0 && _y_size > 0 && _z_size > 0 && _resolution > 0)
        {
            map_size_x = _x_size / _resolution;
            map_size_y = _y_size / _resolution;
            map_size_z = _z_size / _resolution;

            hav_init = true;
        } else {
            ROS_ERROR("map bound is wrong !");
        }
    }

    inline Eigen::Vector3i d2i(const Eigen::Vector3d pos)
    {
        Eigen::Vector3i idx;
        idx(0) = int((pos(0) - _x_l) / _resolution);
        idx(1) = int((pos(1) - _y_l) / _resolution);
        idx(2) = int((pos(2) - _z_l) / _resolution);
        return idx;   
    }

    inline Eigen::Vector3d i2d(const Eigen::Vector3i idx)
    {
        Eigen::Vector3d pos;
        pos[0] = (idx(0) + 0.5) * _resolution + _x_l;
        pos[1] = (idx(1) + 0.5) * _resolution + _y_l;
        pos[2] = (idx(2) + 0.5) * _resolution + _z_l;
        return pos;
    }

    inline void ros_init()
    {
        path_pub = nh.advertise<sensor_msgs::PointCloud2>("/astar/path", 1);
        list_pub = nh.advertise<sensor_msgs::PointCloud2>("/astar/exporationList", 1);
    }

    bool pubList(ros::Publisher &listpub, std::vector<Eigen::Vector3i> &astarList, sensor_msgs::PointCloud2 &pubCloud)
    {
        pcl::PointCloud<pcl::PointXYZ> pubpcl;
        pcl::PointXYZ pt;

        if (astarList.size() == 0)
        {
            std::cout << "list size = 0!!" << std::endl;
            return false;
        }
        else
        {
            for (auto iter = astarList.begin(); iter != astarList.end(); iter++)
            {
                Eigen::Vector3d pos = i2d(*iter);
                pt.x = pos(0);
                pt.y = pos(1);
                pt.z = pos(2);

                pubpcl.points.push_back(pt);
            }

            pubpcl.width = pubpcl.points.size();
            pubpcl.height = 1;
            pubpcl.is_dense = true;

            pcl::toROSMsg(pubpcl, pubCloud);
            pubCloud.header.frame_id = "world";
            listpub.publish(pubCloud);

            return true;
        }
    }

    void displayPath()
    {
        if (!pubList(path_pub, astar_finder->pathList, pathcloud))
        {
            return;
        }
        else
        {
            // std::cout << "find_path successfully" << std::endl;
            std::vector<Eigen::Vector3i> transList;

            for (auto iter = astar_finder->openList.begin(); iter != astar_finder->openList.end(); iter++)
            {
                transList.push_back((*iter)->idx);
            }

            for (auto iter = astar_finder->closeList.begin(); iter != astar_finder->closeList.end(); iter++)
            {
                transList.push_back((*iter)->idx);
            }

            pubList(list_pub, transList, listcloud);
        }
    }

    inline bool check_collision(Eigen::Vector3i idx)
    {
        std::vector<int> indices;
        std::vector<float> distances;

        Eigen::Vector3d pos = i2d(idx);
        pcl::PointXYZ pt(pos(0), pos(1), pos(2));

        int cnt = kdtreeMap.radiusSearch(pt, _resolution / 2, indices, distances);

        if (cnt)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

public:
    astar_ros(ros::NodeHandle &nh_, const Config &cfg)
        : nh(nh_), astar_finder(new Astar)
    {
        hav_start = false;
        have_map = false;
        hav_init = false;

        this->set_param(cfg);
        this->ros_init();

        astar_finder->InitAstar(Eigen::Vector3i(map_size_x, map_size_y, map_size_z), check_collision_t, this);
    };

    ~astar_ros()
    {
        if (this->astar_finder != NULL)
            delete this->astar_finder;
    };

    inline void setCloudMap(pcl::PointCloud<pcl::PointXYZ> &cloud)
    {
        cloudMap = cloud;
    }

    inline void change_map_bd(const Eigen::Vector3d bd_max, const Eigen::Vector3d bd_min)
    {
        _x_l = bd_min(0);
        _x_h = bd_max(0);
        _y_l = bd_min(1);
        _y_h = bd_max(1);
        _z_l = std::max(bd_min(2), floorheight);
        _z_h = std::min(bd_max(2), roofheight);

        _x_size = _x_h - _x_l;
        _y_size = _y_h - _y_l;
        _z_size = _z_h - _z_l;

        if (_x_size > 0 && _y_size > 0 && _z_size > 0 && _resolution > 0)
        {
            map_size_x = _x_size / _resolution;
            map_size_y = _y_size / _resolution;
            map_size_z = _z_size / _resolution;
            astar_finder->setMapSize(Eigen::Vector3i(map_size_x, map_size_y, map_size_z));

            hav_init = true;
        } else {
            ROS_ERROR("map bound is wrong !");
        }
    }

    static bool check_collision_t(astar_ros *th, Eigen::Vector3i idx)
    {
        return th->check_collision(idx);
    }

    inline void setstart(const nav_msgs::Odometry &msg)
    {
        Eigen::Vector3d pos;
        pos(0) = msg.pose.pose.position.x;
        pos(1) = msg.pose.pose.position.y;
        pos(2) = msg.pose.pose.position.z;
        start = d2i(pos);
        hav_start = true;
    }

    inline void setstart(const geometry_msgs::PoseStamped &msg)
    {
        Eigen::Vector3d pos;
        pos(0) = msg.pose.position.x;
        pos(1) = msg.pose.position.y;
        pos(2) = msg.pose.position.z;
        start = d2i(pos);
        hav_start = true;
    }

    inline void setstart(const Eigen::Vector3d pos)
    {
        start = d2i(pos);
        hav_start = true;
    }

    inline void setgoal(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        Eigen::Vector3d pos;
        pos(0) = msg->pose.position.x;
        pos(1) = msg->pose.position.y;
        pos(2) = std::max(msg->pose.position.z, floorheight);
        goal = d2i(pos);

        if (have_map == false)
        {
            kdtreeMap.setInputCloud(cloudMap.makeShared());
            have_map = true;
        }

        if (hav_start && have_map && hav_init)
        {
            astar_finder->pathfinding(start, goal);
            hav_start = false;
            displayPath();
        }
        else
        {
            ROS_WARN("Astar not start because lack of condition!");
            return;
        }
    }

    inline void setgoal(const Eigen::Vector3d pos)
    {
        goal = d2i(pos);

        if (have_map == false)
        {
            kdtreeMap.setInputCloud(cloudMap.makeShared());
            have_map = true;
        }

        if (hav_start && have_map && hav_init)
        {
            astar_finder->pathfinding(start, goal);
            hav_start = false;
            displayPath();
        }
        else
        {
            ROS_WARN("Astar not start because lack of condition!");
            return;
        }
    }

    inline void getpathlist(std::vector<Eigen::Vector3d> &list)
    {
        Eigen::Vector3d pt;
        // std::cout << "pathlist size is = " << astar_finder->pathList.size() << std::endl;
        for (auto iter = astar_finder->pathList.rbegin(); iter != astar_finder->pathList.rend(); iter ++)
        {
            pt = i2d(*iter);
            list.push_back(pt);
        }
        return;
    }
};

#endif