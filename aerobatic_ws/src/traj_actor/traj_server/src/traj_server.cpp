#include <eigen3/Eigen/Dense>
// #include <pcl/point_cloud.h>
// #include <pcl/point_types.h>
// #include <pcl_conversions/pcl_conversions.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <tf/tf.h>
#include <fstream>
#include <tf/transform_datatypes.h>
#include <nav_msgs/Odometry.h>
#include <quadrotor_msgs/PolynomialTrajectory.h>
#include <quadrotor_msgs/PositionCommand.h>
#include <quadrotor_msgs/TrajServerDebug.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Point.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>
#include <sensor_msgs/PointCloud2.h>
#include "traj_server/flatness.hpp"

const int _DIM_x = 0;
const int _DIM_y = 1;
const int _DIM_z = 2;

using namespace std;

int _poly_order_min, _poly_order_max;

class TrajectoryServer
{
private:
    // Subscribers
    ros::Subscriber _odom_sub;
    ros::Subscriber _traj_sub;

    // publishers
    ros::Publisher _cmd_pub;
    ros::Publisher _vis_cmd_pub;
    ros::Publisher _vis_vel_pub;
    ros::Publisher _vis_acc_pub;
    ros::Publisher _vis_traj_pub;
    ros::Publisher _vis_traj_points;
    ros::Publisher _vis_thrust_pub;
    ros::Publisher _traj_debug_pub;
    bool use_yaw_comp;

    // configuration for trajectory
    double slow_speed = 1.0;
    double start_slow_time = 0.0;
    int _n_segment = 0;
    int _traj_id = 0;
    uint32_t _traj_flag = 0;
    Eigen::VectorXd _time;
    Eigen::MatrixXd _coef[3];
    vector<int> _order;

    double last_yaw, last_yaw_dot, time_forward;
    Eigen::Vector3d last_dir, last_dir_dot;
    flatness::FlatnessMap flatmap;

    double _vis_traj_width = 0.2;
    double mag_coeff;
    ros::Time _final_time = ros::TIME_MIN;
    ros::Time _start_time = ros::TIME_MAX;
    // double _start_yaw = 0.0, _final_yaw = 0.0;

    geometry_msgs::Point hover_position;

    // state of the server
    enum ServerState
    {
        INIT = 0,
        TRAJ,
        HOVER
    } state = INIT;
    ;
    nav_msgs::Odometry _odom;
    quadrotor_msgs::PositionCommand _cmd;
    geometry_msgs::PoseStamped _vis_cmd;

    visualization_msgs::Marker _vis_vel, _vis_acc, _vis_thrust;
    visualization_msgs::Marker _vis_traj;
    geometry_msgs::Point xb_pt, dir_xb_pt;
    quadrotor_msgs::TrajServerDebug _traj_debug;

    sensor_msgs::PointCloud2 traj_pts;
    // pcl::PointCloud<pcl::PointXYZ> traj_pts_pcd;
public:
    vector<Eigen::VectorXd> CList;  // Position coefficients vector, used to record all the pre-compute 'n choose k' combinatorial for the bernstein coefficients .
    vector<Eigen::VectorXd> CvList; // Velocity coefficients vector.
    vector<Eigen::VectorXd> CaList; // Acceleration coefficients vector.

    TrajectoryServer(ros::NodeHandle &handle)
    {
        handle.param("slow_speed", slow_speed, 1.0);
        handle.param("start_slow_time", start_slow_time, 0.0);
        handle.param("use_yaw_comp", use_yaw_comp, true);

        _odom_sub =
            handle.subscribe("odometry", 50, &TrajectoryServer::rcvOdometryCallback, this,
                             ros::TransportHints().tcpNoDelay());

        _traj_sub =
            handle.subscribe("trajectory", 2, &TrajectoryServer::rcvTrajectoryCallabck, this);

        _cmd_pub =
            handle.advertise<quadrotor_msgs::PositionCommand>("position_command", 50);

        _vis_cmd_pub =
            handle.advertise<geometry_msgs::PoseStamped>("desired_position", 50);

        _vis_vel_pub =
            handle.advertise<visualization_msgs::Marker>("desired_velocity", 50);

        _vis_acc_pub =
            handle.advertise<visualization_msgs::Marker>("desired_acceleration", 50);

        _vis_thrust_pub =
            handle.advertise<visualization_msgs::Marker>("desired_thrust", 50);

        _vis_traj_pub =
            handle.advertise<visualization_msgs::Marker>("trajectory_vis", 1);

        _traj_debug_pub =
            handle.advertise<quadrotor_msgs::TrajServerDebug>("trajserver_debug", 50);

        last_yaw = 0.0;
        last_yaw_dot = 0.0;
        time_forward = 0.1;
        last_dir = Eigen::Vector3d(1.0, 0.0, 0.0);
        last_dir_dot = Eigen::Vector3d::Zero();

        _vis_traj.header.stamp = ros::Time::now();
        _vis_traj.header.frame_id = "map";

        _vis_traj.ns = "trajectory/trajectory";
        _vis_traj.id = 0;
        _vis_traj.type = visualization_msgs::Marker::SPHERE_LIST;
        _vis_traj.action = visualization_msgs::Marker::ADD;
        _vis_traj.scale.x = _vis_traj_width;
        _vis_traj.scale.y = _vis_traj_width;
        _vis_traj.scale.z = _vis_traj_width;
        _vis_traj.pose.orientation.x = 0.0;
        _vis_traj.pose.orientation.y = 0.0;
        _vis_traj.pose.orientation.z = 0.0;
        _vis_traj.pose.orientation.w = 1.0;
        _vis_traj.color.r = 0.0;
        _vis_traj.color.g = 0.0;
        _vis_traj.color.b = 1.0;
        _vis_traj.color.a = 0.3;
        _vis_traj.points.clear();
    }

    bool cmd_flag = false;
    void rcvOdometryCallback(const nav_msgs::Odometry &odom)
    {
        // ROS_WARN("state = %d",state);
        static ros::Time last_odom_time = ros::Time::now();

        if (odom.child_frame_id == "X" || odom.child_frame_id == "O")
            return;
        // #1. store the odometry
        _odom = odom;
        _vis_cmd.header = _odom.header;
        _vis_cmd.header.frame_id = "world";

        if (state == INIT)
        {
            // ROS_WARN("[TRAJ SERVER] Pub initial pos command");
            _cmd.position = _odom.pose.pose.position;

            _cmd.header.stamp = _odom.header.stamp;
            _cmd.header.frame_id = "world";
            _cmd.trajectory_flag = _traj_flag;

            _cmd.velocity.x = 0.0;
            _cmd.velocity.y = 0.0;
            _cmd.velocity.z = 0.0;

            _cmd.acceleration.x = 0.0;
            _cmd.acceleration.y = 0.0;
            _cmd.acceleration.z = 0.0;
            // _cmd_pub.publish(_cmd);

            _cmd.yaw = last_yaw;
            _cmd.yaw_dot = 0.0;
            _cmd.yaw_dir.x = last_dir(0);
            _cmd.yaw_dir.y = last_dir(1);
            _cmd.yaw_dir.z = last_dir(2);
            _cmd.yaw_dir_dot.x = last_dir_dot(0);
            _cmd.yaw_dir_dot.y = last_dir_dot(1);
            _cmd.yaw_dir_dot.z = last_dir_dot(2);

            _cmd.jerk.x = 0.0;
            _cmd.jerk.y = 0.0;
            _cmd.jerk.z = 0.0;

            _vis_cmd.pose.position.x = _cmd.position.x;
            _vis_cmd.pose.position.y = _cmd.position.y;
            _vis_cmd.pose.position.z = _cmd.position.z;
            // _vis_cmd_pub.publish(_vis_cmd);

            return;
        }

        // change the order between #2 and #3. zxzxzxzx

        // #2. try to calculate the new state
        if (state == TRAJ && (((ros::Time::now() /*odom.header.stamp*/ - _start_time).toSec() / slow_speed / mag_coeff) > (_final_time - _start_time).toSec()))
        {
            state = HOVER;
            hover_position = _odom.pose.pose.position;
            _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_COMPLETED;
        }

        // #3. try to publish command
        pubPositionCommand(last_odom_time);
        last_odom_time = ros::Time::now();
    }

    void rcvTrajectoryCallabck(const quadrotor_msgs::PolynomialTrajectory &traj)
    {
        if (traj.action == quadrotor_msgs::PolynomialTrajectory::ACTION_ADD)
        {
            ROS_WARN("[traj_server] Loading the trajectory.");
            if ((int)traj.trajectory_id < 1)
            {
                ROS_ERROR("[traj_server] The trajectory_id must start from 1"); //. zxzxzxzx
                return;
            }
            if ((int)traj.trajectory_id > 1 && (int)traj.trajectory_id < _traj_id)
                return;

            state = TRAJ;
            _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_READY;
            _traj_id = traj.trajectory_id;
            _n_segment = traj.num_segment;
            _final_time = _start_time = ros::Time::now(); //_odom.header.stamp;
            _time.resize(_n_segment);

            _order.clear();
            for (int idx = 0; idx < _n_segment; ++idx)
            {
                _final_time += ros::Duration(traj.time[idx]);
                _time(idx) = traj.time[idx];
                _order.push_back(traj.order[idx]);
            }

            mag_coeff = traj.mag_coeff;

            int max_order = *max_element(begin(_order), end(_order));

            _coef[_DIM_x] = Eigen::MatrixXd::Zero(max_order + 1, _n_segment);
            _coef[_DIM_y] = Eigen::MatrixXd::Zero(max_order + 1, _n_segment);
            _coef[_DIM_z] = Eigen::MatrixXd::Zero(max_order + 1, _n_segment);

            // ROS_WARN("stack the coefficients");
            int shift = 0;
            for (int idx = 0; idx < _n_segment; ++idx)
            {
                int order = traj.order[idx];

                for (int j = 0; j < (order + 1); ++j)
                {
                    _coef[_DIM_x](j, idx) = traj.coef_x[shift + j];
                    _coef[_DIM_y](j, idx) = traj.coef_y[shift + j];
                    _coef[_DIM_z](j, idx) = traj.coef_z[shift + j];
                }

                shift += (order + 1);
            }
        }
        else if (traj.action == quadrotor_msgs::PolynomialTrajectory::ACTION_ABORT)
        {
            ROS_WARN("[SERVER] Aborting the trajectory.");
            state = HOVER;
            hover_position = _odom.pose.pose.position;
            _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_COMPLETED;
        }
        else if (traj.action == quadrotor_msgs::PolynomialTrajectory::ACTION_WARN_IMPOSSIBLE)
        {
            state = HOVER;
            hover_position = _odom.pose.pose.position;
            _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_IMPOSSIBLE;
        }
    }

    void pubPositionCommand(ros::Time last_odom_time)
    {
        if (state == INIT)
            return;
        if (state == HOVER)
            return;

        if (state == TRAJ)
        {
            _cmd.header.stamp = ros::Time::now(); //_odom.header.stamp;
            _cmd.header.frame_id = "world";
            _cmd.trajectory_flag = _traj_flag;
            _cmd.trajectory_id = _traj_id;

            double t = max(0.0, (_cmd.header.stamp - _start_time).toSec()); // / mag_coeff;
            // set when to start slow the speed
            if (t > start_slow_time)
            {
                t = (t - start_slow_time) / slow_speed + start_slow_time;
            }
            t = min(t, (_final_time - _start_time).toSec());

            double t_f = min(t + time_forward, _time.sum());
            bool use_last_yaw = abs(t_f - _time.sum()) < 1.0e-5 ? true : false;

            // ROS_WARN("[SERVER] the time : %.3lf\n, n = %d, m = %d", t, _n_order, _n_segment);
            for (int idx = 0; idx < _n_segment; ++idx)
            {
                if (t > _time[idx] && idx + 1 < _n_segment)
                {
                    t -= _time[idx];
                    t_f -= _time[idx];
                }
                else
                {
                    double f_idx = idx;
                    while (t_f > _time[f_idx] && f_idx + 1 < _n_segment)
                    {
                        t_f -= _time[f_idx];
                        f_idx++;
                    }
                    t /= _time[idx];
                    t_f /= _time[f_idx];
                    Eigen::Vector3d f_pos = Eigen::Vector3d::Zero();
                    Eigen::Vector3d f_vel = Eigen::Vector3d::Zero();

                    _cmd.position.x = 0.0;
                    _cmd.position.y = 0.0;
                    _cmd.position.z = 0.0;
                    _cmd.velocity.x = 0.0;
                    _cmd.velocity.y = 0.0;
                    _cmd.velocity.z = 0.0;
                    _cmd.acceleration.x = 0.0;
                    _cmd.acceleration.y = 0.0;
                    _cmd.acceleration.z = 0.0;
                    _cmd.jerk.x = 0.0;
                    _cmd.jerk.y = 0.0;
                    _cmd.jerk.z = 0.0;

                    int cur_order = _order[idx];
                    int cur_poly_num = cur_order + 1;

                    for (int i = 0; i < _order[f_idx] + 1; i++)
                    {
                        f_pos(0) += _coef[_DIM_x].col(f_idx)(i) * pow(t_f, i);
                        f_pos(1) += _coef[_DIM_y].col(f_idx)(i) * pow(t_f, i);
                        f_pos(2) += _coef[_DIM_z].col(f_idx)(i) * pow(t_f, i);
                        
                        if (i < (cur_poly_num - 1))
                        {
                            f_vel(0) += (i + 1) * _coef[_DIM_x].col(f_idx)(i + 1) * pow(t_f, i) / _time[f_idx];
                            f_vel(1) += (i + 1) * _coef[_DIM_y].col(f_idx)(i + 1) * pow(t_f, i) / _time[f_idx];
                            f_vel(2) += (i + 1) * _coef[_DIM_z].col(f_idx)(i + 1) * pow(t_f, i) / _time[f_idx];
                        }
                    }

                    for (int i = 0; i < cur_poly_num; i++)
                    {
                        _cmd.position.x += _coef[_DIM_x].col(idx)(i) * pow(t, i);
                        _cmd.position.y += _coef[_DIM_y].col(idx)(i) * pow(t, i);
                        _cmd.position.z += _coef[_DIM_z].col(idx)(i) * pow(t, i);

                        if (i < (cur_poly_num - 1))
                        {
                            _cmd.velocity.x += (i + 1) * _coef[_DIM_x].col(idx)(i + 1) * pow(t, i) / _time[idx];

                            _cmd.velocity.y += (i + 1) * _coef[_DIM_y].col(idx)(i + 1) * pow(t, i) / _time[idx];

                            _cmd.velocity.z += (i + 1) * _coef[_DIM_z].col(idx)(i + 1) * pow(t, i) / _time[idx];
                        }

                        if (i < (cur_poly_num - 2))
                        {
                            _cmd.acceleration.x += (i + 2) * (i + 1) * _coef[_DIM_x].col(idx)(i + 2) * pow(t, i) / _time[idx] / _time[idx];

                            _cmd.acceleration.y += (i + 2) * (i + 1) * _coef[_DIM_y].col(idx)(i + 2) * pow(t, i) / _time[idx] / _time[idx];

                            _cmd.acceleration.z += (i + 2) * (i + 1) * _coef[_DIM_z].col(idx)(i + 2) * pow(t, i) / _time[idx] / _time[idx];
                        }

                        if (i < (cur_poly_num - 3))
                        {
                            _cmd.jerk.x += (i + 3) * (i + 2) * (i + 1) * _coef[_DIM_x].col(idx)(i + 3) * pow(t, i) / _time[idx] / _time[idx] / _time[idx];

                            _cmd.jerk.y += (i + 3) * (i + 2) * (i + 1) * _coef[_DIM_y].col(idx)(i + 3) * pow(t, i) / _time[idx] / _time[idx] / _time[idx];

                            _cmd.jerk.z += (i + 3) * (i + 2) * (i + 1) * _coef[_DIM_z].col(idx)(i + 3) * pow(t, i) / _time[idx] / _time[idx] / _time[idx];
                        }
                    }

                    // calculate yaw
                    Eigen::Vector3d zb, zb_norm, xb;
                    Eigen::Quaterniond ori;
                    Eigen::Matrix3d R;
                    Eigen::Vector3d pos, vel, acc, jer, omg;
                    Eigen::Vector4d quat;
                    Eigen::Vector3d dir, dir_dot;
                    double thr, yaw, yaw_dot;

                    pos(0) = _cmd.position.x;
                    pos(1) = _cmd.position.y;
                    pos(2) = _cmd.position.z;
                    vel(0) = _cmd.velocity.x;
                    vel(1) = _cmd.velocity.y;
                    vel(2) = _cmd.velocity.z;
                    acc(0) = _cmd.acceleration.x;
                    acc(1) = _cmd.acceleration.y;
                    acc(2) = _cmd.acceleration.z;
                    jer(0) = _cmd.jerk.x;
                    jer(1) = _cmd.jerk.y;
                    jer(2) = _cmd.jerk.z;

                    dir = f_pos - pos;
                    dir_dot = f_vel - vel;
                    // dir = vel;
                    // dir_dot = acc;

                    // flatmap.calyawyawdot(vel, acc, jer, dir, dir_dot, yaw, yaw_dot);
                    if (use_last_yaw)
                    {
                        _cmd.yaw_dir.x = last_dir(0);
                        _cmd.yaw_dir.y = last_dir(1);
                        _cmd.yaw_dir.z = last_dir(2);
                        _cmd.yaw_dir_dot.x = last_dir_dot(0);
                        _cmd.yaw_dir_dot.y = last_dir_dot(1);
                        _cmd.yaw_dir_dot.z = last_dir_dot(2);
                    } else {
                        _cmd.yaw_dir.x = dir(0);
                        _cmd.yaw_dir.y = dir(1);
                        _cmd.yaw_dir.z = dir(2);
                        _cmd.yaw_dir_dot.x = dir_dot(0);
                        _cmd.yaw_dir_dot.y = dir_dot(1);
                        _cmd.yaw_dir_dot.z = dir_dot(2);
                        last_dir = dir;
                        last_dir_dot = dir_dot;
                    }

                    if (use_yaw_comp)
                    {
                        flatmap.forward(vel, acc, jer, dir, dir_dot, yaw, yaw_dot, thr, quat, omg);
                        _cmd.yaw = yaw;
                        _cmd.yaw_dot = yaw_dot;
                    } else {
                        flatmap.forward(vel, acc, jer, 0, 0, thr, quat, omg);
                        _cmd.yaw = 0.0;
                        _cmd.yaw_dot = 0.0;
                    }
                    ori = Eigen::Quaterniond(quat(0), quat(1), quat(2), quat(3));
                    R = ori.normalized().toRotationMatrix();
                    xb = Eigen::Vector3d(R(0, 0), R(1, 0), R(2, 0));
                    zb = Eigen::Vector3d(R(0, 2), R(1, 2), R(2, 2));
                    zb_norm = zb.normalized();

                    dir_xb_pt.x = xb(0);
                    dir_xb_pt.y = xb(1);
                    dir_xb_pt.z = xb(2);

                    // bound yaw into [-pi, pi]
                    while (_cmd.yaw > M_PI)
                        _cmd.yaw -= 2 * M_PI;

                    while (_cmd.yaw < -M_PI)
                        _cmd.yaw += 2 * M_PI;

                    // debug
                    _traj_debug.header = _cmd.header;
                    _traj_debug.pos = _cmd.position;
                    _traj_debug.vel = _cmd.velocity;
                    _traj_debug.acc = _cmd.acceleration;
                    _traj_debug.jer = _cmd.jerk;
                    _traj_debug.omg.x = omg(0);
                    _traj_debug.omg.y = omg(1);
                    _traj_debug.omg.z = omg(2);
                    _traj_debug.quat.w = quat(0);
                    _traj_debug.quat.x = quat(1);
                    _traj_debug.quat.y = quat(2);
                    _traj_debug.quat.z = quat(3);
                    _traj_debug.yaw = _cmd.yaw;
                    _traj_debug.yaw_dot = _cmd.yaw_dot;
                    _traj_debug.thr = thr;
                    _traj_debug.force.x = thr * zb_norm(0);
                    _traj_debug.force.y = thr * zb_norm(1);
                    _traj_debug.force.z = thr * zb_norm(2);
                    _traj_debug.acc_norm = acc.norm();
                    _traj_debug.jer_norm = jer.norm();
                    _traj_debug.vel_norm = vel.norm();
                    _traj_debug.omg_norm = sqrt(omg(0) * omg(0) + omg(1) * omg(1));

                    break;
                }
            }
        }

        _cmd_pub.publish(_cmd);
        _traj_debug_pub.publish(_traj_debug);

        _vis_cmd.header = _cmd.header;
        _vis_cmd.pose.position.x = _cmd.position.x;
        _vis_cmd.pose.position.y = _cmd.position.y;
        _vis_cmd.pose.position.z = _cmd.position.z;

        // tf::Quaternion q_ = tf::createQuaternionFromYaw(_cmd.yaw);
        // geometry_msgs::Quaternion odom_quat;
        // tf::quaternionTFToMsg(q_, odom_quat);
        // _vis_cmd.pose.orientation = odom_quat;
        _vis_cmd.pose.orientation.w = 1.0;
        _vis_cmd_pub.publish(_vis_cmd);

        _vis_vel.ns = "vel";
        _vis_vel.id = 0;
        _vis_vel.header.frame_id = "world";
        _vis_vel.type = visualization_msgs::Marker::ARROW;
        _vis_vel.action = visualization_msgs::Marker::ADD;
        _vis_vel.color.a = 1.0;
        _vis_vel.color.r = 0.0;
        _vis_vel.color.g = 1.0;
        _vis_vel.color.b = 0.0;

        _vis_vel.header.stamp = _odom.header.stamp;
        _vis_vel.points.clear();

        geometry_msgs::Point pt;
        pt.x = _cmd.position.x;
        pt.y = _cmd.position.y;
        pt.z = _cmd.position.z;

        _vis_traj.points.push_back(pt);
        _vis_traj_pub.publish(_vis_traj);

        _vis_vel.points.push_back(pt);

        pt.x = _cmd.position.x + _cmd.velocity.x / 8.0;
        pt.y = _cmd.position.y + _cmd.velocity.y / 8.0;
        pt.z = _cmd.position.z + _cmd.velocity.z / 8.0;

        _vis_vel.points.push_back(pt);

        _vis_vel.scale.x = 0.05;
        _vis_vel.scale.y = 0.1;
        _vis_vel.scale.z = 0.1;

        _vis_vel.pose.orientation.w = 1.0;
        _vis_vel_pub.publish(_vis_vel);

        _vis_acc.ns = "acc";
        _vis_acc.id = 0;
        _vis_acc.header.frame_id = "world";
        _vis_acc.type = visualization_msgs::Marker::ARROW;
        _vis_acc.action = visualization_msgs::Marker::ADD;
        _vis_acc.color.a = 1.0;
        _vis_acc.color.r = 1.0;
        _vis_acc.color.g = 0.0;
        _vis_acc.color.b = 0.0;

        _vis_acc.header.stamp = _odom.header.stamp;

        _vis_acc.points.clear();
        pt.x = _cmd.position.x;
        pt.y = _cmd.position.y;
        pt.z = _cmd.position.z;

        _vis_acc.points.push_back(pt);

        pt.x = _cmd.position.x + dir_xb_pt.x / 2.0;
        pt.y = _cmd.position.y + dir_xb_pt.y / 2.0;
        pt.z = _cmd.position.z + dir_xb_pt.z / 2.0;

        _vis_acc.points.push_back(pt);

        _vis_acc.scale.x = 0.05;
        _vis_acc.scale.y = 0.1;
        _vis_acc.scale.z = 0.1;

        _vis_acc.pose.orientation.w = 1.0;
        _vis_acc_pub.publish(_vis_acc);

        _vis_thrust.ns = "thrust";
        _vis_thrust.id = 0;
        _vis_thrust.header.frame_id = "world";
        _vis_thrust.type = visualization_msgs::Marker::ARROW;
        _vis_thrust.action = visualization_msgs::Marker::ADD;

        _vis_thrust.color.a = 1.0;
        _vis_thrust.color.r = 0.0;
        _vis_thrust.color.g = 0.0;
        _vis_thrust.color.b = 1.0;

        _vis_thrust.header.stamp = _odom.header.stamp;

        _vis_thrust.points.clear();
        pt.x = _cmd.position.x;
        pt.y = _cmd.position.y;
        pt.z = _cmd.position.z;

        _vis_thrust.points.push_back(pt);

        pt.x = _cmd.position.x + _cmd.acceleration.x / 12.0;
        pt.y = _cmd.position.y + _cmd.acceleration.y / 12.0;
        pt.z = _cmd.position.z + (_cmd.acceleration.z + 9.81) / 12.0;

        _vis_thrust.points.push_back(pt);

        _vis_thrust.scale.x = 0.05;
        _vis_thrust.scale.y = 0.1;
        _vis_thrust.scale.z = 0.1;

        _vis_thrust.pose.orientation.w = 1.0;
        _vis_thrust_pub.publish(_vis_thrust);
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "gradient_trajectory_server_node");
    ros::NodeHandle handle("~");

    TrajectoryServer server(handle);

    ros::spin();

    return 0;
}