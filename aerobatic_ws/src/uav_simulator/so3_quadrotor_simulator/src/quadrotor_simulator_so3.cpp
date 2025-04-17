#include <Eigen/Geometry>
#include <nav_msgs/Odometry.h>
#include <quadrotor_msgs/SO3Command.h>
#include <quadrotor_simulator/Quadrotor.h>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <uav_utils/geometry_utils.h>
#include <quadrotor_msgs/Motor.h>
#include <quadrotor_msgs/QuadStatus.h>

typedef struct _Control
{
  double rpm[4];
  double des_rpm[4];
  double real_rpm[4];
  double force;
  double M1, M2, M3;
  double in1, in2, in3;
  double eR1, eR2,eR3;
  double eOmg1, eOmg2, eOmg3;
} Control;

typedef struct _Command
{
  float force[3];
  float qx, qy, qz, qw;
  float omg[3];
  float kOmP[3], kOmI[3], kOmD[3];
  float corrections[3];
  float current_yaw;
  bool  use_external_yaw;
} Command;

typedef struct _Disturbance
{
  Eigen::Vector3d f;
  Eigen::Vector3d m;
} Disturbance;

static Command     command;
static Disturbance disturbance;

void stateToOdomMsg(const QuadrotorSimulator::Quadrotor::State& state,
                    nav_msgs::Odometry&                         odom);
void quadToImuMsg(const QuadrotorSimulator::Quadrotor& quad,
                  sensor_msgs::Imu&                    imu);

static Eigen::Vector3d 
BodyRatePIDCtrl(const Eigen::Vector3d &omg_sp,
                const Eigen::Vector3d &omg_rl,
                const Eigen::Vector3d &gain_P,
                const Eigen::Vector3d &gain_I,
                const Eigen::Vector3d &gain_D,
                const Eigen::Vector3d &gain_FF,
                const double &dt)
{
  static Eigen::Vector3d rate_int = Eigen::Vector3d::Zero(); 
  Eigen::Vector3d omg_err = omg_sp - omg_rl;
  Eigen::Vector3d torque;
  torque = gain_P.cwiseProduct(omg_err);
  torque -= gain_D.cwiseProduct(omg_err);
  torque += gain_FF.cwiseProduct(omg_sp);
  torque += rate_int;

  // update rate_int
  for (int i=0; i<3; i++)
  {
    double i_factor = omg_err(i) / (400.0 / 180.0 * M_PI);
    i_factor = std::max(0.0, 1.0-i_factor*i_factor);

    double rate_i = rate_int(i) + i_factor * gain_I(i) * omg_err(i) * dt;
    if (!isinf(rate_i))
    {
      rate_int(i) = std::max(std::min(rate_i, 1.0), -1.0);
    }
  }

  return torque;
}

static Control
getControl(const QuadrotorSimulator::Quadrotor& quad, const Command& cmd, const double &dt)
{
  const double _kf = quad.getPropellerThrustCoefficient();
  const double _km = quad.getPropellerMomentCoefficient();
  const double kf  = _kf - cmd.corrections[0];
  const double km  = _km / _kf * kf;

  const double          d       = quad.getArmLength();
  const Eigen::Matrix3f J       = quad.getInertia().cast<float>();
  const float           I[3][3] = { { J(0, 0), J(0, 1), J(0, 2) },
                                    { J(1, 0), J(1, 1), J(1, 2) },
                                    { J(2, 0), J(2, 1), J(2, 2) } };
  const QuadrotorSimulator::Quadrotor::State state = quad.getState();

  float R13 = state.R(0,2);
  float R23 = state.R(1,2);
  float R33 = state.R(2,2);
  
  float Om1 = state.omega(0);
  float Om2 = state.omega(1);
  float Om3 = state.omega(2);

  float force = 0;
  // if (Psi < 1.0f) // Position control stability guaranteed only when Psi < 1
    force = cmd.force[0] * R13 + cmd.force[1] * R23 + cmd.force[2] * R33;

  float eOm1 = cmd.omg[0] - Om1;
  float eOm2 = cmd.omg[1] - Om2;
  float eOm3 = cmd.omg[2] - Om3;

  // 叉乘 omg X I.omg
  float in1 = Om2 * (I[2][0] * Om1 + I[2][1] * Om2 + I[2][2] * Om3) -
              Om3 * (I[1][0] * Om1 + I[1][1] * Om2 + I[1][2] * Om3);
  float in2 = Om3 * (I[0][0] * Om1 + I[0][1] * Om2 + I[0][2] * Om3) -
              Om1 * (I[2][0] * Om1 + I[2][1] * Om2 + I[2][2] * Om3);
  float in3 = Om1 * (I[1][0] * Om1 + I[1][1] * Om2 + I[1][2] * Om3) -
              Om2 * (I[0][0] * Om1 + I[0][1] * Om2 + I[0][2] * Om3);

  Eigen::Vector3d torque = BodyRatePIDCtrl(Eigen::Vector3d(cmd.omg[0], cmd.omg[1], cmd.omg[2]),
                                           Eigen::Vector3d(Om1, Om2, Om3), 
                                           Eigen::Vector3d(cmd.kOmP[0], cmd.kOmP[1], cmd.kOmP[2]), 
                                           Eigen::Vector3d(cmd.kOmI[0], cmd.kOmI[1], cmd.kOmI[2]), 
                                           Eigen::Vector3d(cmd.kOmD[0], cmd.kOmD[1], cmd.kOmD[2]), 
                                           Eigen::Vector3d(0, 0, 0),
                                           dt);

  double M1 = torque(0) + in1;
  double M2 = torque(1) + in2;
  double M3 = torque(2) + in3;

  // 十字形 M3 YAW, M1, M2 Roll， Pitch 
  float w_sq[4];
  w_sq[0] = force / (4 * kf) - M2 / (2 * d * kf) + M3 / (4 * km);
  w_sq[1] = force / (4 * kf) + M2 / (2 * d * kf) + M3 / (4 * km);
  w_sq[2] = force / (4 * kf) + M1 / (2 * d * kf) - M3 / (4 * km);
  w_sq[3] = force / (4 * kf) - M1 / (2 * d * kf) - M3 / (4 * km);

  Control control;
  control.force = force;
  control.M1 = M1;
  control.M2 = M2;
  control.M3 = M3;

  control.eR1 = 0;
  control.eR2 = 0;
  control.eR3 = 0;

  control.eOmg1 = eOm1;
  control.eOmg2 = eOm2;
  control.eOmg3 = eOm3;

  control.in1 = in1;
  control.in2 = in2;
  control.in3 = in3;

  for (int i=0; i<4; i++)
  {
    if (w_sq[i]<0)
    {
      control.des_rpm[i] = -sqrtf(-w_sq[i]);
    } 
    else 
    {
      control.des_rpm[i] = sqrtf(w_sq[i]);
    }
    control.real_rpm[i] = state.motor_rpm(i);
  }

  for (int i = 0; i < 4; i++)
  {
    if (w_sq[i] < 0)
      w_sq[i] = 0;

    control.rpm[i] = sqrtf(w_sq[i]);
  }
  return control;
}

static void
cmd_callback(const quadrotor_msgs::SO3Command::ConstPtr& cmd)
{
  command.force[0]         = cmd->force.x;
  command.force[1]         = cmd->force.y;
  command.force[2]         = cmd->force.z;
  command.qx               = cmd->orientation.x;
  command.qy               = cmd->orientation.y;
  command.qz               = cmd->orientation.z;
  command.qw               = cmd->orientation.w;
  command.kOmP[0]          = cmd->kOmP[0];
  command.kOmP[1]          = cmd->kOmP[1];
  command.kOmP[2]          = cmd->kOmP[2];
  command.kOmI[0]          = cmd->kOmI[0];
  command.kOmI[1]          = cmd->kOmI[1];
  command.kOmI[2]          = cmd->kOmI[2];
  command.kOmD[0]          = cmd->kOmD[0];
  command.kOmD[1]          = cmd->kOmD[1];
  command.kOmD[2]          = cmd->kOmD[2];
  // command.kOm[2]           = 0.0;
  command.omg[0]           = cmd->omega.x;
  command.omg[1]           = cmd->omega.y;
  command.omg[2]           = cmd->omega.z;
  command.corrections[0]   = cmd->aux.kf_correction;
  command.corrections[1]   = cmd->aux.angle_corrections[0];
  command.corrections[2]   = cmd->aux.angle_corrections[1];
  command.current_yaw      = cmd->aux.current_yaw;
  command.use_external_yaw = cmd->aux.use_external_yaw;
}

static void
force_disturbance_callback(const geometry_msgs::Vector3::ConstPtr& f)
{
  disturbance.f(0) = f->x;
  disturbance.f(1) = f->y;
  disturbance.f(2) = f->z;
}

static void
moment_disturbance_callback(const geometry_msgs::Vector3::ConstPtr& m)
{
  disturbance.m(0) = m->x;
  disturbance.m(1) = m->y;
  disturbance.m(2) = m->z;
}

int
main(int argc, char** argv)
{
  ros::init(argc, argv, "quadrotor_simulator_so3");

  ros::NodeHandle n("~");

  ros::Publisher  odom_pub = n.advertise<nav_msgs::Odometry>("odom", 100);
  ros::Publisher  imu_pub  = n.advertise<sensor_msgs::Imu>("imu", 10);
  ros::Publisher  des_motor_pub = n.advertise<quadrotor_msgs::Motor>("des_rpm", 10);
  ros::Publisher  real_motor_pub = n.advertise<quadrotor_msgs::Motor>("real_rpm", 10);
  ros::Publisher  quad_status_pub = n.advertise<quadrotor_msgs::QuadStatus>("sim_quad_status", 10);
  ros::Subscriber cmd_sub =
    n.subscribe("cmd", 100, &cmd_callback, ros::TransportHints().tcpNoDelay());
  ros::Subscriber f_sub =
    n.subscribe("force_disturbance", 100, &force_disturbance_callback,
                ros::TransportHints().tcpNoDelay());
  ros::Subscriber m_sub =
    n.subscribe("moment_disturbance", 100, &moment_disturbance_callback,
                ros::TransportHints().tcpNoDelay());

  QuadrotorSimulator::Quadrotor quad;
  double                        _init_x, _init_y, _init_z;
  n.param("simulator/init_state_x", _init_x, 0.0);
  n.param("simulator/init_state_y", _init_y, 0.0);
  n.param("simulator/init_state_z", _init_z, 1.0);

  Eigen::Vector3d position = Eigen::Vector3d(_init_x, _init_y, _init_z);
  quad.setStatePos(position);

  double simulation_rate;
  n.param("rate/simulation", simulation_rate, 1000.0);
  ROS_ASSERT(simulation_rate > 0);

  double odom_rate;
  n.param("rate/odom", odom_rate, 100.0);
  const ros::Duration odom_pub_duration(1 / odom_rate);

  std::string quad_name;
  n.param("quadrotor_name", quad_name, std::string("quadrotor"));

  QuadrotorSimulator::Quadrotor::State state = quad.getState();

  ros::Rate    r(simulation_rate);
  const double dt = 1 / simulation_rate;

  Control control;

  nav_msgs::Odometry odom_msg;
  odom_msg.header.frame_id = "/simulator";
  odom_msg.child_frame_id  = "/" + quad_name;

  sensor_msgs::Imu imu;
  imu.header.frame_id = "/simulator";

  quadrotor_msgs::Motor motors;
  motors.header.frame_id = "/simulator";

  quadrotor_msgs::QuadStatus status;
  status.header.frame_id = "/simulator";

  ros::Time next_odom_pub_time = ros::Time::now();
  ros::Time last_sim_time = next_odom_pub_time;
  while (n.ok())
  {
    ros::spinOnce();

    auto last = control;
    control   = getControl(quad, command, (ros::Time::now() - last_sim_time).toSec());
    for (int i = 0; i < 4; ++i)
    {
      //! @bug might have nan when the input is legal
      if (std::isnan(control.rpm[i]))
        control.rpm[i] = last.rpm[i];
    }
    quad.setInput(control.rpm[0], control.rpm[1], control.rpm[2],
                  control.rpm[3]);
    quad.setExternalForce(disturbance.f);
    quad.setExternalMoment(disturbance.m);
    quad.step(dt);

    motors.max_rpm = quad.getMaxRPM();
    motors.min_rpm = quad.getMinRPM();

    ros::Time tnow = ros::Time::now();
    last_sim_time = tnow;

    if (tnow >= next_odom_pub_time)
    {
      next_odom_pub_time += odom_pub_duration;
      odom_msg.header.stamp = tnow;
      status.header.stamp = tnow;
      state                 = quad.getState();
      stateToOdomMsg(state, odom_msg);
      quadToImuMsg(quad, imu);
      odom_pub.publish(odom_msg);
      imu_pub.publish(imu);

      status.motor_rpm_1 = state.motor_rpm(0);
      status.motor_rpm_2 = state.motor_rpm(1);
      status.motor_rpm_3 = state.motor_rpm(2);
      status.motor_rpm_4 = state.motor_rpm(3);

      status.omega_x = state.omega(0);
      status.omega_y = state.omega(1);
      status.omega_z = state.omega(2);

      status.vel_x = state.v(0);
      status.vel_y = state.v(1);
      status.vel_z = state.v(2);

      status.pos_x = state.x(0);
      status.pos_y = state.x(1);
      status.pos_z = state.x(2);

      status.force = control.force;
      status.M1 = control.M1;
      status.M2 = control.M2;
      status.M3 = control.M3;

      status.eOmg1 = control.eOmg1;
      status.eOmg2 = control.eOmg2;
      status.eOmg3 = control.eOmg3;

      status.eR1 = control.eR1;
      status.eR2 = control.eR2;
      status.eR3 = control.eR3;

      status.In1 = control.in1;
      status.In2 = control.in2;
      status.In3 = control.in3;
      
      quad_status_pub.publish(status);

      motors.header.stamp = tnow;
      motors.rpm_1 = control.des_rpm[0];
      motors.rpm_2 = control.des_rpm[1];
      motors.rpm_3 = control.des_rpm[2];
      motors.rpm_4 = control.des_rpm[3];
      des_motor_pub.publish(motors);
      
      motors.rpm_1 = control.real_rpm[0];
      motors.rpm_2 = control.real_rpm[1];
      motors.rpm_3 = control.real_rpm[2];
      motors.rpm_4 = control.real_rpm[3];
      real_motor_pub.publish(motors);
    }

    r.sleep();
  }

  return 0;
}

void
stateToOdomMsg(const QuadrotorSimulator::Quadrotor::State& state,
               nav_msgs::Odometry&                         odom)
{
  odom.pose.pose.position.x = state.x(0);
  odom.pose.pose.position.y = state.x(1);
  odom.pose.pose.position.z = state.x(2);

  Eigen::Quaterniond q(state.R);
  odom.pose.pose.orientation.x = q.x();
  odom.pose.pose.orientation.y = q.y();
  odom.pose.pose.orientation.z = q.z();
  odom.pose.pose.orientation.w = q.w();

  odom.twist.twist.linear.x = state.v(0);
  odom.twist.twist.linear.y = state.v(1);
  odom.twist.twist.linear.z = state.v(2);

  odom.twist.twist.angular.x = state.omega(0);
  odom.twist.twist.angular.y = state.omega(1);
  odom.twist.twist.angular.z = state.omega(2);
}

void
quadToImuMsg(const QuadrotorSimulator::Quadrotor& quad, sensor_msgs::Imu& imu)

{
  QuadrotorSimulator::Quadrotor::State state = quad.getState();
  Eigen::Quaterniond                   q(state.R);
  imu.orientation.x = q.x();
  imu.orientation.y = q.y();
  imu.orientation.z = q.z();
  imu.orientation.w = q.w();

  imu.angular_velocity.x = state.omega(0);
  imu.angular_velocity.y = state.omega(1);
  imu.angular_velocity.z = state.omega(2);

  imu.linear_acceleration.x = quad.getAcc()[0];
  imu.linear_acceleration.y = quad.getAcc()[1];
  imu.linear_acceleration.z = quad.getAcc()[2];
}
