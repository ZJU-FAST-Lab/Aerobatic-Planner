#include <iostream>
#include <so3_control/SO3Control.h>

#include <ros/ros.h>

SO3Control::SO3Control()
    : mass_(0.98), g_(9.81)
{
  acc_.setZero();
  omega_.setZero();
  old_omega_.setZero();
}

void SO3Control::setMass(const double mass)
{
  mass_ = mass;
}

void SO3Control::setGravity(const double g)
{
  g_ = g;
}

void SO3Control::setPosition(const Eigen::Vector3d &position)
{
  pos_ = position;
}

void SO3Control::setVelocity(const Eigen::Vector3d &velocity)
{
  vel_ = velocity;
}

void SO3Control::setQuat(const Eigen::Quaterniond &q)
{
  quat_ = q;
}

void SO3Control::calculateControl(const Eigen::Vector3d &des_pos,
                                  const Eigen::Vector3d &des_vel,
                                  const Eigen::Vector3d &des_acc,
                                  const Eigen::Vector3d &des_jer,
                                  const Eigen::Vector3d &des_dir,
                                  const Eigen::Vector3d &des_dir_dot,
                                  const double des_yaw, const double des_yaw_dot,
                                  const Eigen::Vector3d &kx,
                                  const Eigen::Vector3d &kv,
                                  const Eigen::Vector3d &kR)
{
  Eigen::Vector3d acc_error;

  // pos_err nan->0.0, else bound(pd-pr, -1, 1)
  // vel_err (vd+kp*pos_err - vr) bound in [-1, 1]
  // acc_err kv*vel_err

  // x acceleration
  double x_pos_error = std::isnan(des_pos(0)) ? 0.0 : std::max(std::min(des_pos(0) - pos_(0), 1.0), -1.0);
  double x_vel_error = std::max(std::min((des_vel(0) + kx(0) * x_pos_error) - vel_(0), 1.0), -1.0);
  acc_error(0) = kv(0) * x_vel_error;

  // y acceleration
  double y_pos_error = std::isnan(des_pos(1)) ? 0.0 : std::max(std::min(des_pos(1) - pos_(1), 1.0), -1.0);
  double y_vel_error = std::max(std::min((des_vel(1) + kx(1) * y_pos_error) - vel_(1), 1.0), -1.0);
  acc_error(1) = kv(1) * y_vel_error;

  // z acceleration
  double z_pos_error = std::isnan(des_pos(2)) ? 0.0 : std::max(std::min(des_pos(2) - pos_(2), 1.0), -1.0);
  double z_vel_error = std::max(std::min((des_vel(2) + kx(2) * z_pos_error) - vel_(2), 1.0), -1.0);
  acc_error(2) = kv(2) * z_vel_error;

  // std::cout << "kv = " << kv.transpose() << std::endl;
  force_ = mass_ * g_ * Eigen::Vector3d(0, 0, 1);
  force_ += mass_ * des_acc + acc_error;

  Eigen::Vector4d quat;
  double thr;
  
  if (forward(des_vel, des_acc + acc_error / mass_, des_jer, des_dir, des_dir_dot, yaw_, yaw_dot_, thr, quat, omega_))
  // if (forward(des_vel, des_acc, des_jer, des_dir, des_dir_dot, yaw_, yaw_dot_, thr, quat, omega_))
  {
    orientation_ = Eigen::Quaterniond(quat(0), quat(1), quat(2), quat(3));

    const Eigen::Quaterniond qerr = quat_.inverse() * orientation_;
    Eigen::Matrix3d trans_R(qerr);
    const Eigen::Vector3d temp_bodyrates = omega_;
    omega_ = trans_R * temp_bodyrates;

    Eigen::AngleAxisd rotation_vector(qerr); // debug
    Eigen::Vector3d axis = rotation_vector.axis();

    omega_.x() += kR(0) * rotation_vector.angle() * axis(0);
    omega_.y() += kR(1) * rotation_vector.angle() * axis(1);
    omega_.z() += kR(2) * rotation_vector.angle() * axis(2);

    if (omega_.z() > M_PI)
      omega_.z() = M_PI;
    if (omega_.z() < -M_PI)
      omega_.z() = -M_PI;
    old_omega_.x() = omega_.x();
    old_omega_.y() = omega_.y();
    old_omega_.z() = omega_.z();
  }
  else
  {
    omega_.x() = old_omega_.x();
    omega_.y() = old_omega_.y();
    omega_.z() = old_omega_.z();
  }
}

bool SO3Control::forward(const Eigen::Vector3d &vel,
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
  double mass = 0.98, grav = 9.81, dh = 0, dv = 0, cp = 0, veps = 0.0001;
  double w0, w1, w2, dw0, dw1, dw2;
  double v0, v1, v2, a0, a1, a2, v_dot_a;
  double z0, z1, z2, dz0, dz1, dz2;
  double x0, x1, x2, dx0, dx1, dx2;
  double qw, qx, qy, dqw, dqx, dqy;
  double cp_term, w_term, dh_over_m;
  double zu_sqr_norm, zu_norm, zu0, zu1, zu2;
  double zu_sqr0, zu_sqr1, zu_sqr2, zu01, zu12, zu02;
  double ng00, ng01, ng02, ng11, ng12, ng22, ng_den;
  double dw_term, dz_term0, dz_term1, dz_term2, f_term0, f_term1, f_term2;
  double tilt_den, tilt0, tilt1, tilt2, c_half_psi, s_half_psi;
  double dtilt_den, dtilt0, dtilt1, dtilt2;
  double c_psi, s_psi, omg_den, omg_term;
  double ddir_xb_norm, dir_dot_zb, dir_xb0, dir_xb1, dir_xb2;
  double dir_xb_norm, ddir_xb0, ddir_xb1, ddir_xb2, ddir_dot_zb;
  double xb_dot_dir_xb, dxb_dot_dir_xb;
  double xb_cross_dir_xb0, xb_cross_dir_xb1, xb_cross_dir_xb2;
  double temp, dtemp, sign;
  v0 = vel(0);
  v1 = vel(1);
  v2 = vel(2);
  a0 = acc(0);
  a1 = acc(1);
  a2 = acc(2);
  // w_term 函数 delta(||v||)， w表示 delta() * V
  cp_term = sqrt(v0 * v0 + v1 * v1 + v2 * v2 + veps);
  w_term = 1.0 + cp * cp_term;
  w0 = w_term * v0;
  w1 = w_term * v1;
  w2 = w_term * v2;
  // zu = a + db / m * delta(|v|) + ge3
  dh_over_m = dh / mass;
  zu0 = a0 + dh_over_m * w0;
  zu1 = a1 + dh_over_m * w1;
  zu2 = a2 + dh_over_m * w2 + grav;
  // 除以模长得到 z
  zu_sqr0 = zu0 * zu0;
  zu_sqr1 = zu1 * zu1;
  zu_sqr2 = zu2 * zu2;
  zu01 = zu0 * zu1;
  zu12 = zu1 * zu2;
  zu02 = zu0 * zu2;
  zu_sqr_norm = zu_sqr0 + zu_sqr1 + zu_sqr2;
  zu_norm = sqrt(zu_sqr_norm);
  // real zb = [z0, z1, z2]^T
  z0 = zu0 / zu_norm;
  z1 = zu1 / zu_norm;
  z2 = zu2 / zu_norm;
  // 归一化函数的导数
  ng_den = zu_sqr_norm * zu_norm;
  ng00 = (zu_sqr1 + zu_sqr2) / ng_den;
  ng01 = -zu01 / ng_den;
  ng02 = -zu02 / ng_den;
  ng11 = (zu_sqr0 + zu_sqr2) / ng_den;
  ng12 = -zu12 / ng_den;
  ng22 = (zu_sqr0 + zu_sqr1) / ng_den;
  // 算 delta( ||v|| ) * v 的导数
  v_dot_a = v0 * a0 + v1 * a1 + v2 * a2;
  dw_term = cp * v_dot_a / cp_term;
  dw0 = w_term * a0 + dw_term * v0;
  dw1 = w_term * a1 + dw_term * v1;
  dw2 = w_term * a2 + dw_term * v2;
  // 算 zu 的导数
  dz_term0 = jer(0) + dh_over_m * dw0;
  dz_term1 = jer(1) + dh_over_m * dw1;
  dz_term2 = jer(2) + dh_over_m * dw2;
  // 最终 zb 的导数 链式法则
  dz0 = ng00 * dz_term0 + ng01 * dz_term1 + ng02 * dz_term2;
  dz1 = ng01 * dz_term0 + ng11 * dz_term1 + ng12 * dz_term2;
  dz2 = ng02 * dz_term0 + ng12 * dz_term1 + ng22 * dz_term2;

  // 计算推力 f = zb * (ma + dv delta()v + mge3)
  f_term0 = mass * a0 + dv * w0;
  f_term1 = mass * a1 + dv * w1;
  f_term2 = mass * (a2 + grav) + dv * w2;
  thr = z0 * f_term0 + z1 * f_term1 + z2 * f_term2;

  // 计算姿态q
  tilt_den = sqrt(2.0 * (1.0 + z2));
  tilt0 = 0.5 * tilt_den;
  tilt1 = -z1 / tilt_den;
  tilt2 = z0 / tilt_den;
  // if (z2 + 1.0 < 1.0e-5)
  //   dtilt_den = 0.0;
  // else
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

  // dir_xb = dir - dir.dot(zb) * zb;
  dir_dot_zb = dir(0) * z0 + dir(1) * z1 + dir(2) * z2;
  dir_xb0 = dir(0) - dir_dot_zb * z0;
  dir_xb1 = dir(1) - dir_dot_zb * z1;
  dir_xb2 = dir(2) - dir_dot_zb * z2;
  temp = dir_xb0 * dir_xb0 + dir_xb1 * dir_xb1 + dir_xb2 * dir_xb2;
  dir_xb_norm = sqrt(temp);
  ddir_dot_zb = z0 * ddir(0) + dir(0) * dz0 + z1 * ddir(1) + dir(1) * dz1 + z2 * ddir(2) + dir(2) * dz2;
  ddir_xb0 = ddir(0) - z0 * ddir_dot_zb - dir_dot_zb * dz0;
  ddir_xb1 = ddir(1) - z1 * ddir_dot_zb - dir_dot_zb * dz1;
  ddir_xb2 = ddir(2) - z2 * ddir_dot_zb - dir_dot_zb * dz2;
  dtemp = 2 * dir_xb0 * ddir_xb0 + 2 * dir_xb1 * ddir_xb1 + 2 * dir_xb2 * ddir_xb2;
  // if (temp < 1.0e-5)
  //   ddir_xb_norm = 0.0;
  // else
  ddir_xb_norm = dtemp / (2.0 * dir_xb_norm);

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
  if (1.0 - xb_dot_dir_xb < 1.0e-5 || xb_dot_dir_xb + 1.0 < 1.0e-5)
  {
    dpsi = 0.0;
    return false;
  }
  else
    dpsi = -sign * (dxb_dot_dir_xb / sqrt(1.0 - xb_dot_dir_xb * xb_dot_dir_xb));

  c_half_psi = cos(0.5 * psi);
  s_half_psi = sin(0.5 * psi);
  quat(0) = tilt0 * c_half_psi;
  quat(1) = tilt1 * c_half_psi + tilt2 * s_half_psi;
  quat(2) = tilt2 * c_half_psi - tilt1 * s_half_psi;
  quat(3) = tilt0 * s_half_psi;

  // 计算角速度omega
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

const Eigen::Vector3d &
SO3Control::getComputedForce(void)
{
  return force_;
}

const double &
SO3Control::getComputedYaw(void)
{
  return yaw_;
}

const double &
SO3Control::getComputedYawDot(void)
{
  return yaw_dot_;
}

const Eigen::Quaterniond &
SO3Control::getComputedOrientation(void)
{
  return orientation_;
}

const Eigen::Vector3d &
SO3Control::getComputedOmega(void)
{
  return omega_;
}

void SO3Control::setAcc(const Eigen::Vector3d &acc)
{
  acc_ = acc;
}
