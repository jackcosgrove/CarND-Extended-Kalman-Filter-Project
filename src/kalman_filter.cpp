#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::InitPredict(MatrixXd &F_in, MatrixXd &Q_in) {
  F_ = F_in;
  Q_ = Q_in;
}

void KalmanFilter::InitUpdate(MatrixXd &H_in, MatrixXd &R_in, MatrixXd &I_in) {
  H_ = H_in;
  R_ = R_in;
  I_ = I_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  
  MatrixXd Ht = H_.transpose();
  
  MatrixXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  
  VectorXd h(3);
  
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float rho = sqrt(px*px+py*py);
  
  float phi = atan2(py, px);
  if (fabs(px) < 0.001)
    phi = atan(py/0.001);
  
  float rho_dot = 0.0;
  if (fabs(rho) > 0.001)
    rho_dot = (px*vx + py*vy)/rho;
  
  h << rho, phi, rho_dot;
  VectorXd y = z - h;
  
  while (y(1) > M_PI) {
    y(1) -= 2 * M_PI;
  }
  while (y(1) < -M_PI) {
    y(1) += 2 * M_PI;
  }
  
  MatrixXd Ht = H_.transpose();
  
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}
