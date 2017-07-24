#include "kalman_filter.h"
#include <math.h>


using Eigen::MatrixXd;
using Eigen::VectorXd;

float cyclicAngle(float angle){
    if (angle > PI)
        angle -= 2 * PI;
    else if (angle < -PI)
        angle += 2 * PI;
    return angle;
}

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
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
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    VectorXd z_pred = VectorXd(3);
    cout<<LARGENUM<<endl;
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    float rhodot;
    float rho = sqrt(px * px + py * py);
    float phi = atan2(py,px);
    //cyclic angle
    phi = cyclicAngle(phi);

    float sqrtsq = sqrt(px * px + py * py);
    if (sqrtsq != 0){
    	rhodot = (px * vx + py * vy) / sqrtsq;
    }
    else {
        rhodot = LARGENUM;
    }
    z_pred << rho, phi, rhodot;

    VectorXd y = z - z_pred;

    //cyclic angle
    y(1) = cyclicAngle(y(1));

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
