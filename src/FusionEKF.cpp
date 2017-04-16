#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  P_ = MatrixXd(4, 4);
  F_ = MatrixXd(4, 4);
  I_ = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
        0, 1, 0, 0;
  
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;
  
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  I_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
  
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);
      
      ekf_.x_ << rho * cos(phi), rho * sin(phi), rho_dot * cos(phi), rho_dot * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
    }
      
    previous_timestamp_ = measurement_pack.timestamp_;
    
    ekf_.P_ = P_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float noise_ax = 9.0;
  float noise_ay = 9.0;
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  
  previous_timestamp_ = measurement_pack.timestamp_;
  
  F_(0, 2) = dt;
  F_(1, 3) = dt;
  
  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;
  
  MatrixXd Q(4, 4);
  Q << dt4/4*noise_ax, 0, dt3/2*noise_ax, 0,
        0, dt4/4*noise_ay, 0, dt3/2*noise_ay,
        dt3/2*noise_ax, 0, dt2*noise_ax, 0,
        0, dt3/2*noise_ay, 0, dt2*noise_ay;
  
  float px = ekf_.x_(0);
  float py = ekf_.x_(1);
  if (px == 0 && py == 0) {
    ekf_.x_(0) = ekf_.x_(1) = 0.001;
  }
  
  ekf_.InitPredict(F_, Q);

  ekf_.Predict();
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    
    MatrixXd Hj = tools.CalculateJacobian(ekf_.x_);
    
    ekf_.InitUpdate(Hj, R_radar_, I_);
    
    VectorXd z(3);
    z << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1),
          measurement_pack.raw_measurements_(2);
    
    ekf_.UpdateEKF(z);
    
  } else {
    // Laser updates
    
    ekf_.InitUpdate(H_laser_, R_laser_, I_);
    
    VectorXd z(2);
    z << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1);
    
    ekf_.Update(z);
    
  }
  
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
