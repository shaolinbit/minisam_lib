/**
 * @file    GenericStereoFactor.h
 * @brief   A non-linear factor for stereo measurements
 */

#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../geometry/StereoCamera.h"

namespace minisam {

/**
 * A Generic Stereo Factor
 * @addtogroup SLAM
 */
class GenericStereoFactor: public NoiseModelFactor2
{
private:

  // Keep a copy of measurement and calibration for I/O
  StereoPoint2* measured_;                      ///< the measurement
  Cal3_S2Stereo* K_;                ///< shared pointer to calibration
  Pose3* body_P_sensor_;        ///< The pose of the sensor in the body frame


public:


  /**
   * Default constructor
   */
  GenericStereoFactor() : K_(new Cal3_S2Stereo(444, 555, 666, 777, 888, 1.0)),measured_(NULL),body_P_sensor_(NULL)
  {

  }

  /**
   * Constructor
   * @param measured is the Stereo Point measurement (u_l, u_r, v). v will be identical for left & right for rectified stereo pair
   * @param model is the noise model in on the measurement
   * @param poseKey the pose variable key
   * @param landmarkKey the landmark variable key
   * @param K the constant calibration
   * @param body_P_sensor is the transform from body to sensor frame (default identity)
   */
  GenericStereoFactor(StereoPoint2* measured, GaussianNoiseModel*  model,
      int poseKey, int landmarkKey, Cal3_S2Stereo* K,
      Pose3* body_P_sensor=NULL) :
    NoiseModelFactor2(model, poseKey, landmarkKey), measured_(measured), K_(K), body_P_sensor_(body_P_sensor)
     {}


  /** Virtual destructor */
  virtual ~GenericStereoFactor()
  {
  if(measured_!=NULL)
  {
      delete measured_;
      measured_=NULL;
  }
  }

  /// @return a deep copy of this factor
  virtual NoiseModelFactor* clone() const {
    return new GenericStereoFactor(measured_,noiseModel_,key1(),key2(),K_,body_P_sensor_);

    }

    /** h(x)-z */
   virtual minivector
    evaluateError(const minimatrix* pose3, const minimatrix* point3) const
   {
    Pose3 ppose(pose3);
    minivector ppoint(point3);

    try {
      if(body_P_sensor_!=NULL) {
          StereoCamera stereoCam(ppose.compose(*body_P_sensor_), *K_);
          return (stereoCam.project(ppoint) - (*measured_)).vector();

      } else {
        StereoCamera stereoCam(ppose, *K_);
        return (stereoCam.project(ppoint) - (*measured_)).vector();
      }
    } catch(runtime_error& e) {
      if (DEBUGSTATERE)
      std::cout << e.what() << ": Landmark "<< this->key2()<<
          " moved behind camera " << this->key1() << std::endl;
        throw e;
    }
   // return Vector3::Constant(2.0 * K_->fx());
   minivector result(3,2.0*(K_->fx()));
   return result;
  }



  /** h(x)-z */
   virtual minivector
    evaluateError(const minimatrix* pose3, const minimatrix* point3, minimatrix &H1, minimatrix &H2) const
   {
         Pose3 ppose(pose3);
    minivector ppoint(point3);
    try {
      if(body_P_sensor_!=NULL) {
          minimatrix H0;
          minimatrix H1temp;
          StereoCamera stereoCam(ppose.compose(*body_P_sensor_, &H0), *K_);
          StereoPoint2 reprojectionError(stereoCam.project(ppoint, &H1temp, &H2) - (*measured_));
          //*H1 = *H1 * H0;
          miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,H1temp,H0,0.0,&H1);
          return reprojectionError.vector();
      } else {
        StereoCamera stereoCam(ppose, *K_);
        return (stereoCam.project(ppoint, &H1, &H2) - (*measured_)).vector();
      }
    } catch(runtime_error& e) {
      //if (H1) *H1 = Matrix::Zero(3,6);
      minimatrix_resize(&H1,3,6);
      minimatrix_set_zero(&H1);
     // if (H2) *H2 = Z_3x3;
      minimatrix_resize(&H2,3,3);
      minimatrix_set_zero(&H2);
      if (DEBUGSTATERE)
      std::cout << e.what() << ": Landmark "<< this->key2()<<
          " moved behind camera " << this->key1() << std::endl;
        throw e;
    }
   // return Vector3::Constant(2.0 * K_->fx());
    minivector result(3,2.0*(K_->fx()));
   return result;
  }

  /** return the measured */
   StereoPoint2* measured()  {
    return measured_;
  }

  /** return the calibration object */
  inline  Cal3_S2Stereo* calibration(){
    return K_;
  }


};


}; // \ namespace minisam
