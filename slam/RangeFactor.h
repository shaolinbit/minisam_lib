/**
 *  @file  RangeFactor.H
 **/

#pragma once
#include "../nonlinear/NonlinearFactor.h"
#include "../geometry/Pose2.h"

namespace minisam {

  /**
   * Binary factor for a range measurement
   * @addtogroup SLAM
   */
  class RangeFactor2D: public NoiseModelFactor2{
  private:

    double measured_; /** measurement */
    Pose2* body_P_sensor_; ///< The pose of the sensor in the body frame


  public:

    RangeFactor2D() {} /* Default constructor */

    RangeFactor2D(int poseKey, int pointKey, double measured,
        GaussianNoiseModel* model, Pose2* body_P_sensor =NULL) :
          NoiseModelFactor2(model, poseKey, pointKey), measured_(measured), body_P_sensor_(body_P_sensor) {
    }

    virtual ~RangeFactor2D() {}

    /// @return a deep copy of this factor
    virtual NoiseModelFactor* clone() const {
      return new RangeFactor2D(key1(),key2(),measured_,noiseModel_,body_P_sensor_);
       }

    /** h(x)-z */
   // minivector evaluateError(const POSE& pose, const POINT& point,
   //     boost::optional<Matrix&> H1 = boost::none, boost::optional<Matrix&> H2 = boost::none) const
    virtual minivector
    evaluateError(const minimatrix* pose, const minimatrix* point, minimatrix &H1, minimatrix &H2) const
    {
      double hx;
      Pose2 ppose(pose);
      minivector ppoint(point);
      if(body_P_sensor_!=NULL) {
       // if(H1) {
          minimatrix H0;//3*3
          Pose2* pcomp=ppose.compose(*body_P_sensor_, &H0);
          minimatrix H1temp(1,3);
          hx = pcomp->range(ppoint, &H1temp, &H2);//1*3,1*2
          //hx = ppose.compose(*body_P_sensor_, &H0).range(point, H1, H2);
         minimatrix_resize(&H1,1,3);
         miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,H1temp,H0,0.0,&H1);
         // *H1 = *H1 * H0;
         delete pcomp;


      //  } else {
       //   hx = ppose.compose(*body_P_sensor_).range(point, H1, H2);
      //  }
      } else {
        hx = ppose.range(ppoint, &H1, &H2);
      }
      //return (Vector(1) << hx - measured_);
      return minivector(1,hx-measured_);
    }
     virtual minivector
    evaluateError(const minimatrix* pose, const minimatrix* point) const
    {
        /*
        cout<<"RangeFactor key"<<endl;
        cout<<symbolChr(key1())<<endl;
        cout<<symbolIndex(key1())<<endl;
        cout<<symbolChr(key2())<<endl;
        cout<<symbolIndex(key2())<<endl;
        */

      double hx;
      Pose2 ppose(pose);
      minivector ppoint(point);
      if(body_P_sensor_!=NULL) {
        Pose2* pcomp=ppose.compose(*body_P_sensor_);
        hx=pcomp->range(ppoint);
        delete pcomp;
        //  hx = ppose.compose(*body_P_sensor_).range(point, H1, H2);
      } else {
        hx = ppose.range(ppoint);
      }
      //return (Vector(1) << hx - measured_);
      return minivector(1,hx-measured_);
    }

    /** return the measured */
    double measured() const {
      return measured_;
    }

  }; // RangeFactor

} // namespace gtsam
