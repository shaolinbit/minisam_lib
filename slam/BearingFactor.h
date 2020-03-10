/**
 *  @file  BearingFactor.H
 **/

#pragma once

#include "../nonlinear/NonlinearFactor.h"

namespace minisam {

  /**
   * Binary factor for a bearing measurement
   * @addtogroup SLAM
   */

  class BearingFactor2D: public NoiseModelFactor2
  {
   public:
    Rot2* measured_; /** measurement */


    /** default constructor for serialization/testing only */
    BearingFactor2D() {}

    /** primary constructor */
    BearingFactor2D(int poseKey, int pointKey, Rot2* measured,
        GaussianNoiseModel* model) :
          NoiseModelFactor2(model, poseKey, pointKey), measured_(measured) {
    }

    virtual ~BearingFactor2D()
     {
        if(measured_!=NULL)
        {
        delete measured_;
        measured_=NULL;
        }
    }

    /// @return a deep copy of this factor
    virtual NoiseModelFactor* clone() const
    {
      return new BearingFactor2D(key1(),key2(),measured_,noiseModel_);
    }

    /** h(x)-z -> between(z,h(x)) for Rot2 manifold */
     virtual minivector evaluateError(const minimatrix* pose, const minimatrix* point) const
    {
      Pose2 ppose(pose);
      minivector ppoint(point);
      Rot2 hx = ppose.bearing(ppoint);
      minimatrix mbe=measured_->between(&hx);
      return minivector(1,Rot2::Logmap(Rot2(&mbe)));
    }

    virtual minivector evaluateError(const minimatrix* pose, const minimatrix* point, minimatrix &H1, minimatrix &H2) const
    {
      Pose2 ppose(pose);
      minivector ppoint(point);
      Rot2 hx = ppose.bearing(ppoint, &H1, &H2);
      minimatrix mbe=measured_->between(&hx);
      return minivector(1,Rot2::Logmap(Rot2(&mbe)));
    }


    /** return the measured */
    Rot2* measured()
    {
      return measured_;
    }


  }; // BearingFactor

} // namespace minisam
