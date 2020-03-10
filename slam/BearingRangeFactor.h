/**
 * @file BearingRangeFactor.h
 * @brief a single factor contains both the bearing and the range to prevent handle to pair bearing and range factors
 */

#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../geometry/Pose2.h"
#include "../geometry/Pose3.h"

namespace minisam {

  /**
   * Binary factor for a 2d bearing measurement
   * @addtogroup SLAM
   */
  class BearingRangeFactor2d: public NoiseModelFactor2
  {

  public:

    // the measurement
    Rot2* measuredBearing_;
    double measuredRange_;

  public:

    BearingRangeFactor2d() {} /* Default constructor */
    BearingRangeFactor2d(int poseKey, int pointKey,  Rot2* measuredBearing, const double measuredRange,
        GaussianNoiseModel* model) :
          NoiseModelFactor2(model, poseKey, pointKey), measuredBearing_(measuredBearing), measuredRange_(measuredRange)
    {
    }

    virtual ~BearingRangeFactor2d()
    {
     if(measuredBearing_!=NULL)
        {
        delete measuredBearing_;
        measuredBearing_=NULL;
        }

    }

    /// @return a deep copy of this factor
    virtual NoiseModelFactor* clone() const
     {
       BearingRangeFactor2d* newfactor=new BearingRangeFactor2d(key1(),key2(),measuredBearing_,measuredRange_,noiseModel_);
        return newfactor;
    }


    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x,
                                       std::vector<minimatrix> &H) const
    {
        std::map<int,minimatrix*>::const_iterator itb1=x.find(key1());
        std::map<int,minimatrix*>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second,
                             *(H.begin()),*(H.begin()+1));
    }

    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x) const
    {
        std::map<int,minimatrix*>::const_iterator itb1=x.find(key1());
        std::map<int,minimatrix*>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second);
    }



    /** h(x)-z -> between(z,h(x)) for Rot manifold */
    virtual minivector evaluateError(const minimatrix* pose, const minimatrix* point) const
     {
      Pose2 ppose(pose);
      minivector ppoint(point);

      Rot2 y1 = ppose.bearing(ppoint);

      minimatrix measured_bw_y1=measuredBearing_->between(&y1);
      Rot2 rmeasured_bw_y1(&measured_bw_y1);
      double e1=Rot2::Logmap(rmeasured_bw_y1);

      double y2 =ppose.range(ppoint);
      double e2 =y2 - measuredRange_;

      minivector result(2);
      result.data[0]=e1;result.data[1]=e2;

      return result;
    }
      /** h(x)-z -> between(z,h(x)) for Rot manifold */
    virtual minivector  evaluateError(const minimatrix* pose,const minimatrix* point,
                                      minimatrix& H1,minimatrix& H2) const {
      minimatrix H11, H21, H12, H22;

      Pose2 ppose(pose);
      minivector ppoint(point);

      Rot2 y1 = ppose.bearing(ppoint, &H11, &H12);//1*3,1*2

      minimatrix measured_bw_y1=measuredBearing_->between(&y1);
      Rot2 rmeasured_bw_y1(&measured_bw_y1);

      double e1 = Rot2::Logmap(rmeasured_bw_y1);

      double y2 = ppose.range(ppoint, &H21, &H22);//1*3,1*2
      double e2 =y2 - measuredRange_;


      minimatrix_resize(&H1,2,3);
      H1.data[0]=H11.data[0];H1.data[1]=H11.data[1];H1.data[2]=H11.data[2];
      H1.data[3]=H21.data[0];H1.data[4]=H21.data[1];H1.data[5]=H21.data[2];
      minimatrix_resize(&H2,2,2);
      H2.data[0]=H12.data[0];H2.data[1]=H12.data[1];
      H2.data[2]=H22.data[0];H2.data[3]=H22.data[1];

      minivector result(2);
      result.data[0]=e1;result.data[1]=e2;

      return result;
    }

    /** return the measured */
    const std::pair<Rot2*, double> measured() const {
      return std::make_pair(measuredBearing_, measuredRange_);
    }

  }; // BearingRangeFactor2d
}
