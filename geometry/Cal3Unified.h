
/**
 * @file Cal3Unified.h
 * @brief Unified Calibration Model, see Mei07icra for details
 */

/**
 * @addtogroup geometry
 */

#pragma once

#include "../geometry/Cal3DS2_Base.h"

namespace minisam
{

/**
 * @brief Calibration of a omni-directional camera with mirror + lens radial distortion
 * @addtogroup geometry
 * \nosubgrouping
 *
 * Similar to Cal3DS2, does distortion but has additional mirror parameter xi
 * K = [ fx s u0 ; 0 fy v0 ; 0 0 1 ]
 * Pn = [ P.x / (1 + xi * \sqrt{P.x^2 + P.y^2 + 1}), P.y / (1 + xi * \sqrt{P.x^2 + P.y^2 + 1})]
 * rr = Pn.x^2 + Pn.y^2
 * \hat{pn} = (1 + k1*rr + k2*rr^2 ) pn + [ 2*k3 pn.x pn.y + k4 (rr + 2 Pn.x^2) ;
 *                      k3 (rr + 2 Pn.y^2) + 2*k4 pn.x pn.y  ]
 * pi = K*pn
 */
class  Cal3Unified : public minivector
{

public:

    /// @name Standard Constructors
    /// @{

    /// Default Constructor with only unit focal length
    Cal3Unified() :minivector(10)
    {
        data[0]=1.0;
        data[1]=1.0;
        data[2]=0.0;
        data[3]=0.0;
        data[4]=0.0;
        data[5]=0.0;
        data[6]=0.0;
        data[7]=0.0;
        data[8]=0.0;
        data[9]=0.0;

    }

    Cal3Unified(double fx, double fy, double s, double u0, double v0,
                double k1, double k2, double p1 = 0.0, double p2 = 0.0, double xi = 0.0) :
        minivector(10)
    {
        data[0]=fx;
        data[1]=fy;
        data[2]=s;
        data[3]=u0;
        data[4]=v0;
        data[5]=k1;
        data[6]=k2;
        data[7]=p1;
        data[8]=p2;
        data[9]=xi;
    }

    virtual ~Cal3Unified() {}

    /// @}
    /// @name Advanced Constructors
    /// @{

    Cal3Unified(const minivector &v) ;

    /// @}


    /// @name Standard Interface
    /// @{

    /// mirror parameter
    /// focal length x
    inline double fx() const
    {
        return data[0];
    }

    /// focal length x
    inline double fy() const
    {
        return data[1];
    }

    /// skew
    inline double skew() const
    {
        return data[2];
    }

    /// image center in x
    inline double px() const
    {
        return data[3];
    }

    /// image center in y
    inline double py() const
    {
        return data[4];
    }

    /// First distortion coefficient
    inline double k1() const
    {
        return data[5];
    }

    /// Second distortion coefficient
    inline double k2() const
    {
        return data[6];
    }

    /// First tangential distortion coefficient
    inline double p1() const
    {
        return data[7];
    }

    /// Second tangential distortion coefficient
    inline double p2() const
    {
        return data[8];
    }

    inline double xi() const
    {
        return data[9];
    }

    /**
     * convert intrinsic coordinates xy to image coordinates uv
     * @param p point in intrinsic coordinates
     * @param Dcal optional 2*10 Jacobian wrpt Cal3Unified parameters
     * @param Dp optional 2*2 Jacobian wrpt intrinsic coordinates
     * @return point in image coordinates
     */
    minivector uncalibrate(const minivector& p,
                           minimatrix* Dcal=NULL,
                           minimatrix* Dp=NULL) const;
    minivector baseuncalibrate(const minivector& p,
                               minimatrix* Dcal=NULL,
                               minimatrix* Dp=NULL) const;
    /// Conver a pixel coordinate to ideal coordinate
    minivector calibrate(const minivector& p, const double tol=1e-5) const;

    /// Convert a 3D point to normalized unit plane
    minivector spaceToNPlane(const minivector& p) const;

    /// Convert a normalized unit plane point to 3D space
    minivector nPlaneToSpace(const minivector& p) const;

    /// @}
    /// @name Manifold
    /// @{

    /// Given delta vector, update calibration
    virtual minimatrix* Retract(const minimatrix* mpose);


    /// Given a different calibration, calculate update to obtain it
    virtual minimatrix LocalCoordinates(const minimatrix* mpose) const;

    /// Return dimensions of calibration manifold object
    size_t dim() const
    {
        return 10 ;
    }


    /// Return all parameters as a vector
    minivector vector() const ;

    minimatrix D2d_calibration(const minivector& p) const;

    /// @}
};


};

