/**
 * @file Cal3DS2.h
 * @brief Calibration of a camera with radial distortion
 */

#pragma once

#include "../mat/MatCal.h"
#include "../mat/Matrix.h"

namespace minisam
{

/**
 * @brief Calibration of a camera with radial distortion
 * @addtogroup geometry
 * \nosubgrouping
 *
 * Uses same distortionmodel as OpenCV, with
 * http://docs.opencv.org/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html
 * but using only k1,k2,p1, and p2 coefficients.
 * K = [ fx s u0 ; 0 fy v0 ; 0 0 1 ]
 * rr = Pn.x^2 + Pn.y^2
 * \hat{Pn} = (1 + k1*rr + k2*rr^2 ) Pn + [ 2*p1 Pn.x Pn.y + p2 (rr + 2 Pn.x^2) ;
 *                      p1 (rr + 2 Pn.y^2) + 2*p2 Pn.x Pn.y  ]
 * pi = K*Pn
 */
class Cal3DS2_Base:public minivector
{

public:

    /// @name Standard Constructors
    /// @{

    /// Default Constructor with only unit focal length
    Cal3DS2_Base() :minivector(9)
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
    }

    Cal3DS2_Base(double fx, double fy, double s, double u0, double v0,
                 double k1, double k2, double p1 = 0.0, double p2 = 0.0) :minivector(9)
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
    }

    virtual ~Cal3DS2_Base() {}

    /// @}
    /// @name Advanced Constructors
    /// @{

    Cal3DS2_Base(const minivector &v) ;

    /// @}

    /// @name Standard Interface
    /// @{

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

    /// return calibration matrix -- not really applicable
    minimatrix K() const;

    /// return distortion parameter vector
    minivector k() const
    {
        return minivector_subvector(*this,5,4);
    }
    /// Return all parameters as a vector
    minivector vector() const;

    virtual minimatrix* Retract(const minimatrix* mpose)
    {
        minivector d(mpose);
        minivector result(9);
        minivector_add(&result,d,vector());
        return new Cal3DS2_Base(result);
    }
    /// Given a different calibration, calculate update to obtain it
    virtual minimatrix LocalCoordinates(const minimatrix* mpose,minimatrix* H1=NULL,minimatrix* H2=NULL) const
    {
        Cal3DS2_Base bb(mpose);
        minivector fb=bb.vector();
        minivector result(9);
        minivector_sub(&result,fb,vector());
        return result;
    }

    /**
     * convert intrinsic coordinates xy to (distorted) image coordinates uv
     * @param p point in intrinsic coordinates
     * @param Dcal optional 2*9 Jacobian wrpt Cal3DS2 parameters
     * @param Dp optional 2*2 Jacobian wrpt intrinsic coordinates
     * @return point in (distorted) image coordinates
     */
    minivector uncalibrate(const minivector& p,
                           minimatrix* Dcal = NULL,//2*9
                           minimatrix* Dp = NULL//2*2
                          ) const ;

    /// Convert (distorted) image coordinates uv to intrinsic coordinates xy
    minivector calibrate(const minivector& p, const double tol=1e-5) const;

    /// Derivative of uncalibrate wrpt intrinsic coordinates
    minimatrix D2d_intrinsic(const minivector& p) const ;

    /// Derivative of uncalibrate wrpt the calibration parameters
    minimatrix D2d_calibration(const minivector& p) const ;


    /// @}
    /// @name Clone
    /// @{

    /// @return a deep copy of this object
    virtual Cal3DS2_Base* clone() const
    {
        return new Cal3DS2_Base(*this);
    }

    /// @}

};
minimatrix D2dcalibration(double x, double y, double xx,
                          double yy, double xy, double rr, double r4, double pnx, double pny,
                          const minimatrix& DK);
minimatrix D2dintrinsic(double x, double y, double rr,
                        double g, double k1, double k2, double p1, double p2,
                        const minimatrix& DK);
}

