#ifndef CAL3_S2STEREO_H
#define CAL3_S2STEREO_H


/**
 * @file   Cal3_S2Stereo.h
 * @brief  The most common 5DOF 3D->2D calibration + Stereo baseline
 * @author Chris Beall
 */


#pragma once

#include "../geometry/Cal3_S2.h"
#include <iosfwd>


namespace minisam
{

class Cal3_S2Stereo:public minivector
{

public:

    /// @name Standard Constructors
    /// @{

    /// default calibration leaves coordinates unchanged
    Cal3_S2Stereo() :minivector(6)
    {
        data[0]=1.0;
        data[1]=1.0;
        data[2]=0.0;
        data[3]=0.0;
        data[4]=0.0;
        data[5]=1.0;
    }

    /// constructor from doubles
    Cal3_S2Stereo(double fx, double fy, double s, double u0, double v0, double b) :minivector(6)
    {
        data[0]=fx;
        data[1]=fy;
        data[2]=s;
        data[3]=u0;
        data[4]=v0;
        data[5]=b;
    }

    /// constructor from vector
    Cal3_S2Stereo(const minivector &d):minivector(d)
    {
    }

    /// easy constructor; field-of-view in degrees, assumes zero skew
    Cal3_S2Stereo(double fov, int w, int h, double b) :minivector(6)
    {
        double a = fov * M_PI / 360.0; // fov/2 in radians
        data[0]=0.5*(double) w /tan(a);//fx_ = (double) w / (2.0 * tan(a)); //
        data[1]=data[0];// fy_ = fx_;
        data[2]=0;//s_(0),
        data[3]=0.5*(double)w;//u0_((double) w / 2.0),
        data[4]=0.5*(double)h;// v0_((double) h / 2.0)s
        data[5]=b;
    }


    /// @}


    /// @name Standard Interface
    /// @{

    /// return calibration, same for left and right
    const Cal3_S2& calibration() const
    {
        minivector result=minivector_subvector(*this,0,5);
        return static_cast<Cal3_S2>(result);
    }

    /// return calibration matrix K, same for left and right
    minimatrix matrix() const
    {
        minivector result=minivector_subvector(*this,0,5);
        return Cal3_S2(&result).matrix();
    }

    /// focal length x
    inline double fx() const
    {
        return minivector_get(this,0);
    }

    /// focal length x
    inline double fy() const
    {
        return minivector_get(this,1);
    }

    /// skew
    inline double skew() const
    {
        return minivector_get(this,2);
    }

    /// image center in x
    inline double px() const
    {
        return minivector_get(this,3);
    }

    /// image center in y
    inline double py() const
    {
        return minivector_get(this,4);
    }

    /// return baseline
    inline double baseline() const
    {
        return minivector_get(this,5);
    }

    /// vectorized form (column-wise)
    minivector vector() const
    {
        return *this;
    }

    /// @}
    /// @name Manifold
    /// @{

    /// return DOF, dimensionality of tangent space
    inline int dim() const
    {
        return 6;
    }


    virtual minimatrix* Retract(const minimatrix* mpose) const
    {
        minivector d(mpose);
        return new Cal3_S2Stereo(fx() + minivector_get(&d,0),
                                 fy() +  minivector_get(&d,1),
                                 skew() +  minivector_get(&d,2),
                                 px() + minivector_get(&d,3),
                                 py() +  minivector_get(&d,4),
                                 baseline() +  minivector_get(&d,5));
    }

    /// Unretraction for the calibration
    virtual minimatrix LocalCoordinates(const minimatrix* T2)
    {
        minivector t2v(T2);
        minivector result(2);
        minivector_sub(&result,t2v,vector());
        return t2v;
    }
    /// @}


};
};
#endif // CAL3_S2STEREO_H
