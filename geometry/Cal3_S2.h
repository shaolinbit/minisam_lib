#ifndef CAL3_S2_H
#define CAL3_S2_H



/**
 * @brief The most common 5DOF 3D->2D calibration
 * @addtogroup geometry
 * \nosubgrouping
 */

#pragma once
#include "../mat/Matrix.h"
#include "../mat/MatCal.h"

namespace minisam
{

class Cal3_S2:public minivector
{
public:

    /// @name Standard Constructors
    /// @{

    /// Create a default calibration that leaves coordinates unchanged
    Cal3_S2() :minivector(5)
    {
        data[0]=1.0;
        data[1]=1.0;
        data[2]=0.0;
        data[3]=0.0;
        data[4]=0.0;
    }

    /// constructor from doubles
    Cal3_S2(double fx, double fy, double s, double u0, double v0) :minivector(5)
    {
        data[0]=fx;
        data[1]=fy;
        data[2]=s;
        data[3]=u0;
        data[4]=v0;
    }

    /// constructor from vector
    Cal3_S2(const minivector &d):minivector(5)
    {
        minivector_memcpy(this,d);
    }

    /**
     * Easy constructor, takes fov in degrees, asssumes zero skew, unit aspect
     * @param fov field of view in degrees
     * @param w image width
     * @param h image height
     */
    Cal3_S2(double fov, int w, int h);

    /// @}
    /// @}
    /// @name Testable
    /// @{

    /// Output stream operator
    friend std::ostream &operator<<(std::ostream &os, const Cal3_S2& cal);

    /// @}
    /// @name Standard Interface
    /// @{

    /// focal length x
    inline double fx() const
    {
        return minivector_get(this,0);
    }

    /// focal length y
    inline double fy() const
    {
        return minivector_get(this,1);
    }

    /// aspect ratio
    inline double aspectRatio() const
    {
        return fx()/fy();
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


    /// vectorized form (column-wise)
    minivector vector() const;

    /// return calibration matrix K
    minimatrix K() const;

    /** @deprecated The following function has been deprecated, use K above */
    minimatrix matrix() const;

    /**
     * convert intrinsic coordinates xy to image coordinates uv, fixed derivaitves
     * @param p point in intrinsic coordinates
     * @param Dcal optional 2*5 Jacobian wrpt Cal3_S2 parameters
     * @param Dp optional 2*2 Jacobian wrpt intrinsic coordinates
     * @return point in image coordinates
     */
    minivector uncalibrate(const minivector& p);
    minivector uncalibrate(const minivector& p, minimatrix* Dcal=NULL,
                           minimatrix* Dp=NULL) const;
    /**
     * convert image coordinates uv to intrinsic coordinates xy
     * @param p point in image coordinates
     * @param Dcal optional 2*5 Jacobian wrpt Cal3_S2 parameters
     * @param Dp optional 2*2 Jacobian wrpt intrinsic coordinates
     * @return point in intrinsic coordinates
     */
    minivector calibrate(const minivector& p) const;
    minivector calibrate(const minivector& p, minimatrix* Dcal,
                         minimatrix* Dp) const;


    /// "Between", subtracts calibrations. between(p,q) == compose(inverse(p),q)

    virtual minimatrix between(const minimatrix* q) const
    {
         return Cal3_S2(minimatrix_get(q,0,0)-fx(), minimatrix_get(q,1,0)-fy(),
                         minimatrix_get(q,2,0)-skew(), minimatrix_get(q,3,0)-px(), minimatrix_get(q,4,0)-py());
    }

    virtual minimatrix between(const minimatrix* q,minimatrix& H1,minimatrix& H2) const
    {
            minimatrix_resize(&H1,dimension,dimension);
            minimatrix_set_neg_identity(&H1);
            minimatrix_resize(&H2,dimension,dimension);
            minimatrix_set_identity(&H2);
        return Cal3_S2(minimatrix_get(q,0,0)-fx(), minimatrix_get(q,1,0)-fy(),
                         minimatrix_get(q,2,0)-skew(), minimatrix_get(q,3,0)-px(), minimatrix_get(q,4,0)-py());

    }


    /// @}
    /// @name Manifold
    /// @{

    /// return DOF, dimensionality of tangent space
    inline int dim() const
    {
        return 5;
    }

    virtual minimatrix* Retract(const minimatrix* d) const
    {
        return new Cal3_S2(fx() + d->data[0], fy() +  d->data[d->prd],
                           skew()+  d->data[2*d->prd], px() +  d->data[3*d->prd], py() +  d->data[4*d->prd]);
    }

        virtual minimatrix LocalCoordinates(const minimatrix* T2) const;
    /// @}

};
};
#endif // CAL3_S2_h
