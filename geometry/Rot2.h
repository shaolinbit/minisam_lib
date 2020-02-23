#ifndef  ROT2_H
#define  ROT2_H

/**
 * @file Rot2.h
 * @brief 2D rotation
 * @date
 * @author
 */
#pragma once
#include "../mat/Matrix.h"
#include "../mat/MatCal.h"
#include <assert.h>
/**
 * Rotation matrix
 * NOTE: the angle theta is in radians unless explicitly stated
 * @addtogroup geometry
 * \nosubgrouping
 */
namespace minisam
{

class Rot2:public minivector
{
public:
    /** normalize to make sure cos and sin form unit vector */
    Rot2& normalize();

    /** private constructor from cos/sin */
    inline Rot2(double c, double s) : minivector(2)
    {
        dimension=1;
        data[0]=c;
        data[1]=s;
    }


    /// @name Constructors and named constructors
    /// @{

    /** default constructor, zero rotation */
    Rot2() : minivector(2)
    {
        dimension=1;
        data[0]=1.0;
        data[1]=0.0;
    }

    /// Constructor from angle in radians == exponential map at identity
    Rot2(double theta): minivector(2)
    {
        dimension=1;
        data[0]=cos(theta);
        data[1]=sin(theta);
    }

    Rot2(const Rot2& rr) : minivector(2)
    {
        dimension=1;
        data[0]=rr.data[0];
        data[1]=rr.data[1];
    }

    Rot2(const minivector& rr) : minivector(2)
    {
        assert(rr.size1==2);
        data[0]=minivector_get(&rr,0);
        data[1]=minivector_get(&rr,1);
        dimension=1;
    }
    Rot2(const minivector* m_memory)
    {
        size1=m_memory->size1;
        size2=1;
        prd=m_memory->prd;
        data=m_memory->data;
        owner=0;
        dimension=m_memory->dimension;
    }

     Rot2(const minimatrix* m_memory)
    {
    size1=m_memory->size1;
    size2=m_memory->size2;
    prd=m_memory->prd;
    data=m_memory->data;
    owner=0;
    dimension=1;
    }



    /// Named constructor from angle in radians
    static Rot2 fromAngle(double theta)
    {
        return Rot2(theta);
    }

    /// Named constructor from angle in degrees
    static Rot2 fromDegrees(double theta)
    {
        const double degree = M_PI / 180;
        return fromAngle(theta * degree);
    }

    /// Named constructor from cos(theta),sin(theta) pair, will *not* normalize!
    static Rot2 fromCosSin(double c, double s);

    /**
     * Named constructor with derivative
     * Calculate relative bearing to a landmark in local coordinate frame
     * @param d 2D location of landmark
     * @param H optional reference for Jacobian
     * @return 2D rotation \f$ \in SO(2) \f$
     */
    static Rot2 relativeBearing(const minivector& d,minimatrix* H=NULL);
    /** Named constructor that behaves as atan2, i.e., y,x order (!) and normalizes */
    static Rot2 atan2(double y, double x);

    /// @}
    /// @name Group
    /// @{

    /** identity */
    inline static Rot2 identity()
    {
        return Rot2();
    }

    /** The inverse rotation - negative angle */
    Rot2 inverse() const
    {
        return Rot2(data[0], -data[1]);
    }

     Rot2 multiply(const Rot2& R) const
    {
        //return fromCosSin(c_ * R.c_ - s_ * R.s_, s_ * R.c_ + c_ * R.s_);
        return fromCosSin(data[0] * R.data[0] - data[1] * R.data[1],
        data[1] * R.data[0] + data[0] * R.data[1]);
    }

    /// @}
    /// @name Lie Group
    /// @{

    /// Exponential map at identity - create a rotation from canonical coordinates
    static Rot2 Expmap(const minivector& v);
    Rot2* ExpmapP(const minivector& v);
    static Rot2 Expmap(const minivector& v,minimatrix* H);
    /// Log map at identity - return the canonical coordinates of this rotation
    static double Logmap(const Rot2& r);
    static double Logmap(const Rot2& r, minimatrix* H);
    /** Calculate Adjoint map */
    minimatrix AdjointMap() const
    {
        minimatrix ajm(1,1);
        ajm.data[0]=1.0;
        return ajm;
    }

    /// Left-trivialized derivative of the exponential map
    static minimatrix ExpmapDerivative()
    {
        minimatrix exmd(1,1);
        exmd.data[0]=1.0;
        return exmd;
    }

    /// Left-trivialized derivative inverse of the exponential map
    static minimatrix LogmapDerivative()
    {
        minimatrix logmd(1,1);
        logmd.data[0]=1.0;
        return logmd;
    }

    // Chart at origin simply uses exponential map and its inverse
    struct ChartAtOrigin
    {
        static Rot2 retract(const minivector& v)
        {
            return Expmap(v);
        }
        static Rot2 retract(const minivector& v,minimatrix* H)
        {
            return Expmap(v,H);
        }
        static double Local(const Rot2& r)
        {
            return Logmap(r);
        }
        static double Local(const Rot2& r,minimatrix* H)
        {
            return Logmap(r,H);
        }
    };
    double local(const Rot2&r);
    /// @}
    /// @name Group Action on Point2
    /// @{

    /**
     * rotate point from rotated coordinate frame to world \f$ p^w = R_c^w p^c \f$
     */
    minivector rotate(const minivector& p) const;
    minivector rotate(const minivector& p,minimatrix* H1,
                      minimatrix* H2) const;
    /** syntactic sugar for rotate
    inline minivector operator*(const minivector& p) const
    }*/
    inline minivector multiplyvector(const minivector& p) const
    {
        return rotate(p);
    }

    virtual minimatrix LocalCoordinates(const minimatrix* mpose) const
    {
        Rot2 g(mpose);
        double re= Rot2::ChartAtOrigin::Local(g);
        minivector result(1);
        result.data[0]=re;
        return result;
    }

    virtual minimatrix* Retract(const minimatrix* mpose)
    {
        minivector v(mpose);
        return ExpmapP(v);
    }

    virtual minimatrix between(const minimatrix* mpose) const
    {
        Rot2 X2(mpose);
        Rot2 result=this->inverse().multiply(X2);
        return result;

    }
    virtual minimatrix between(const minimatrix* mpose,minimatrix& H1,minimatrix& H2) const
    {
        Rot2 X2(mpose);
        Rot2 result=this->inverse().multiply(X2);
        minimatrix_resize(&H1,1,1);
        minimatrix_memcpy(&H1,result.inverse().AdjointMap());
        minimatrix_scale(&H1,-1.0);
        minimatrix_resize(&H2,1,1);
        minimatrix_set_identity(&H2);

        return result;

    }


    /**
     * rotate point from world to rotated frame \f$ p^c = (R_c^w)^T p^w \f$
     */

    minivector unrotate(const minivector& p) const;
    minivector unrotate(const minivector& p,
                        minimatrix* H1,
                        minimatrix* H2) const;
    /// @}
    /// @name Standard Interface
    /// @{

    /// Creates a unit vector as a Point2
    inline minivector unit() const
    {
        minivector p2(2);
        p2.data[0]=data[0];
        p2.data[1]=data[1];

        return p2;
    }

    /** return angle (RADIANS) */
    double theta() const
    {
        return std::atan2(data[1], data[0]);
    }

    /** return angle (DEGREES) */
    double degrees() const
    {
        const double degree = M_PI / 180;
        return theta() / degree;
    }

    /** return cos */
    inline double c() const
    {
        return data[0];
    }

    /** return sin */
    inline double s() const
    {
        return data[1];
    }

    /** return 2*2 rotation matrix */
    minimatrix matrix() const;
    /** return 2*2 transpose (inverse) rotation matrix   */
    minimatrix transpose() const;

    friend std::ostream &operator<<(std::ostream &os, const Rot2& p);

    ///@}

};
};
#endif // ROT2_H
