
/**
 * @file StereoPoint2.h
 * @brief A 2D stereo point (uL,uR,v)
 */

#pragma once

#include "../mat/MatCal.h"
#include "../mat/Matrix.h"

namespace minisam
{

/**
 * A 2D stereo point, v will be same for rectified images
 * @addtogroup geometry
 * \nosubgrouping
 */
class StereoPoint2:public minivector
{

public:
    /// @name Standard Constructors
    /// @{

    /** default constructor */
    StereoPoint2() :minivector(3)
    {
        minivector_set_zero(this);
    }

    /** constructor */
    StereoPoint2(double uL, double uR, double v) :minivector(3)//uL_(uL), uR_(uR), v_(v)
    {
        data[0]=uL;
        data[1]=uR;
        data[2]=v;

    }

    /// construct from 3D vector
    explicit StereoPoint2(const minivector& v) :minivector(3)
    {
        minivector_memcpy(this,v);
    }

    /// @}

    /// @name Group
    /// @{

    /// identity
    inline static StereoPoint2 identity()
    {
        return StereoPoint2();
    }

    /// inverse
    StereoPoint2 operator-() const
    {
        return StereoPoint2(-uL(), -uR(), -v());
    }

    /// add vector on right
    inline StereoPoint2 operator +(const minivector& d) const
    {
        return StereoPoint2(uL() +minivector_get(&d,0), uR() + minivector_get(&d,1), v() + minivector_get(&d,2));
    }

    /// add
    inline StereoPoint2 operator +(const StereoPoint2& b) const
    {
        return StereoPoint2(uL() + b.uL(), uR() + b.uR(), v() + b.v());

    }

    /// subtract
    inline StereoPoint2 operator -(const StereoPoint2& b) const
    {
        return StereoPoint2(uL() - b.uL(), uR() - b.uR(), v() - b.v());
    }

    /// @}
    /// @name Standard Interface
    /// @{

    /// equality
    inline bool operator ==(const StereoPoint2& q) const
    {
        return uL()== q.uL() && uR()==q.uR() && v() == q.v();
    }

    /// get uL
    inline double uL() const
    {
        return minivector_get(this,0);
    }

    /// get uR
    inline double uR() const
    {
        return minivector_get(this,1);
    }

    /// get v
    inline double v() const
    {
        return minivector_get(this,2);
    }

    /** convert to vector */
    minivector vector() const
    {
        return minivector(*this);
    }

    /** convenient function to get a Point2 from the left image */
    minivector point2() const
    {
        minivector result(2);
        result.data[0]=uL();
        result.data[1]=v();
        return result;
    }

    /** convenient function to get a Point2 from the right image */
    minivector right() const
    {
        minivector result(2);
        result.data[0]=uR();
        result.data[1]=v();
        return result;
    }


    inline StereoPoint2 inverse() const
    {
        return StereoPoint2()- (*this);
    }
    inline StereoPoint2 compose(const StereoPoint2& p1) const
    {
        return *this + p1;
    }


    virtual  minimatrix between(const minimatrix* mpose) const
    {
        minivector result(3);
        result.data[0]=minimatrix_get(mpose,0,0)-minivector_get(this,0);
        result.data[1]=minimatrix_get(mpose,1,0)-minivector_get(this,1);
        result.data[2]=minimatrix_get(mpose,2,0)-minivector_get(this,2);
        return result;
    }

    virtual  minimatrix between(const minimatrix* mpose,minimatrix& H1,minimatrix& H2) const
    {
        minivector result(3);
        result.data[0]=minimatrix_get(mpose,0,0)-minivector_get(this,0);
        result.data[1]=minimatrix_get(mpose,1,0)-minivector_get(this,1);
        result.data[2]=minimatrix_get(mpose,2,0)-minivector_get(this,2);

        minimatrix_resize(&H1,dimension,dimension);
        minimatrix_set_neg_identity(&H1);
        minimatrix_resize(&H2,dimension,dimension);
        minimatrix_set_identity(&H2);

        return result;
    }

    virtual minimatrix LocalCoordinates(const minimatrix* mpose,minimatrix* H1=NULL,minimatrix* H2=NULL) const
    {
        StereoPoint2 result;
        result.data[0]=minimatrix_get(mpose,0,0)-minivector_get(this,0);
        result.data[1]=minimatrix_get(mpose,1,0)-minivector_get(this,1);
        result.data[2]=minimatrix_get(mpose,2,0)-minivector_get(this,2);
        return result;
    }

    virtual minimatrix* Retract(const minimatrix* mpose)
    {
        // minivector v(mpose);
        minivector composeresult(3);
        composeresult.data[0]=data[0]+minimatrix_get(mpose,0,0);
        composeresult.data[1]=data[1]+minimatrix_get(mpose,1,0);
        composeresult.data[2]=data[2]+minimatrix_get(mpose,2,0);

        return new StereoPoint2(composeresult);
    }

    static inline minivector Logmap(const StereoPoint2& p)
    {
        return p.vector();
    }
    static inline StereoPoint2 Expmap(const minivector& d)
    {
        return StereoPoint2(minivector_get(&d,0),minivector_get(&d,1), minivector_get(&d,2));
    }

};

};
