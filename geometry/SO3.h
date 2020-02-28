#ifndef SO3_H
#define SO3_H


/**
 * @file    SO3.h
 * @brief   3*3 matrix representation of SO(3)
 * @author
 * @date
 */

#pragma once

#include "../mat/Matrix.h"
#include "../mat/MatCal.h"

#include <cmath>
#include <iosfwd>

namespace minisam
{

/**
 *  True SO(3), i.e., 3*3 matrix subgroup
 *  We guarantee (all but first) constructors only generate from sub-manifold.
 *  However, round-off errors in repeated composition could move off it...
 */
class SO3:public minimatrix
{
public:

    /// @name Constructors
    /// @{

    /// Constructor from AngleAxisd

    SO3() :minimatrix(3,3)
    {
        minimatrix_set_identity(this);
        dimension=3;
    }


    /// Constructor from Eigen Matrix
    SO3(const minimatrix& R):minimatrix(3,3)
    {
        minimatrix_memcpy(this,R);
        dimension=3;
    }

    SO3(const SO3& obj):minimatrix(3,3)
    {
        minimatrix_memcpy(this,obj);
        dimension=3;
    }



    ~SO3()
    {

    }


    /// Static, named constructor TODO think about relation with above
    static SO3 AxisAngle(const minivector& axis, double theta);

    /// @}

    /// @name Group
    /// @{

    /// identity rotation for group operation
    static SO3 identity();

    /// inverse of a rotation = transpose
    SO3 inverse() const;

    /// @}
    /// @name Lie Group
    /// @{

    /**
     * Exponential map at identity - create a rotation from canonical coordinates
     * \f$ [R_x,R_y,R_z] \f$ using Rodrigues' formula
     */
    static SO3 Expmap(const minivector& omega);
    static SO3 Expmap(const minivector& omega,minimatrix* H);
    /// Derivative of Expmap
    static minimatrix ExpmapDerivative(const minivector& omega);

    /**
     * Log map at identity - returns the canonical coordinates
     * \f$ [R_x,R_y,R_z] \f$ of this rotation
     */
    static minivector Logmap(const SO3& R);//, ChartJacobian H = boost::none);
    static minivector Logmap(const SO3& R, minimatrix* H);
    /// Derivative of Logmap
    static minimatrix LogmapDerivative(const minivector& omega);

    minimatrix AdjointMap() const
    {
        return *this;
    }

       minivector multiplyvector(const minivector& p) const
    {

        minivector result(3);
        miniblas_dgemv (blasNoTrans, 1.0, *this, p, 0.0, &result);
        return result;
    }

    SO3& operator=(const SO3& obj)
    {
        minimatrix_memcpy(this,obj);
        return *this;

    }


     minimatrix multiplymatrix(const  minimatrix& p) const
    {
        minimatrix result(3,3);
        miniblas_dgemm (blasNoTrans,blasNoTrans, 1.0, *this, p, 0.0, &result);
        return result;
    }


    minimatrix* multiplymatrix_pointer(const  minimatrix& p) const
    {
        minimatrix* result=new minimatrix(3,3);
        miniblas_dgemm (blasNoTrans,blasNoTrans, 1.0, *this, p, 0.0, result);
        return result;
    }

     minimatrix scaledouble(double p) const
    {
        minimatrix result(3,3);
        minimatrix_memcpy(&result,*this);
        minimatrix_scale(&result,p);

        return result;
    }


    // Chart at origin
    struct ChartAtOrigin
    {
        static SO3 retract(const minivector& omega)
        {
            return Expmap(omega);
        }
        static SO3 retract(const minivector& omega, minimatrix* H)
        {
            return Expmap(omega, H);
        }
        static minivector Local(const SO3& R)
        {
            return Logmap(R);
        }
        static minivector Local(const SO3& R, minimatrix* H)
        {
            return Logmap(R, H);
        }
    };
    virtual minimatrix* Retract(const minimatrix* mpose);
    virtual minimatrix LocalCoordinates(const minimatrix* mpose,minimatrix* H1=NULL,minimatrix* H2=NULL) const;
    virtual minimatrix between(const minimatrix* mpose) const;
    virtual minimatrix between(const minimatrix* mpose,minimatrix& H1,minimatrix& H2) const;
    /// @}
};

// This namespace exposes two functors that allow for saving computation when
// exponential map and its derivatives are needed at the same location in so<3>
// The second functor also implements dedicated methods to apply dexp and/or inv(dexp)

/// Functor implementing Exponential map
class ExpmapFunctor
{
protected:
    double theta2;
    minimatrix W, K, KK;
    bool nearZero;
    double theta, sin_theta, one_minus_cos;  // only defined if !nearZero

    void init(bool nearZeroApprox = false);

public:
    /// Constructor with element of Lie algebra so(3)
    ExpmapFunctor(double theta_2,const minivector& omega, bool nearZeroApprox = false);

    /// Constructor with axis-angle
    ExpmapFunctor(const minivector& axis, double angle, bool nearZeroApprox = false);

    ~ExpmapFunctor()
    {
    }

    /// Rodrigues formula
    SO3 expmap() const;
};

/// Functor that implements Exponential map *and* its derivatives
class DexpFunctor : public ExpmapFunctor
{
    minivector omega;
    double a, b;
    minimatrix dexp_;

public:
    /// Constructor with element of Lie algebra so(3)
    DexpFunctor(const minivector& omega, bool nearZeroApprox = false);

    // NOTE(luca): Right Jacobian for Exponential map in SO(3) - equation
    // (10.86) and following equations in G.S. Chirikjian, "Stochastic Models,
    // Information Theory, and Lie Groups", Volume 2, 2008.
    //   expmap(omega + v) \approx expmap(omega) * expmap(dexp * v)
    // This maps a perturbation v in the tangent space to
    // a perturbation on the manifold Expmap(dexp * v) */
    minimatrix dexp() const
    {
        return dexp_;
    }
    ~DexpFunctor()
    {
    }

    /// Multiplies with dexp(), with optional derivatives
    minivector applyDexp(const minivector& v) const;
    minivector applyDexp(const minivector& v,minimatrix* H1) const;
    minivector applyDexp(const minivector& v,minimatrix* H1,
                         minimatrix* H2) const;

    /// Multiplies with dexp().inverse(), with optional derivatives
    minivector applyInvDexp(const minivector& v) const;
    minivector applyInvDexp(const minivector& v,
                            minimatrix* H1,
                            minimatrix* H2 ) const;
};
};

#endif // SO3_H
