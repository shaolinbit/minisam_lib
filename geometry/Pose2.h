#ifndef POSE2_H
#define POSE2_H


/**
 * @file  Pose2.h
 * @brief 2D Pose
 */


#pragma once

#include "../geometry/Rot2.h"
#include "../mat/Matrix.h"
#include "../mat/MatCal.h"
#include <iostream>

using namespace std;
namespace minisam
{

/**
 * A 2D pose (Point2,Rot2)
 * @addtogroup geometry
 * \nosubgrouping
 */
class Pose2:public minivector
{
public:
    /// @name Standard Constructors
    /// @{

    /** default constructor = origin */
    Pose2();

    ~Pose2();

    /** copy constructor */
    Pose2(const Pose2 &pose);
    /**
    * construct from (x,y,theta)
    * @param x x coordinate
    * @param y y coordinate
    * @param theta angle with positive X-axis
    */
    Pose2(double x, double y, double theta);

    /** construct from rotation and translation */
    Pose2(double theta, const minivector &t);

    /** construct from r,t */
    Pose2(const Rot2 &r, const minivector &t);
    /** Constructor from 3*3 matrix */
    Pose2(const minimatrix &T);

    Pose2(const minimatrix* T);

    /// @}
    /// @name Advanced Constructors
    /// @{

    /** Construct from canonical coordinates \f$ [T_x,T_y,\theta] \f$ (Lie algebra) */
    Pose2(const minivector &v);

    /// @}
    /// @name Group
    /// @{

    /// identity for group operation
    inline static Pose2 identity()
    {
        return Pose2();
    }

    /// inverse
    Pose2 inverse() const;


    inline Pose2 multiply(const Pose2 &p2) const
    {
        minivector b=r().multiplyvector(p2.t());
        minivector_add(&b,t());
        Pose2 result(r().multiply(p2.r()),b);

        if(DEBUGSTATE)
        {
        cout<<"this and pose2"<<endl;
        cout<<*this<<endl;
        cout<<p2<<endl;

        std::cout<<"result"<<std::endl;
        std::cout<<result<<std::endl;
        }
        return result;
    }

    inline Pose2* multiplypointer(const Pose2 &p2) const
    {
        minivector b=r().multiplyvector(p2.t());
        minivector_add(&b,t());
        Pose2* result=new Pose2(r().multiply(p2.r()),b);
        return result;
    }

    /// @}
    /// @name Lie Group
    /// @{

    ///Exponential map at identity - create a rotation from canonical coordinates \f$ [T_x,T_y,\theta] \f$
    static Pose2 Expmap(const minivector &xi);
    static Pose2 Expmap(const minivector &xi, minimatrix *H);

    ///Log map at identity - return the canonical coordinates \f$ [T_x,T_y,\theta] \f$ of this rotation
    static minivector Logmap(const Pose2 &p);
    static minivector Logmap(const Pose2 &p, minimatrix *H);
    /**
    * Calculate Adjoint map
    * Ad_pose is 3*3 matrix that when applied to twist xi \f$ [T_x,T_y,\theta] \f$, returns Ad_pose(xi)
    */
    minimatrix AdjointMap() const;
    inline minivector Adjoint(const minivector &xi) const;


    /**
    * Compute the [ad(w,v)] operator for SE2 as in [Kobilarov09siggraph], pg 19
    */
    static minimatrix adjointMap(const minivector &v);

    /**
    * wedge for SE(2):
    * @param xi 3-dim twist (v,omega) where
    *  omega is angular velocity
    *  v (vx,vy) = 2D velocity
    * @return xihat, 3*3 element of Lie algebra that can be exponentiated
    */
    static inline minimatrix wedge(double vx, double vy, double w)
    {
        minimatrix m(3,3);
        m.data[0]=0;
        m.data[1]= -w;
        m.data[2]=vx;
        m.data[3]=w;
        m.data[4]=0.0;
        m.data[5]=vy;
        m.data[6]=0.0;
        m.data[7]= 0.0;
        m.data[8]=0.0;
        return m;
    }

    Pose2* compose(const Pose2& g, minimatrix* H1=NULL,
                  minimatrix* H2 = NULL) const;

    /// Derivative of Expmap
    static minimatrix ExpmapDerivative(const minivector &v);

    /// Derivative of Logmap
    static minimatrix LogmapDerivative(const Pose2 &v);

    // Chart at origin, depends on compile-time flag SLOW_BUT_CORRECT_EXPMAP
    struct ChartAtOrigin
    {
        static Pose2 retract(const minivector &v);
        static Pose2 retract(const minivector &v, minimatrix *H);
        static minivector Local(const Pose2 &r);
        static minivector Local(const Pose2 &r, minimatrix *H);
    };

    virtual minimatrix* Retract(const minimatrix* mpose);
    virtual minimatrix  between(const minimatrix* mpose) const ;
    virtual minimatrix between(const minimatrix* mpose,minimatrix& H1,minimatrix& H2) const;
    virtual minimatrix LocalCoordinates(const minimatrix *pg,minimatrix* H1=NULL,minimatrix* H2=NULL) const;

    /** Return point coordinates in pose coordinate frame */
    minivector transform_to(const minivector &point,
                            minimatrix *H1,
                            minimatrix *H2) const;

    minivector transform_to(const minivector &point) const;
    /** Return point coordinates in global frame */
    minivector transform_from(const minivector &point,
                              minimatrix *H1,
                              minimatrix *H2) const;
    minivector transform_from(const minivector &point) const;

    /** syntactic sugar for transform_from */
    minivector multiplyvector(const minivector &point) const
    {
        return transform_from(point);
    }

    /// @}
    /// @name Standard Interface
    /// @{

    /// get x
    virtual double x() const
    {
        return minivector_get(this,2);
    }

    /// get y
    virtual double y() const
    {
        return minivector_get(this,3);
    }

    /// get theta
    inline double theta() const
    {
        return std::atan2(data[1], data[0]);
    }

    /// translation
    inline minivector t() const
    {
        return minivector_subvector(*this,2,2);
    }

    /// rotation
    inline Rot2 r() const
    {
        return Rot2(minivector_get(this,0),minivector_get(this,1));
    }

    /// translation
    inline minivector translation() const
    {
        return minivector_subvector(*this,2,2);
    }

    /// rotation
    inline Rot2 rotation() const
    {
        minivector bb=minivector_subvector(*this,0,2);
        return Rot2(bb);
    }

    //// return transformation matrix
    minimatrix matrix() const;

    /**
    * Calculate bearing to a landmark
    * @param point 2D location of landmark
    * @return 2D rotation \f$ \in SO(2) \f$
    */
    Rot2 bearing(const minivector &point) const;
    Rot2 bearing(const minivector &point,
                 minimatrix *H1,
                 minimatrix *H2) const;

    /**
    * Calculate bearing to another pose
    * @param point SO(2) location of other pose
    * @return 2D rotation \f$ \in SO(2) \f$
    */
    Rot2 bearing(const Pose2 &pose) const;
    Rot2 bearing(const Pose2 &pose,
                 minimatrix *H1,
                 minimatrix *H2) const;

    /**
    * Calculate range to a landmark
    * @param point 2D location of landmark
    * @return range (double)
    */
    double range(const minivector &point) const;
    double range(const minivector &point,
                 minimatrix *H1,
                 minimatrix *H2) const;

    /**
    * Calculate range to another pose
    * @param point 2D location of other pose
    * @return range (double)

    */
    double range(const Pose2 &point) const;
    double range(const Pose2 &point,
                 minimatrix *H1,
                 minimatrix *H2) const;

    /// @}
    /// @name Advanced Interface
    /// @{

    /**
    * Return the start and end indices (inclusive) of the translation component of the
    * exponential map parameterization
    * @return a pair of [start, end] indices into the tangent space vector
    */
    inline static std::pair<int, int> translationInterval()
    {
        return std::make_pair(0, 1);
    }

    /**
    * Return the start and end indices (inclusive) of the rotation component of the
    * exponential map parameterization
    * @return a pair of [start, end] indices into the tangent space vector
    */
    static std::pair<int, int> rotationInterval();

    friend std::ostream &operator<<(std::ostream &os, const Pose2 &p);
    ///@}
}; // Pose2

};
#endif // POSE2_H
