#ifndef POSE2_H
#define POSE2_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file  Pose2.h
 * @brief 2D Pose
 * @author: Frank Dellaert
 * @author: Richard Roberts
 */


#pragma once

#include "../geometry/Rot2.h"
#include "../base/Matrix.h"
#include "../base/MatCal.h"
#include <iostream>

namespace minisam
{

/**
 * A 2D pose (Point2,Rot2)
 * @addtogroup geometry
 * \nosubgrouping
 */
class Pose2
{
private:
    Rot2 r_;
    Eigen::Vector2d t_;

public:
    /// @name Standard Constructors
    /// @{

    /** default constructor = origin */
    Pose2();

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
    Pose2(double theta, const Eigen::Vector2d &t);

    /** construct from r,t */
    Pose2(const Rot2 &r, const Eigen::Vector2d &t);
    /** Constructor from 3*3 matrix */
    Pose2(const Eigen::Matrix3d &T);

    /// @}
    /// @name Advanced Constructors
    /// @{

    /** Construct from canonical coordinates \f$ [T_x,T_y,\theta] \f$ (Lie algebra) */
    Pose2(const Eigen::VectorXd &v);

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

    /// compose syntactic sugar
    inline Pose2 operator*(const Pose2 &p2) const
    {
        return Pose2(r_ * p2.r(), t_ + r_ * p2.t());
    }
    inline Pose2* multiply(const Pose2 &p2) const
    {
        return new Pose2(r_ * p2.r(), t_ + r_ * p2.t());
    }

    /// @}
    /// @name Lie Group
    /// @{

    ///Exponential map at identity - create a rotation from canonical coordinates \f$ [T_x,T_y,\theta] \f$
    static Pose2 Expmap(const Eigen::Vector3d &xi);
    static Pose2 Expmap(const Eigen::Vector3d &xi, Eigen::MatrixXd *H);

    ///Log map at identity - return the canonical coordinates \f$ [T_x,T_y,\theta] \f$ of this rotation
    static Eigen::Vector3d Logmap(const Pose2 &p);
    static Eigen::Vector3d Logmap(const Pose2 &p, Eigen::MatrixXd *H);
    /**
    * Calculate Adjoint map
    * Ad_pose is 3*3 matrix that when applied to twist xi \f$ [T_x,T_y,\theta] \f$, returns Ad_pose(xi)
    */
    Eigen::Matrix3d AdjointMap() const;
    inline Eigen::Vector3d Adjoint(const Eigen::Vector3d &xi) const
    {
        return AdjointMap() * xi;
    }

    /**
    * Compute the [ad(w,v)] operator for SE2 as in [Kobilarov09siggraph], pg 19
    */
    static Eigen::Matrix3d adjointMap(const Eigen::Vector3d &v);

    /**
    * wedge for SE(2):
    * @param xi 3-dim twist (v,omega) where
    *  omega is angular velocity
    *  v (vx,vy) = 2D velocity
    * @return xihat, 3*3 element of Lie algebra that can be exponentiated
    */
    static inline Eigen::Matrix3d wedge(double vx, double vy, double w)
    {
        Eigen::Matrix3d m;
        m << 0., -w, vx,
        w, 0., vy,
        0., 0., 0.;
        return m;
    }

    /// Derivative of Expmap
    static Eigen::Matrix3d ExpmapDerivative(const Eigen::Vector3d &v);

    /// Derivative of Logmap
    static Eigen::Matrix3d LogmapDerivative(const Pose2 &v);

    // Chart at origin, depends on compile-time flag SLOW_BUT_CORRECT_EXPMAP
    struct ChartAtOrigin
    {
        static Pose2 retract(const Eigen::Vector3d &v);
        static Pose2 retract(const Eigen::Vector3d &v, Eigen::MatrixXd *H);
        static Eigen::VectorXd Local(const Pose2 &r);
        static Eigen::VectorXd Local(const Pose2 &r, Eigen::MatrixXd *H);
    };

    Pose2 retract(const Eigen::VectorXd &v);
    Pose2* retractpointer(const Eigen::VectorXd &v);


    Pose2 between(const Pose2 &X2, Eigen::MatrixXd &H1, Eigen::MatrixXd &H2) const;
    Pose2 between(const Pose2 &X2) const;
    Eigen::VectorXd LocalCoordinates(const Pose2 &pg) const;
    Eigen::VectorXd LocalCoordinates(const Pose2 *pg) const;

    /** Return point coordinates in pose coordinate frame */
    Eigen::Vector2d transform_to(const Eigen::Vector2d &point,
                                 Eigen::MatrixXd *H1,
                                 Eigen::Matrix2d *H2) const;

    Eigen::Vector2d transform_to(const Eigen::Vector2d &point) const;
    /** Return point coordinates in global frame */
    Eigen::Vector2d transform_from(const Eigen::Vector2d &point,
                                   Eigen::MatrixXd *H1,
                                   Eigen::Matrix2d *H2) const;
    Eigen::Vector2d transform_from(const Eigen::Vector2d &point) const;

    /** syntactic sugar for transform_from */
    inline Eigen::Vector2d operator*(const Eigen::Vector2d &point) const
    {
        return transform_from(point);
    }

    /// @}
    /// @name Standard Interface
    /// @{

    /// get x
    inline double x() const
    {
        return t_(0);
    }

    /// get y
    inline double y() const
    {
        return t_(1);
    }

    /// get theta
    inline double theta() const
    {
        return r_.theta();
    }

    /// translation
    inline const Eigen::Vector2d &t() const
    {
        return t_;
    }

    /// rotation
    inline const Rot2 &r() const
    {
        return r_;
    }

    /// translation
    inline const Eigen::Vector2d &translation() const
    {
        return t_;
    }

    /// rotation
    inline const Rot2 &rotation() const
    {
        return r_;
    }

    //// return transformation matrix
    Eigen::Matrix3d matrix() const;

    /**
    * Calculate bearing to a landmark
    * @param point 2D location of landmark
    * @return 2D rotation \f$ \in SO(2) \f$
    */
    Rot2 bearing(const Eigen::Vector2d &point) const;
    Rot2 bearing(const Eigen::Vector2d &point,
                 Eigen::MatrixXd *H1,
                 Eigen::MatrixXd *H2) const;

    /**
    * Calculate bearing to another pose
    * @param point SO(2) location of other pose
    * @return 2D rotation \f$ \in SO(2) \f$
    */
    Rot2 bearing(const Pose2 &pose) const;
    Rot2 bearing(const Pose2 &pose,
                 Eigen::MatrixXd *H1,
                 Eigen::MatrixXd *H2) const;

    /**
    * Calculate range to a landmark
    * @param point 2D location of landmark
    * @return range (double)
    */
    double range(const Eigen::Vector2d &point) const;
    double range(const Eigen::Vector2d &point,
                 Eigen::MatrixXd *H1,
                 Eigen::MatrixXd *H2) const;

    /**
    * Calculate range to another pose
    * @param point 2D location of other pose
    * @return range (double)

    */
    double range(const Pose2 &point) const;
    double range(const Pose2 &point,
                 Eigen::MatrixXd *H1,
                 Eigen::MatrixXd *H2) const;

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

/**
 * Calculate pose between a vector of 2D point correspondences (p,q)
 * where q = Pose2::transform_from(p) = t + R*p
 */
typedef std::pair<Eigen::Vector2d, Eigen::Vector2d> Point2Pair;
void align(const std::vector<Point2Pair> &pairs, Pose2 *alignresult);

std::map<int, Pose2> Pose2ValuesRetract(const std::map<int, Pose2> &LinPose2,
                                        const std::map<int, Eigen::VectorXd> &vectorvalues2);

std::map<int, Eigen::VectorXd> Pose2VectorValuesZero(std::map<int, Pose2*> vvz);

void Pose2ValuesUpdate(std::map<int, Pose2*> *vectorvalues1, const std::map<int, Pose2*> &vectorvalues2);
};
#endif // POSE2_H
