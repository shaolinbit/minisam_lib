#ifndef  POSE3_H
#define  POSE3_H


/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *@file  Pose3.h
 *@brief 3D Pose
 */

#pragma once
#include "../geometry/Rot3.h"
#include "../gmfconfig.h"

namespace minisam
{

class Pose2;
// forward declare

/**
 * A 3D pose (R,t) : (Rot3,Point3)
 * @addtogroup geometry
 * \nosubgrouping
 */
class  Pose3
{
public:

    Rot3 R_; ///< Rotation gRp, between global and pose frame
    Eigen::Vector3d t_; ///< Translation gTp, from global origin to pose frame origin

public:

    /// @name Standard Constructors
    /// @{

    /** Default constructor is origin */
     Pose3();

    /** Copy constructor */
     Pose3(const Pose3& pose);

    /** Construct from R,t */
     Pose3(const Rot3& R, const Eigen::Vector3d& t);

    /** Construct from Pose2 */
     explicit Pose3(const Pose2& pose2);

    /** Constructor from 4*4 matrix */
     Pose3(const Eigen::Matrix3d &T);

    /// Named constructor with derivatives
    static Pose3 Create(const Rot3& R, const Eigen::Vector3d& t);
    static Pose3 Create(const Rot3& R, const Eigen::Vector3d& t,
                        Eigen::MatrixXd* H1,
                        Eigen::MatrixXd* H2);
    /**
     *  Create Pose3 by aligning two point pairs
     *  A pose aTb is estimated between pairs (a_point, b_point) such that a_point = aTb * b_point
     *  Meant to replace the deprecated function 'align', which orders the pairs the opposite way.
     *  Note this allows for noise on the points but in that case the mapping will not be exact.
     */
    static void Align(const std::vector<std::pair<Eigen::Vector3d,
                      Eigen::Vector3d>>& abPointPairs,Pose3* alignresult);

    /// assert equality up to a tolerance
    bool equals(const Pose3& pose, double tol = 1e-9) const;

    /// @}
    /// @name Group
    /// @{

    /// identity for group operation
    static Pose3 identity();

    /// inverse transformation with derivatives
    Pose3 inverse() const;

    /// compose syntactic sugar
     Pose3 operator*(const Pose3& T) const
    {
        return Pose3(R_ * T.R_, t_ + R_ * T.t_);
    }
    Pose3* multiply(const Pose3& T) const;

    void setT(const Eigen::Vector3d& outt);
    void setRot(const Rot3& outR);


    /// @}
    /// @name Lie Group
    /// @{

    /// Exponential map at identity - create a rotation from canonical coordinates \f$ [R_x,R_y,R_z,T_x,T_y,T_z] \f$
    static Pose3 Expmap(const Eigen::VectorXd& xi);
    static Pose3 Expmap(const Eigen::VectorXd& xi, Eigen::MatrixXd* H);

    /// Log map at identity - return the canonical coordinates \f$ [R_x,R_y,R_z,T_x,T_y,T_z] \f$ of this rotation
    static Eigen::VectorXd Logmap(const Pose3& p);
    static Eigen::VectorXd Logmap(const Pose3& p, Eigen::MatrixXd* H);
    /**
     * Calculate Adjoint map, transforming a twist in the this pose's (i.e, body) frame to the world spatial frame
     * Ad_pose is 6*6 matrix that when applied to twist xi \f$ [R_x,R_y,R_z,T_x,T_y,T_z] \f$, returns Ad_pose(xi)
     */
    Eigen::MatrixXd AdjointMap() const; /// FIXME Not tested - marked as incorrect

    /**
     * Apply this pose's AdjointMap Ad_g to a twist \f$ \xi_b \f$, i.e. a body-fixed velocity, transforming it to the spatial frame
     * \f$ \xi^s = g*\xi^b*g^{-1} = Ad_g * \xi^b \f$
     */
    Eigen::VectorXd Adjoint(const Eigen::VectorXd& xi_b) const
    {
        return AdjointMap() * xi_b;
    } /// FIXME Not tested - marked as incorrect

    /**
     * Compute the [ad(w,v)] operator as defined in [Kobilarov09siggraph], pg 11
     * [ad(w,v)] = [w^, zero3; v^, w^]
     * Note that this is the matrix representation of the adjoint operator for se3 Lie algebra,
     * aka the Lie bracket, and also the derivative of Adjoint map for the Lie group SE3.
     *
     * Let \f$ \hat{\xi}_i \f$ be the se3 Lie algebra, and \f$ \hat{\xi}_i^\vee = \xi_i = [\omega_i,v_i] \in \mathbb{R}^6\f$ be its
     * vector representation.
     * We have the following relationship:
     * \f$ [\hat{\xi}_1,\hat{\xi}_2]^\vee = ad_{\xi_1}(\xi_2) = [ad_{(\omega_1,v_1)}]*\xi_2 \f$
     *
     * We use this to compute the discrete version of the inverse right-trivialized tangent map,
     * and its inverse transpose in the discrete Euler Poincare' (DEP) operator.
     *
     */
    static Eigen::MatrixXd adjointMap(const Eigen::VectorXd  &xi);

    /**
     * Action of the adjointMap on a Lie-algebra vector y, with optional derivatives
     */
    static Eigen::VectorXd  adjoint(const Eigen::VectorXd  &xi, const Eigen::VectorXd  &y);
    static Eigen::VectorXd  adjoint(const Eigen::VectorXd  &xi, const Eigen::VectorXd  &y,
                                    Eigen::MatrixXd* H);
    /**
     * The dual version of adjoint action, acting on the dual space of the Lie-algebra vector space.
     */
    static Eigen::VectorXd adjointTranspose(const Eigen::VectorXd& xi, const Eigen::VectorXd& y);

    static Eigen::VectorXd adjointTranspose(const Eigen::VectorXd& xi, const Eigen::VectorXd& y,
                                            Eigen::MatrixXd* H);

    /// Derivative of Expmap
    static Eigen::MatrixXd ExpmapDerivative(const Eigen::VectorXd& xi);

    /// Derivative of Logmap
    static Eigen::MatrixXd LogmapDerivative(const Pose3& xi);

    // Chart at origin, depends on compile-time flag POSE3_EXPMAP
    struct ChartAtOrigin
    {
        static Pose3 retract(const Eigen::VectorXd& v);
        static Pose3 retract(const Eigen::VectorXd& v, Eigen::MatrixXd* H);
        static Eigen::VectorXd Local(const Pose3& r);
        static Eigen::VectorXd Local(const Pose3& r, Eigen::MatrixXd*  H);
    };


    Pose3 retract(const Eigen::VectorXd& v);

    Pose3* retractpointer(const Eigen::VectorXd& v);

    Eigen::VectorXd LocalCoordinates(const Pose3& pg) const;

    Eigen::VectorXd LocalCoordinates(const Pose3* pg) const;

    Pose3 compose(const Pose3& g, Eigen::MatrixXd* H1=NULL,
                  Eigen::MatrixXd* H2 = NULL) const;
    Pose3 between(const Pose3& X2,Eigen::MatrixXd& H1,Eigen::MatrixXd& H2) const;
    Pose3 between(const Pose3& X2) const;
    /**
     * wedge for Pose3:
     * @param xi 6-dim twist (omega,v) where
     *  omega = (wx,wy,wz) 3D angular velocity
     *  v (vx,vy,vz) = 3D velocity
     * @return xihat, 4*4 element of Lie algebra that can be exponentiated
     */
    static Eigen::MatrixXd wedge(double wx, double wy, double wz, double vx, double vy,
                                 double vz);

    /// @}
    /// @name Group Action on Point3
    /// @{

    /**
     * @brief takes point in Pose coordinates and transforms it to world coordinates
     * @param p point in Pose coordinates
     * @param Dpose optional 3*6 Jacobian wrpt this pose
     * @param Dpoint optional 3*3 Jacobian wrpt point
     * @return point in world coordinates
     */
    Eigen::Vector3d transform_from(const Eigen::Vector3d& p) const;
    Eigen::Vector3d transform_from(const Eigen::Vector3d& p,Eigen::MatrixXd* Dpose,
                                   Eigen::Matrix3d* Dpoint) const;

    /** syntactic sugar for transform_from */
     inline Eigen::Vector3d operator*(const Eigen::Vector3d& p) const
    {
        return transform_from(p);
    }

    /**
     * @brief takes point in world coordinates and transforms it to Pose coordinates
     * @param p point in world coordinates
     * @param Dpose optional 3*6 Jacobian wrpt this pose
     * @param Dpoint optional 3*3 Jacobian wrpt point
     * @return point in Pose coordinates
     */
    Eigen::Vector3d transform_to(const Eigen::Vector3d& p) const;

    Eigen::Vector3d transform_to(const Eigen::Vector3d& p, Eigen::MatrixXd* Dpose,
                                 Eigen::MatrixXd* Dpoint) const;

    Eigen::Vector3d transform_toDPoint(const Eigen::Vector3d& p,
                                       Eigen::Matrix3d* Dpoint) const;
    /// @}
    /// @name Standard Interface
    /// @{

    /// get rotation
    const Rot3& rotation() const;
    const Rot3& rotation(Eigen::MatrixXd* H) const;
    /// get translation
    const Eigen::Vector3d& translation() const;
    const Eigen::Vector3d& translation(Eigen::MatrixXd*  H) const;
    /// get x
    double x() const;

    /// get y
    double y() const;

    /// get z
    double z() const;

    /** convert to 4*4 matrix */
    Eigen::MatrixXd matrix() const;

    /** receives a pose in local coordinates and transforms it to world coordinates
    * @deprecated: This is actually equivalent to transform_from, so it is WRONG! Use
    * transform_pose_to instead. */
    Pose3 transform_to(const Pose3& pose) const;

    /** receives a pose in world coordinates and transforms it to local coordinates */
    Pose3 transform_pose_to(const Pose3& pose) const;
    Pose3 transform_pose_to(const Pose3& pose, Eigen::MatrixXd* H1,
                            Eigen::MatrixXd* H2) const;

    /**
     * Calculate range to a landmark
     * @param point 3D location of landmark
     * @return range (double)
     */
    double range(const Eigen::Vector3d& point) const;
    double range(const Eigen::Vector3d& point, Eigen::MatrixXd* H1,
                 Eigen::MatrixXd* H2) const;
    /**
     * Calculate range to another pose
     * @param pose Other SO(3) pose
     * @return range (double)
     */
    double range(const Pose3& pose) const;
    double range(const Pose3& pose, Eigen::MatrixXd* H1,
                 Eigen::MatrixXd* H2) const;

    /**
     * Calculate bearing to a landmark
     * @param point 3D location of landmark
     * @return bearing (Unit3)
     */
    Unit3 bearing(const Eigen::Vector3d& point) const;
    Unit3 bearing(const Eigen::Vector3d& point, Eigen::MatrixXd* H1,
                  Eigen::MatrixXd* H2) const;
    /**
     * Calculate bearing to another pose
     * @param other 3D location and orientation of other body. The orientation
     * information is ignored.
     * @return bearing (Unit3)
     */
    Unit3 bearing(const Pose3& pose) const;

    Unit3 bearing(const Pose3& pose, Eigen::MatrixXd* H1,
                  Eigen::MatrixXd* H2) const;
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
        return std::make_pair(3, 5);
    }

    /**
     * Return the start and end indices (inclusive) of the rotation component of the
     * exponential map parameterization
     * @return a pair of [start, end] indices into the tangent space vector
     */
    static std::pair<int, int> rotationInterval()
    {
        return std::make_pair(0, 2);
    }

    /// Output stream operator

     friend std::ostream &operator<<(std::ostream &os, const Pose3& p);
    /// @}

};
// Pose3 class

/**
 * wedge for Pose3:
 * @param xi 6-dim twist (omega,v) where
 *  omega = 3D angular velocity
 *  v = 3D velocity
 * @return xihat, 4*4 element of Lie algebra that can be exponentiated
 */
inline Eigen::MatrixXd wedge(const Eigen::VectorXd& xi)
{
    return Pose3::wedge(xi(0), xi(1), xi(2), xi(3), xi(4), xi(5));
}

/**
 * Calculate pose between a vector of 3D point correspondences (b_point, a_point)
 * where a_point = Pose3::transform_from(b_point) = t + R*b_point
 * @deprecated: use Pose3::Align with point pairs ordered the opposite way
 */
void align(const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& baPointPairs,Pose3* alignresult);
std::map<int,Pose3> Pose3ValuesRetract(const std::map<int,Pose3>& LinPose3,
                                       const std::map<int,Eigen::VectorXd>& vectorvalues2);
std::map<int,Pose3> DPose3ValuesRetract(const std::map<int,Pose3>& LinPose3,
                                        const std::map<int,Eigen::VectorXd>& vectorvalues2);

std::map<int,Eigen::VectorXd> Pose3VectorValuesZero(std::map<int,Pose3> vvz);
void Pose3ValuesUpdate(std::map<int,Pose3*>* vectorvalues1,const std::map<int,Pose3*>& vectorvalues2);

};
#endif // POSE3_H
