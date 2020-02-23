#ifndef  POSE3_H
#define  POSE3_H


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
class  Pose3:
#ifdef USE_QUATERNIONS
    public minivector
#else
    public minimatrix
#endif
{
public:

    /// @name Standard Constructors
    /// @{

    /** Default constructor is origin */
    Pose3();
    ~Pose3();

    /** Copy constructor */
    Pose3(const Pose3& pose);
    Pose3(const Pose3* pose);
    Pose3(const minimatrix* pose);

    /** Construct from R,t */
    Pose3(const Rot3& R, const minivector& t);

    /** Construct from Pose2 */
    explicit Pose3(const Pose2& pose2);

    /** Constructor from 3*4 matrix */
    Pose3(const minimatrix &T);

    /// Named constructor with derivatives
    static Pose3 Create(const Rot3& R, const minivector& t);
    static Pose3 Create(const Rot3& R, const minivector& t,
                        minimatrix* H1,
                        minimatrix* H2);

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
    Pose3 multiply(const Pose3& T) const;

    Pose3& operator=(const Pose3& obj);


    Pose3* multiply_pointer(const Pose3& T) const;

    void setT(const minivector& outt);
    void setRot(const Rot3& outR);


    /// @}
    /// @name Lie Group
    /// @{

    /// Exponential map at identity - create a rotation from canonical coordinates \f$ [R_x,R_y,R_z,T_x,T_y,T_z] \f$
    static Pose3 Expmap(const minivector& xi);
    static Pose3 Expmap(const minivector& xi, minimatrix* H);

    /// Log map at identity - return the canonical coordinates \f$ [R_x,R_y,R_z,T_x,T_y,T_z] \f$ of this rotation
    static minivector Logmap(const Pose3& p);
    static minivector Logmap(const Pose3& p, minimatrix* H);
    /**
     * Calculate Adjoint map, transforming a twist in the this pose's (i.e, body) frame to the world spatial frame
     * Ad_pose is 6*6 matrix that when applied to twist xi \f$ [R_x,R_y,R_z,T_x,T_y,T_z] \f$, returns Ad_pose(xi)
     */
    minimatrix AdjointMap() const; /// FIXME Not tested - marked as incorrect

    /**
     * Apply this pose's AdjointMap Ad_g to a twist \f$ \xi_b \f$, i.e. a body-fixed velocity, transforming it to the spatial frame
     * \f$ \xi^s = g*\xi^b*g^{-1} = Ad_g * \xi^b \f$
     */
    minivector Adjoint(const minivector& xi_b) const;
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
    static minimatrix adjointMap(const minivector  &xi);

    /**
     * Action of the adjointMap on a Lie-algebra vector y, with optional derivatives
     */
    static minivector  adjoint(const minivector  &xi, const minivector  &y);
    static minivector  adjoint(const minivector  &xi, const minivector  &y,
                               minimatrix* H);
    /**
     * The dual version of adjoint action, acting on the dual space of the Lie-algebra vector space.
     */
    static minivector adjointTranspose(const minivector& xi, const minivector& y);

    static minivector adjointTranspose(const minivector& xi, const minivector& y,
                                       minimatrix* H);

    /// Derivative of Expmap
    static minimatrix ExpmapDerivative(const minivector& xi);

    /// Derivative of Logmap
    static minimatrix LogmapDerivative(const Pose3& xi);

    // Chart at origin, depends on compile-time flag POSE3_EXPMAP
    struct ChartAtOrigin
    {
        static Pose3 retract(const minivector& v);
        static Pose3 retract(const minivector& v, minimatrix* H);
        static minivector Local(const Pose3& r);
        static minivector Local(const Pose3& r, minimatrix*  H);
    };


    Pose3 retract(const minivector& v);

    virtual minimatrix* Retract(const minimatrix* mpose);

    virtual minimatrix LocalCoordinates(const minimatrix* pg) const;

    Pose3 compose(const Pose3& g, minimatrix* H1=NULL,
                  minimatrix* H2 = NULL) const;
    virtual minimatrix between(const minimatrix* mpose) const;
    virtual minimatrix between(const minimatrix* mpose,minimatrix& H1,minimatrix& H2) const;
    /**
     * wedge for Pose3:
     * @param xi 6-dim twist (omega,v) where
     *  omega = (wx,wy,wz) 3D angular velocity
     *  v (vx,vy,vz) = 3D velocity
     * @return xihat, 4*4 element of Lie algebra that can be exponentiated
     */
    static minimatrix wedge(double wx, double wy, double wz, double vx, double vy,
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
    minivector transform_from(const minivector& p) const;
    minivector transform_from(const minivector& p,minimatrix* Dpose,
                              minimatrix* Dpoint) const;

    /** syntactic sugar for transform_from */
    inline minivector multiplyvector(const minivector& p) const;

    /**
     * @brief takes point in world coordinates and transforms it to Pose coordinates
     * @param p point in world coordinates
     * @param Dpose optional 3*6 Jacobian wrpt this pose
     * @param Dpoint optional 3*3 Jacobian wrpt point
     * @return point in Pose coordinates
     */
    minivector transform_to(const minivector& p) const;

    minivector transform_to(const minivector& p, minimatrix* Dpose,
                            minimatrix* Dpoint) const;

    minivector transform_toDPoint(const minivector& p,
                                  minimatrix* Dpoint) const;
    /// @}
    /// @name Standard Interface
    /// @{

    /// get rotation
    Rot3 rotation() const;
    Rot3 rotation(minimatrix* H) const;
    /// get translation
     minivector translation() const;
     minivector translation(minimatrix*  H) const;

    minivector translationP(minimatrix&  H) const;
    /// get x
    virtual double x() const;

    /// get y
    virtual double y() const;

    /// get z
    double z() const;

    /** convert to 4*4 matrix */
    minimatrix matrix() const;

    /** receives a pose in local coordinates and transforms it to world coordinates
    * @deprecated: This is actually equivalent to transform_from, so it is WRONG! Use
    * transform_pose_to instead. */
    Pose3 transform_to(const Pose3& pose) const;

    /** receives a pose in world coordinates and transforms it to local coordinates */
    Pose3 transform_pose_to(const Pose3& pose) const;
    Pose3 transform_pose_to(const Pose3& pose, minimatrix* H1,
                            minimatrix* H2) const;

    /**
     * Calculate range to a landmark
     * @param point 3D location of landmark
     * @return range (double)
     */
    double range(const minivector& point) const;
    double range(const minivector& point, minimatrix* H1,
                 minimatrix* H2) const;
    /**
     * Calculate range to another pose
     * @param pose Other SO(3) pose
     * @return range (double)
     */
    double range(const Pose3& pose) const;
    double range(const Pose3& pose, minimatrix* H1,
                 minimatrix* H2) const;

    /**
     * Calculate bearing to a landmark
     * @param point 3D location of landmark
     * @return bearing (Unit3)
     */
    Unit3 bearing(const minivector& point) const;
    Unit3 bearing(const minivector& point, minimatrix* H1,
                  minimatrix* H2) const;
    /**
     * Calculate bearing to another pose
     * @param other 3D location and orientation of other body. The orientation
     * information is ignored.
     * @return bearing (Unit3)
     */
    Unit3 bearing(const Pose3& pose) const;

    Unit3 bearing(const Pose3& pose, minimatrix* H1,
                  minimatrix* H2) const;
    /// @}
    /// @name Advanced Interface
    /// @{

    /**
     * Return the start and end indices (inclusive) of the translation component of the
     * exponential map parameterization
     * @return a pair of [start, end] indices into the tangent space vector
     */
    inline static std::pair<int, int> translationInterval();

    /**
     * Return the start and end indices (inclusive) of the rotation component of the
     * exponential map parameterization
     * @return a pair of [start, end] indices into the tangent space vector
     */
    static std::pair<int, int> rotationInterval();

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
inline minimatrix wedge(const minivector& xi)
{
    return Pose3::wedge(xi.data[0], xi.data[1], xi.data[2], xi.data[3], xi.data[4], xi.data[5]);
}
};
#endif // POSE3_H
