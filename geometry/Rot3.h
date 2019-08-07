#ifndef ROT3_H
#define ROT3_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    Rot3.h
 * @brief   3D rotation represented as a rotation matrix or quaternion
 * @author  Alireza Fathi
 * @author  Christian Potthast
 * @author  Frank Dellaert
 * @author  Richard Roberts
 * @author  Luca Carlone
 */

#pragma once

#include "../geometry/Unit3.h"
#include "../geometry/Quaternion.h"
#include "../geometry/SO3.h"
#include "../gmfconfig.h"

// Get USE_QUATERNIONS macro

// You can override the default coordinate mode using this flag
#ifndef ROT3_DEFAULT_COORDINATES_MODE
#ifdef USE_QUATERNIONS
// Exponential map is very cheap for quaternions
#define ROT3_DEFAULT_COORDINATES_MODE Rot3::EXPMAP
#else
// If user doesn't require ROT3_EXPMAP in cmake when building
#ifndef ROT3_EXPMAP
// For rotation matrices, the Cayley transform is a fast retract alternative
#define ROT3_DEFAULT_COORDINATES_MODE Rot3::CAYLEY
#else
#define ROT3_DEFAULT_COORDINATES_MODE Rot3::EXPMAP
#endif
#endif
#endif
namespace minisam
{


/**
 * @brief A 3D rotation represented as a rotation matrix if the preprocessor
 * symbol USE_QUATERNIONS is not defined, or as a quaternion if it
 * is defined.
 * @addtogroup geometry
 * \nosubgrouping
 */
class Rot3
{
public:

#ifdef USE_QUATERNIONS
    /** Internal Eigen Quaternion */
    QQuaternion quaternion_;
#else
    Eigen::Matrix3d rot_;
#endif

public:

    /// @name Constructors and named constructors
    /// @{

    /** default constructor, unit rotation */
    Rot3();

    /**
     * Constructor from *columns*
     * @param r1 X-axis of rotated frame
     * @param r2 Y-axis of rotated frame
     * @param r3 Z-axis of rotated frame
     */
    Rot3(const Eigen::Vector3d& col1, const Eigen::Vector3d& col2, const Eigen::Vector3d& col3);

    /** constructor from a rotation matrix, as doubles in *row-major* order !!! */
    Rot3(double R11, double R12, double R13,
         double R21, double R22, double R23,
         double R31, double R32, double R33);

    /**
     * Constructor from a rotation matrix
     * Version for generic matrices. Need casting to Matrix3
     * in quaternion mode, since Eigen's quaternion doesn't
     * allow assignment/construction from a generic matrix.
     * See: http://stackoverflow.com/questions/27094132/cannot-understand-if-this-is-circular-dependency-or-clang#tab-top
     */

    inline explicit Rot3(const Eigen::Matrix3d& R)
    {
#ifdef USE_QUATERNIONS
        quaternion_=QQuaternion(R);
#else
        rot_ = R;
#endif
    }


    /** Constructor from a quaternion.  This can also be called using a plain
     * Vector, due to implicit conversion from Vector to Quaternion
     * @param q The quaternion
     */
    Rot3(const QQuaternion& q);
    Rot3(double x, double y, double z, double w) : Rot3(QQuaternion(x, y, z, w)) {}

    /** Virtual destructor */
    virtual ~Rot3() {}

    /* Static member function to generate some well known rotations */

    /// Rotation around X axis as in http://en.wikipedia.org/wiki/Rotation_matrix, counterclockwise when looking from unchanging axis.
    static Rot3 Rx(double t);

    /// Rotation around Y axis as in http://en.wikipedia.org/wiki/Rotation_matrix, counterclockwise when looking from unchanging axis.
    static Rot3 Ry(double t);

    /// Rotation around Z axis as in http://en.wikipedia.org/wiki/Rotation_matrix, counterclockwise when looking from unchanging axis.
    static Rot3 Rz(double t);

    /// Rotations around Z, Y, then X axes as in http://en.wikipedia.org/wiki/Rotation_matrix, counterclockwise when looking from unchanging axis.
    static Rot3 RzRyRx(double x, double y, double z);

    /// Rotations around Z, Y, then X axes as in http://en.wikipedia.org/wiki/Rotation_matrix, counterclockwise when looking from unchanging axis.
    inline static Rot3 RzRyRx(const Eigen::VectorXd& xyz)
    {
        assert(xyz.size() == 3);
        return RzRyRx(xyz(0), xyz(1), xyz(2));
    }

    /// Positive yaw is to right (as in aircraft heading). See ypr
    static Rot3 Yaw  (double t);
    /// Positive pitch is up (increasing aircraft altitude).See ypr
    static Rot3 Pitch(double t);

    //// Positive roll is to right (increasing yaw in aircraft).
    static Rot3 Roll (double t);

    /**
     * Returns rotation nRb from body to nav frame.
     * For vehicle coordinate frame X forward, Y right, Z down:
     * Positive yaw is to right (as in aircraft heading).
     * Positive pitch is up (increasing aircraft altitude).
     * Positive roll is to right (increasing yaw in aircraft).
     * Tait-Bryan system from Spatial Reference Model (SRM) (x,y,z) = (roll,pitch,yaw)
     * as described in http://www.sedris.org/wg8home/Documents/WG80462.pdf.
     *
     * For vehicle coordinate frame X forward, Y left, Z up:
     * Positive yaw is to left (as in aircraft heading).
     * Positive pitch is down (decreasing aircraft altitude).
     * Positive roll is to right (decreasing yaw in aircraft).
     */
    static Rot3 Ypr(double y, double p, double r);

    /** Create from Quaternion coefficients */
    static Rot3 Quaternion(double w, double x, double y, double z);
    /**
     * Convert from axis/angle representation
     * @param  axisw is the rotation axis, unit length
     * @param   angle rotation angle
     * @return incremental rotation
     */
    static Rot3 AxisAngle(const Eigen::Vector3d& axis, double angle);

    /**
     * Convert from axis/angle representation
     * @param   axis is the rotation axis
     * @param   angle rotation angle
     * @return incremental rotation
     */
    static Rot3 AxisAngle(const Unit3& axis, double angle);
    /**
     * Rodrigues' formula to compute an incremental rotation
     * @param w a vector of incremental roll,pitch,yaw
     * @return incremental rotation
     */
    static Rot3 Rodrigues(const Eigen::Vector3d& w);
    /**
     * Rodrigues' formula to compute an incremental rotation
     * @param wx Incremental roll (about X)
     * @param wy Incremental pitch (about Y)
     * @param wz Incremental yaw (about Z)
     * @return incremental rotation
     */
    static Rot3 Rodrigues(double wx, double wy, double wz);

    /// Determine a rotation to bring two vectors into alignment, using the rotation axis provided
    static Rot3 AlignPair(const Unit3& axis, const Unit3& a_p, const Unit3& b_p);

    /// Calculate rotation from two pairs of homogeneous points using two successive rotations
    static Rot3 AlignTwoPairs(const Unit3& a_p, const Unit3& b_p,  //
                              const Unit3& a_q, const Unit3& b_q);

    /// @}
    /// @name Testable
    /// @{

    /** equals with an tolerance */
    bool equals(const Rot3& p, double tol = 1e-9) const;

    /// @}
    /// @name Group
    /// @{

    /// identity rotation for group operation
    inline static Rot3 identity()
    {
        return Rot3();
    }

    /// Syntatic sugar for composing two rotations
    Rot3 operator*(const Rot3& R2) const;

    /// inverse of a rotation, TODO should be different for M/Q
    Rot3 inverse() const;

    /**
     * Conjugation: given a rotation acting in frame B, compute rotation c1Rc2 acting in a frame C
     * @param cRb rotation from B frame to C frame
     * @return c1Rc2 = cRb * b1Rb2 * cRb'
     */
    Rot3 conjugate(const Rot3& cRb) const;
    /// @}
    /// @name Manifold
    /// @{

    /**
     * The method retract() is used to map from the tangent space back to the manifold.
     * Its inverse, is localCoordinates(). For Lie groups, an obvious retraction is the
     * exponential map, but this can be expensive to compute. The following Enum is used
     * to indicate which method should be used.  The default
     * is determined by ROT3_DEFAULT_COORDINATES_MODE, which may be set at compile time,
     * and itself defaults to Rot3::CAYLEY, or if USE_QUATERNIONS is defined,
     * to Rot3::EXPMAP.
     */
    enum CoordinatesMode
    {
        EXPMAP, ///< Use the Lie group exponential map to retract
#ifndef USE_QUATERNIONS
        CAYLEY, ///< Retract and localCoordinates using the Cayley transform.
#endif
    };

#ifndef USE_QUATERNIONS

    // Cayley chart around origin
    struct CayleyChart
    {
        static Rot3 retract(const Eigen::Vector3d& v);
        static Rot3 retract(const Eigen::Vector3d& v, Eigen::Matrix3d* H);
        static Eigen::Vector3d Local(const Rot3& r);
        static Eigen::Vector3d Local(const Rot3& r, Eigen::Matrix3d* H);
    };

    /// Retraction from R^3 to Rot3 manifold using the Cayley transform
    Rot3 retractCayley(const Eigen::VectorXd& omega) const;

    /// Inverse of retractCayley
    Eigen::Vector3d localCayley(const Rot3& other) const;



#endif
    Rot3 between(const Rot3& g, Eigen::Matrix3d* H1=NULL,
                 Eigen::Matrix3d* H2 =NULL) const;

    /// @}
    /// @name Lie Group
    /// @{

    /**
     * Exponential map at identity - create a rotation from canonical coordinates
     * \f$ [R_x,R_y,R_z] \f$ using Rodrigues' formula
     */
    static Rot3 Expmap(const Eigen::Vector3d& v);
    static Rot3 Expmap(const Eigen::Vector3d& v, Eigen::Matrix3d* H);


    /**
     * Log map at identity - returns the canonical coordinates
     * \f$ [R_x,R_y,R_z] \f$ of this rotation
     */
    static Eigen::Vector3d Logmap(const Rot3& R);

    static Eigen::Vector3d Logmap(const Rot3& R, Eigen::Matrix3d* H);

    /// Derivative of Expmap
    static Eigen::Matrix3d ExpmapDerivative(const Eigen::Vector3d& x);

    /// Derivative of Logmap
    static Eigen::Matrix3d LogmapDerivative(const Eigen::Vector3d& x);

    /** Calculate Adjoint map */
    Eigen::Matrix3d AdjointMap() const
    {
        return matrix();
    }

    // Chart at origin, depends on compile-time flag ROT3_DEFAULT_COORDINATES_MODE
    struct ChartAtOrigin
    {
        static Rot3 retract(const Eigen::Vector3d& v);
        static Rot3 retract(const Eigen::Vector3d& v, Eigen::Matrix3d* H);
        static Eigen::Vector3d Local(const Rot3& r);
        static Eigen::Vector3d Local(const Rot3& r, Eigen::Matrix3d* H);
    };

    Rot3 compose(const Rot3& g, Eigen::Matrix3d* H1=NULL,Eigen::Matrix3d* H2=NULL)
    {
        if(H1!=NULL)
            *H1=g.inverse().AdjointMap();
        if(H2!=NULL)
            *H2= Eigen::MatrixXd::Identity(3,3);
        return *this*g;
    }
    Rot3 expmap(const Eigen::VectorXd& v,Eigen::Matrix3d* H1=NULL,Eigen::Matrix3d* H2=NULL)
    {
        Eigen::Matrix3d D_g_v;
        Rot3 g;
        if(H2!=NULL)
        {
            g=this->Expmap(v,&D_g_v);
        }
        else
        {
            g=this->Expmap(v);
        }
        Rot3 h=this->compose(g,NULL);

        if(H1!=NULL)
            *H1=g.inverse().AdjointMap();
        if(H2!=NULL)
            *H2=D_g_v;
        return h;
    }
    Rot3 retract(const Eigen::VectorXd& v)
    {
        return Rot3::ChartAtOrigin::retract(v);
    }
    Eigen::Vector3d LocalCoordinates(const Rot3& g)
    {
        return Rot3::ChartAtOrigin::Local(g);
    }

    /// @}
    /// @name Group Action on Point3
    /// @{

    /**
     * rotate point from rotated coordinate frame to world \f$ p^w = R_c^w p^c \f$
     */
    Eigen::Vector3d rotatePoint(const Eigen::Vector3d& p, Eigen::Matrix3d* H1=NULL,
                                Eigen::Matrix3d* H2=NULL) const;

    /// rotate point from rotated coordinate frame to world = R*p
    Eigen::Vector3d operator*(const Eigen::Vector3d& p) const;

    /// rotate point from world to rotated frame \f$ p^c = (R_c^w)^T p^w \f$
    //Eigen::Vector3d unrotatePoint(const Eigen::Vector3d& p) const;
    Eigen::Vector3d unrotatePoint(const Eigen::Vector3d& p, Eigen::Matrix3d* H1=NULL,
                                  Eigen::Matrix3d* H2=NULL) const;

    /// @}
    /// @name Group Action on Unit3
    /// @{

    /// rotate 3D direction from rotated coordinate frame to world frame
    Unit3 rotateUnit(const Unit3& p) const;
    Unit3 rotateUnit(const Unit3& p, Eigen::MatrixXd* HR,
                     Eigen::Matrix3d* Hp ) const;

    /// unrotate 3D direction from world frame to rotated coordinate frame
    Unit3 unrotateUnit(const Unit3& p) const;

    Unit3 unrotateUnit(const Unit3& p, Eigen::MatrixXd* HR,
                       Eigen::Matrix2d* Hp ) const;

    /// rotate 3D direction from rotated coordinate frame to world frame
    Unit3 operator*(const Unit3& p) const;

    /// @}
    /// @name Standard Interface
    /// @{

    /** return 3*3 rotation matrix */
    Eigen::Matrix3d matrix() const;

    /**
     * Return 3*3 transpose (inverse) rotation matrix
     */
    Eigen::Matrix3d transpose() const;

    /// @deprecated, this is base 1, and was just confusing
    Eigen::Vector3d  column(int index) const;

    Eigen::Vector3d  r1() const; ///< first column
    Eigen::Vector3d  r2() const; ///< second column
    Eigen::Vector3d  r3() const; ///< third column

    /**
     * Use RQ to calculate xyz angle representation
     * @return a vector containing x,y,z s.t. R = Rot3::RzRyRx(x,y,z)
     */
    Eigen::Vector3d  xyz() const;

    /**
     * Use RQ to calculate yaw-pitch-roll angle representation
     * @return a vector containing ypr s.t. R = Rot3::Ypr(y,p,r)
     */
    Eigen::Vector3d  ypr() const;

    /**
     * Use RQ to calculate roll-pitch-yaw angle representation
     * @return a vector containing ypr s.t. R = Rot3::Ypr(y,p,r)
     */
    Eigen::Vector3d rpy() const;

    /**
     * Accessor to get to component of angle representations
     * NOTE: these are not efficient to get to multiple separate parts,
     * you should instead use xyz() or ypr()
     * TODO: make this more efficient
     */
    inline double roll() const
    {
        return ypr()(2);
    }

    /**
     * Accessor to get to component of angle representations
     * NOTE: these are not efficient to get to multiple separate parts,
     * you should instead use xyz() or ypr()
     * TODO: make this more efficient
     */
    inline double pitch() const
    {
        return ypr()(1);
    }

    /**
     * Accessor to get to component of angle representations
     * NOTE: these are not efficient to get to multiple separate parts,
     * you should instead use xyz() or ypr()
     * TODO: make this more efficient
     */
    inline double yaw() const
    {
        return ypr()(0);
    }

    /// @}
    /// @name Advanced Interface
    /// @{

    /** Compute the quaternion representation of this rotation.
     * @return The quaternion
     */
    QQuaternion toQuaternion() const;

    /**
     * Converts to a generic matrix to allow for use with matlab
     * In format: w x y z
     */
    Eigen::VectorXd quaternion() const;

    /**
     * @brief Spherical Linear intERPolation between *this and other
     * @param s a value between 0 and 1
     * @param other final point of iterpolation geodesic on manifold
     */
    Rot3 slerp(double t, const Rot3& other) const;

    inline Rot3& operator=(const Rot3& rObj)
    {

#ifdef USE_QUATERNIONS
        quaternion_=rObj.quaternion_;
#else
        rot_=rObj.matrix();
#endif // USE_QUATERNIONS
        return *this;
    }

    /// Output stream operator
    friend std::ostream &operator<<(std::ostream &os, const Rot3& p);

    /// @}


};

/**
 * [RQ] receives a 3 by 3 matrix and returns an upper triangular matrix R
 * and 3 rotation angles corresponding to the rotation matrix Q=Qz'*Qy'*Qx'
 * such that A = R*Q = R*Qz'*Qy'*Qx'. When A is a rotation matrix, R will
 * be the identity and Q is a yaw-pitch-roll decomposition of A.
 * The implementation uses Givens rotations and is based on Hartley-Zisserman.
 * @param A 3 by 3 matrix A=RQ
 * @return an upper triangular matrix R
 * @return a vector [thetax, thetay, thetaz] in radians.
 */
std::pair<Eigen::Matrix3d,Eigen::Vector3d> RQ(const Eigen::Matrix3d& A);
};
#endif // ROT3_H
