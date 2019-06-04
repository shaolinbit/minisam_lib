/*
 * @file Unit3.h
 * @date Feb 02, 2011
 * @author
 * @brief Develop a Unit3 class - basically a point on a unit sphere
 */

#pragma once

#include "../base/Matrix.h"
#include <string>


/// Represents a 3D point on a unit sphere.
class Unit3
{

private:

    Eigen::Vector3d p_; ///< The location of the point on the unit sphere
    Eigen::MatrixXd* B_;
    Eigen::MatrixXd* H_B_;

public:

    enum
    {
        dimension = 2
    };

    /// @name Constructors
    /// @{

    /// Default constructor
    Unit3()
    {
        // p_=Eigen::VectorXd(3);
        p_<<1.0, 0.0, 0.0;
    }

    /// Construct from point
    explicit Unit3(const Eigen::Vector3d& p)
    {
        p_=p.normalized();
    }

    /// Construct from x,y,z
    Unit3(double x, double y, double z)
    {
        // p_=Eigen::VectorXd(3);
        p_<<x,y,z;
        p_.normalize();
    }

    /// Construct from 2D point in plane at focal length f
    /// Unit3(p,1) can be viewed as normalized homogeneous coordinates of 2D point
    explicit Unit3(const Eigen::Vector2d& p, double f)
    {
        //p_=Eigen::VectorXd(3);
        p_<<p(0), p(1), f;
        p_.normalize();
    }

    /// Copy constructor
    Unit3(const Unit3& u)
    {
        p_ = u.p_;
    }

    /// Copy assignment
    Unit3& operator=(const Unit3 & u)
    {
        p_ = u.p_;
        return *this;
    }

    /// Named constructor from Point3 with optional Jacobian
    static Unit3 FromPoint3(const Eigen::Vector3d& point);
    static Unit3 FromPoint3(const Eigen::Vector3d& point, //
                            Eigen::MatrixXd* H);

    /// Random direction, using boost::uniform_on_sphere
    //static Unit3 Random(boost::mt19937 & rng);

    /// @}

    /// @name Testable
    /// @{

    friend std::ostream& operator<<(std::ostream& os, const Unit3& pair);

    /// The print fuction
// void print(const std::string& s = std::string()) const;

    /// The equals function with tolerance
    //bool equals(const Unit3& s, double tol = 1e-9) const {
    //  return equal_with_abs_tol(p_, s.p_, tol);
    //}
    /// @}

    /// @name Other functionality
    /// @{

    /**
     * Returns the local coordinate frame to tangent plane
     * It is a 3*2 matrix [b1 b2] composed of two orthogonal directions
     * tangent to the sphere at the current direction.
     * Provides derivatives of the basis with the two basis vectors stacked up as a 6x1.
     */
    //const Matrix32& basis(OptionalJacobian<6, 2> H = boost::none) const;
    const Eigen::MatrixXd& basis() const;
    const Eigen::MatrixXd& basis(Eigen::MatrixXd* H) const;
    /// Return skew-symmetric associated with 3D point on unit sphere
    Eigen::Matrix3d skew() const;

    /// Return unit-norm Point3
    Eigen::Vector3d point3() const;
    Eigen::Vector3d point3(Eigen::MatrixXd* H) const;
    /// Return unit-norm Vector
    Eigen::Vector3d unitVector() const;
    Eigen::Vector3d unitVector(Eigen::MatrixXd* H) const;
    /// Return scaled direction as Point3
    friend Eigen::Vector3d operator*(double s, const Unit3& d)
    {
        //return Point3(s * d.p_);
        Eigen::Vector3d p3=s*d.p_;
        return p3;
    }

    /// Return dot product with q
    double dot(const Unit3& q) const;
    double dot(const Unit3& q, Eigen::MatrixXd* H1, //
               Eigen::MatrixXd* H2) const;

    /// Signed, vector-valued error between two directions
    /// @deprecated, errorVector has the proper derivatives, this confusingly has only the second.
    Eigen::Vector2d error(const Unit3& q) const;
    Eigen::Vector2d error(const Unit3& q, Eigen::Matrix2d* H_q) const;
    /// Signed, vector-valued error between two directions
    /// NOTE(hayk): This method has zero derivatives if this (p) and q are orthogonal.
    Eigen::Vector2d errorVector(const Unit3& q) const;

    Eigen::Vector2d errorVector(const Unit3& q, Eigen::Matrix2d* H_p, //
                                Eigen::Matrix2d* H_q) const;

    /// Distance between two directions
    double distance(const Unit3& q) const;
    double distance(const Unit3& q, Eigen::MatrixXd* H) const;
    /// Cross-product between two Unit3s
    Unit3 cross(const Unit3& q) const
    {
        return Unit3(p_.cross(q.p_));
    }

    /// Cross-product w Point3
    Eigen::Vector3d cross(const Eigen::Vector3d& q) const
    {
        return point3().cross(q);
    }

    /// @}

    /// @name Manifold
    /// @{

    /// Dimensionality of tangent space = 2 DOF
    inline static int Dim()
    {
        return 2;
    }

    /// Dimensionality of tangent space = 2 DOF
    inline int dim() const
    {
        return 2;
    }

    enum CoordinatesMode
    {
        EXPMAP, ///< Use the exponential map to retract
        RENORM ///< Retract with vector addition and renormalize.
    };

    /// The retract function
    Unit3 retract(const Eigen::Vector2d& v) const;

    /// The local coordinates function
    Eigen::Vector2d localCoordinates(const Unit3& s) const;

    /// @}

};

