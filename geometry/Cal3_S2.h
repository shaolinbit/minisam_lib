#ifndef CAL3_S2_H
#define CAL3_S2_H

/**
 * @file   Cal3_S2.h
 * @brief  The most common 5DOF 3D->2D calibration
 * @author
 */

/**
 * @addtogroup geometry
 */

#pragma once
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Cholesky>

/**
 * @brief The most common 5DOF 3D->2D calibration
 * @addtogroup geometry
 * \nosubgrouping
 */
class Cal3_S2
{
private:
    double fx_, fy_, s_, u0_, v0_;

public:
    enum { dimension = 5 };
    //typedef boost::shared_ptr<Cal3_S2> shared_ptr; ///< shared pointer to calibration object

    /// @name Standard Constructors
    /// @{

    /// Create a default calibration that leaves coordinates unchanged
    Cal3_S2() :
        fx_(1), fy_(1), s_(0), u0_(0), v0_(0)
    {
    }

    /// constructor from doubles
    Cal3_S2(double fx, double fy, double s, double u0, double v0) :
        fx_(fx), fy_(fy), s_(s), u0_(u0), v0_(v0)
    {
    }

    /// constructor from vector
    Cal3_S2(const Eigen::VectorXd &d) :
        fx_(d(0)), fy_(d(1)), s_(d(2)), u0_(d(3)), v0_(d(4))
    {
    }

    /**
     * Easy constructor, takes fov in degrees, asssumes zero skew, unit aspect
     * @param fov field of view in degrees
     * @param w image width
     * @param h image height
     */
    Cal3_S2(double fov, int w, int h);

    /// @}
    /// @name Advanced Constructors
    /// @{

    /// load calibration from location (default name is calibration_info.txt)
    Cal3_S2(const std::string &path);

    /// @}
    /// @name Testable
    /// @{

    /// Output stream operator
    friend std::ostream &operator<<(std::ostream &os, const Cal3_S2& cal);

    /// print with optional string
// void print(const std::string& s = "Cal3_S2") const;

    /// Check if equal up to specified tolerance
    //bool equals(const Cal3_S2& K, double tol = 10e-9) const;

    /// @}
    /// @name Standard Interface
    /// @{

    /// focal length x
    inline double fx() const
    {
        return fx_;
    }

    /// focal length y
    inline double fy() const
    {
        return fy_;
    }

    /// aspect ratio
    inline double aspectRatio() const
    {
        return fx_/fy_;
    }

    /// skew
    inline double skew() const
    {
        return s_;
    }

    /// image center in x
    inline double px() const
    {
        return u0_;
    }

    /// image center in y
    inline double py() const
    {
        return v0_;
    }

    /// return the principal point
    Eigen::VectorXd principalPoint() const
    {
        Eigen::VectorXd Point2(2);
        Point2<<u0_,v0_;
        return Point2;
        //return Point2(u0_, v0_);
    }

    /// vectorized form (column-wise)
    Eigen::VectorXd vector() const
    {
        Eigen::VectorXd v;
        v << fx_, fy_, s_, u0_, v0_;
        return v;
    }

    /// return calibration matrix K
    Eigen::MatrixXd K() const
    {
        Eigen::MatrixXd K(3,3);
        K <<  fx_, s_, u0_, 0.0, fy_, v0_, 0.0, 0.0, 1.0;
        return K;
    }

    /** @deprecated The following function has been deprecated, use K above */
    Eigen::MatrixXd matrix() const
    {
        return K();
    }

    /// return inverted calibration matrix inv(K)
    Eigen::MatrixXd matrix_inverse() const
    {
        const double fxy = fx_ * fy_, sv0 = s_ * v0_, fyu0 = fy_ * u0_;
        Eigen::MatrixXd K_inverse(3,3);
        K_inverse << 1.0 / fx_, -s_ / fxy, (sv0 - fyu0) / fxy, 0.0,
                  1.0 / fy_, -v0_ / fy_, 0.0, 0.0, 1.0;
        return K_inverse;
    }

    /**
     * convert intrinsic coordinates xy to image coordinates uv, fixed derivaitves
     * @param p point in intrinsic coordinates
     * @param Dcal optional 2*5 Jacobian wrpt Cal3_S2 parameters
     * @param Dp optional 2*2 Jacobian wrpt intrinsic coordinates
     * @return point in image coordinates
     */
    //Point2 uncalibrate(const Point2& p, OptionalJacobian<2,5> Dcal = boost::none,
//      OptionalJacobian<2,2> Dp = boost::none) const;
    Eigen::VectorXd uncalibrate(const Eigen::VectorXd& p);
    Eigen::VectorXd uncalibrate(const Eigen::VectorXd& p, Eigen::MatrixXd* Dcal,
                                Eigen::MatrixXd* Dp) const;
    /**
     * convert image coordinates uv to intrinsic coordinates xy
     * @param p point in image coordinates
     * @param Dcal optional 2*5 Jacobian wrpt Cal3_S2 parameters
     * @param Dp optional 2*2 Jacobian wrpt intrinsic coordinates
     * @return point in intrinsic coordinates
     */
    Eigen::VectorXd calibrate(const Eigen::VectorXd& p) const;
    Eigen::VectorXd calibrate(const Eigen::VectorXd& p, Eigen::MatrixXd* Dcal,
                              Eigen::MatrixXd* Dp) const;

    /**
     * convert homogeneous image coordinates to intrinsic coordinates
     * @param p point in image coordinates
     * @return point in intrinsic coordinates
     */
    Eigen::VectorXd calibrate3(const Eigen::VectorXd& p) const;

    /// "Between", subtracts calibrations. between(p,q) == compose(inverse(p),q)
    inline Cal3_S2 between(const Cal3_S2& q,
                           Eigen::MatrixXd* H1,
                           Eigen::MatrixXd* H2) const
    {
        *H1 =-Eigen::MatrixXd::Identity(5,5);
        *H2 =Eigen::MatrixXd::Identity(5,5);
        return Cal3_S2(q.fx_-fx_, q.fy_-fy_, q.s_-s_, q.u0_-u0_, q.v0_-v0_);
    }

    inline Cal3_S2 between(const Cal3_S2& q) const
    {
        //if(H1) *H1 = -I_5x5;
        //if(H2) *H2 =  I_5x5;
        return Cal3_S2(q.fx_-fx_, q.fy_-fy_, q.s_-s_, q.u0_-u0_, q.v0_-v0_);
    }


    /// @}
    /// @name Manifold
    /// @{

    /// return DOF, dimensionality of tangent space
    inline int dim() const
    {
        return 5;
    }

    /// return DOF, dimensionality of tangent space
    static int Dim()
    {
        return 5;
    }

    /// Given 5-dim tangent vector, create new calibration
    inline Cal3_S2 retract(const Eigen::VectorXd& d) const
    {
        return Cal3_S2(fx_ + d(0), fy_ + d(1), s_ + d(2), u0_ + d(3), v0_ + d(4));
    }

    /// Unretraction for the calibration
    Eigen::VectorXd localCoordinates(const Cal3_S2& T2) const
    {
        return T2.vector() - vector();
    }

};

#endif // CAL3_S2_h
