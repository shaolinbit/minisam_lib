#ifndef CAL3_S2STEREO_H
#define CAL3_S2STEREO_H


/**
 * @file   Cal3_S2Stereo.h
 * @brief  The most common 5DOF 3D->2D calibration + Stereo baseline
 * @author
 */

#pragma once

#include "../geometry/Cal3_S2.h"
#include <iosfwd>



/**
 * @brief The most common 5DOF 3D->2D calibration, stereo version
 * @addtogroup geometry
 * \nosubgrouping
 */
class Cal3_S2Stereo
{
private:

    Cal3_S2 K_;
    double b_;

public:

    enum { dimension = 6 };
    // typedef boost::shared_ptr<Cal3_S2Stereo> shared_ptr;  ///< shared pointer to stereo calibration object

    /// @name Standard Constructors
    /// @

    /// default calibration leaves coordinates unchanged
    Cal3_S2Stereo() :
        K_(1, 1, 0, 0, 0), b_(1.0)
    {
    }

    /// constructor from doubles
    Cal3_S2Stereo(double fx, double fy, double s, double u0, double v0, double b) :
        K_(fx, fy, s, u0, v0), b_(b)
    {
    }

    /// constructor from vector
    Cal3_S2Stereo(const Eigen::VectorXd &d): K_(d(0), d(1), d(2), d(3), d(4)), b_(d(5)) {}

    /// easy constructor; field-of-view in degrees, assumes zero skew
    Cal3_S2Stereo(double fov, int w, int h, double b) :
        K_(fov, w, h), b_(b)
    {
    }

    /// @}
    /// @name Testable
    /// @{

    // void print(const std::string& s = "") const;

    /// Check if equal up to specified tolerance
    // bool equals(const Cal3_S2Stereo& other, double tol = 10e-9) const;

    /// @}
    /// @name Standard Interface
    /// @{

    /// return calibration, same for left and right
    const Cal3_S2& calibration() const
    {
        return K_;
    }

    /// return calibration matrix K, same for left and right
    Eigen::MatrixXd matrix() const
    {
        return K_.matrix();
    }

    /// focal length x
    inline double fx() const
    {
        return K_.fx();
    }

    /// focal length x
    inline double fy() const
    {
        return K_.fy();
    }

    /// skew
    inline double skew() const
    {
        return K_.skew();
    }

    /// image center in x
    inline double px() const
    {
        return K_.px();
    }

    /// image center in y
    inline double py() const
    {
        return K_.py();
    }

    /// return the principal point
    Eigen::Vector2d principalPoint() const
    {
        return K_.principalPoint();
    }

    /// return baseline
    inline double baseline() const
    {
        return b_;
    }

    /// vectorized form (column-wise)
    Eigen::VectorXd vector() const
    {
        Eigen::VectorXd v(6);
        v << K_.vector(), b_;
        return v;
    }

    /// @}
    /// @name Manifold
    /// @{

    /// return DOF, dimensionality of tangent space
    inline int dim() const
    {
        return 6;
    }

    /// return DOF, dimensionality of tangent space
    static int Dim()
    {
        return 6;
    }

    /// Given 6-dim tangent vector, create new calibration
    inline Cal3_S2Stereo retract(const Eigen::VectorXd& d) const
    {
        return Cal3_S2Stereo(K_.fx() + d(0), K_.fy() + d(1), K_.skew() + d(2), K_.px() + d(3), K_.py() + d(4), b_ + d(5));
    }

    /// Unretraction for the calibration
    Eigen::VectorXd localCoordinates(const Cal3_S2Stereo& T2) const
    {
        return T2.vector() - vector();
    }



};

#endif // CAL3_S2STEREO_H
