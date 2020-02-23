/**
 * @file Cal3Bundler.h
 * @brief Calibration used by Bundler
 */

#pragma once

#include "../mat/MatCal.h"
#include "../mat/Matrix.h"


namespace minisam
{

/**
 * @brief Calibration used by Bundler
 * @addtogroup geometry
 * \nosubgrouping
 */
class  Cal3Bundler:public minivector
{

public:

    /// @name Standard Constructors
    /// @{

    /// Default constructor
    Cal3Bundler();

    /**
     *  Constructor
     *  @param f focal length
     *  @param k1 first radial distortion coefficient (quadratic)
     *  @param k2 second radial distortion coefficient (quartic)
     *  @param u0 optional image center (default 0), considered a constant
     *  @param v0 optional image center (default 0), considered a constant
     */
    Cal3Bundler(double f, double k1, double k2, double u0 = 0, double v0 = 0);

    virtual ~Cal3Bundler() {}

    /// @}

    /// @name Standard Interface
    /// @{

    minimatrix K() const; ///< Standard 3*3 calibration matrix
    minivector k() const; ///< Radial distortion parameters (4 of them, 2 0)

    minivector vector() const;

    /// focal length x
    inline double fx() const
    {
        return data[0];
    }

    /// focal length y
    inline double fy() const
    {
        return data[0];
    }

    /// distorsion parameter k1
    inline double k1() const
    {
        return data[1];
    }

    /// distorsion parameter k2
    inline double k2() const
    {
        return data[2];
    }

    /// get parameter u0
    inline double u0() const
    {
        return data[3];
    }

    /// get parameter v0
    inline double v0() const
    {
        return data[4];
    }


    /**
     * @brief: convert intrinsic coordinates xy to image coordinates uv
     * Version of uncalibrate with derivatives
     * @param p point in intrinsic coordinates
     * @param Dcal optional 2*3 Jacobian wrpt CalBundler parameters
     * @param Dp optional 2*2 Jacobian wrpt intrinsic coordinates
     * @return point in image coordinates
     */
    minivector uncalibrate(const minivector& p,minimatrix* Dcal=NULL,
                           minimatrix* Dp=NULL
                          ) const;

    /// Conver a pixel coordinate to ideal coordinate
    minivector calibrate(const minivector& pi, const double tol = 1e-5) const;

    /// @deprecated might be removed in next release, use uncalibrate
    minimatrix D2d_intrinsic(const minivector& p) const;

    /// @deprecated might be removed in next release, use uncalibrate
    minimatrix D2d_calibration(const minivector& p) const;

    /// @deprecated might be removed in next release, use uncalibrate
    minimatrix D2d_intrinsic_calibration(const minivector& p) const;
    /// @}
    /// @name Manifold
    /// @{

    /// Update calibration with tangent space delta
    virtual minimatrix* Retract(const minimatrix* d);
    /// Calculate local coordinates to another calibration
    virtual minimatrix LocalCoordinates(const minimatrix* T2) const;
    /// dimensionality
    virtual size_t dim() const
    {
        return 3;
    }


};

};
