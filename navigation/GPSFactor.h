#ifndef GPSFACTOR_H
#define GPSFACTOR_H
/**
 *  @file   GPSFactor.h
 *  @brief  Header file for GPS factor
 **/
#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../navigation/NavState.h"
#include "../geometry/Pose3.h"

namespace minisam
{
/**
 * Prior on position in a Cartesian frame.
 * Possibilities include:
 *   ENU: East-North-Up navigation frame at some local origin
 *   NED: North-East-Down navigation frame at some local origin
 *   ECEF: Earth-centered Earth-fixed, origin at Earth's center
 * See Farrell08book or e.g. http://www.dirsig.org/docs/new/coordinates.html
 * @addtogroup Navigation
 */
class  GPSFactor: public NoiseModelFactor1
{

private:


    minimatrix* nT_; ///< Position measurement in cartesian coordinates

public:


    /** default constructor - only use for serialization */
    GPSFactor()
    {
        nT_=new minivector(3);
        minimatrix_set_zero(nT_);
    }

    virtual ~GPSFactor()
    {
        if(nT_!=NULL)
            delete nT_;
    }

    /**
     * @brief Constructor from a measurement in a Cartesian frame.
     * Use GeographicLib to convert from geographic (latitude and longitude) coordinates
     * @param key of the Pose3 variable that will be constrained
     * @param gpsIn measurement already in correct coordinates
     * @param model Gaussian noise model
     */
    GPSFactor(int key, const minimatrix& gpsIn,  GaussianNoiseModel* model) :
        NoiseModelFactor1(model, key)
    {
        nT_=new minivector(3);
        minimatrix_memcpy(nT_,gpsIn);
    }

    /// @return a deep copy of this factor
    virtual NoiseModelFactor* clone()
    {
        GPSFactor* ngs=new GPSFactor(this->keys()[0],*(this->nT_),this->noiseModel_);
        return ngs;
    }

    /// vector of errors
    virtual minivector evaluateError(const minimatrix* x) const;
    virtual minivector evaluateError(const minimatrix* x, minimatrix &H) const;

    inline const minivector & measurementIn() const
    {
        return minivector(*nT_);
    }


    /**
     *  Convenience function to estimate state at time t, given two GPS
     *  readings (in local NED Cartesian frame) bracketing t
     *  Assumes roll is zero, calculates yaw and pitch from NED1->NED2 vector.
     */
    static std::pair<Pose3, minivector> EstimateState(double t1, const minivector& NED1,
            double t2, const minivector& NED2, double timestamp);
};

/**
 * Version of GPSFactor for NavState
 * @addtogroup Navigation
 */
class  GPSFactor2: public NoiseModelFactor1
{

private:
    minimatrix* nT_; ///< Position measurement in cartesian coordinates
public:

    /// default constructor - only use for serialization
    GPSFactor2()
    {
        nT_=new minivector(3);
        minimatrix_set_zero(nT_);
    }

    virtual ~GPSFactor2()
    {

    }

    /// Constructor from a measurement in a Cartesian frame.
    GPSFactor2(int key, const minimatrix& gpsIn, GaussianNoiseModel*  model) :
        NoiseModelFactor1(model, key)//, nT_(gpsIn)
    {
        nT_=new minivector(3);
        minimatrix_memcpy(nT_,gpsIn);
    }

    /// @return a deep copy of this factor
    virtual NoiseModelFactor* clone()
    {
        GPSFactor2* ngs=new GPSFactor2(this->keys()[0],*(this->nT_),this->noiseModel_);
        return ngs;
    }


    virtual minivector evaluateError(const minimatrix* x) const;

    virtual minivector evaluateError(const minimatrix* x, minimatrix &H) const;

    inline const minivector & measurementIn() const
    {
        return minivector(*nT_);
    }
};

} /// namespace minisam
#endif
