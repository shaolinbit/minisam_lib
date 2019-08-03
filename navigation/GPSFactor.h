/**
 *  @file   GPSFactor.h
 *  @author
 *  @brief  Header file for GPS factor
 *  @date
 **/
#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../navigation/NavState.h"
#include "../geometry/Pose3.h"

namespace minisam {
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


    Eigen::Vector3d nT_; ///< Position measurement in cartesian coordinates

public:


    /** default constructor - only use for serialization */
    GPSFactor(): nT_(0, 0, 0) {}

    virtual ~GPSFactor() {}

    /**
     * @brief Constructor from a measurement in a Cartesian frame.
     * Use GeographicLib to convert from geographic (latitude and longitude) coordinates
     * @param key of the Pose3 variable that will be constrained
     * @param gpsIn measurement already in correct coordinates
     * @param model Gaussian noise model
     */
    GPSFactor(int key, const Eigen::Vector3d& gpsIn,  SharedNoiseModel* model) :
        NoiseModelFactor1(model, key,1), nT_(gpsIn)
    {
    }

    /// @return a deep copy of this factor
    virtual NonlinearFactor* clone()
    {
        GPSFactor* ngs=new GPSFactor(this->keys()[0],this->nT_,this->noiseModel_);
        return ngs;
    }

    /// vector of errors
    virtual Eigen::VectorXd evaluateError(const Pose3& p) const;

    virtual Eigen::VectorXd evaluateError(const Pose3& p,
                                  Eigen::MatrixXd& H) const;

    inline const Eigen::Vector3d & measurementIn() const
    {
        return nT_;
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x,std::vector<Eigen::MatrixXd>& H) const
    {
        std::map<int,Pose3>::const_iterator itb=x.find(key());
        return evaluateError(itb->second,*(H.begin()));
    }

    //nonsense for virtual;
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Eigen::VectorXd>& x)const
    {
        Eigen::VectorXd uw(6);
        uw.setZero();
        return uw;
    }


    /**
     *  Convenience function to estimate state at time t, given two GPS
     *  readings (in local NED Cartesian frame) bracketing t
     *  Assumes roll is zero, calculates yaw and pitch from NED1->NED2 vector.
     */
    static std::pair<Pose3, Eigen::Vector3d> EstimateState(double t1, const Eigen::Vector3d& NED1,
            double t2, const Eigen::Vector3d& NED2, double timestamp);
};

/**
 * Version of GPSFactor for NavState
 * @addtogroup Navigation
 */
class  GPSFactor2: public NoiseModelFactor1
{

private:
    Eigen::Vector3d nT_; ///< Position measurement in cartesian coordinates
public:

    /// default constructor - only use for serialization
    GPSFactor2():nT_(Eigen::Vector3d::Zero()) {}

    virtual ~GPSFactor2() {}

    /// Constructor from a measurement in a Cartesian frame.
    GPSFactor2(int key, const Eigen::Vector3d& gpsIn, SharedNoiseModel*  model) :
        NoiseModelFactor1(model, key), nT_(gpsIn)
    {
    }

    /// @return a deep copy of this factor
    virtual NonlinearFactor* clone()
    {
        GPSFactor2* ngs=new GPSFactor2(this->keys()[0],this->nT_,this->noiseModel_);
        return ngs;
    }


    Eigen::VectorXd evaluateError(const NavState& p,
                                  Eigen::MatrixXd* H) const;

    inline const Eigen::Vector3d & measurementIn() const
    {
        return nT_;
    }
};

} /// namespace minisam
