#ifndef PROJECTIONFACTOR_H
#define PROJECTIONFACTOR_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file ProjectionFactor.h
 * @brief Basic bearing factor from 2D measurement
 * @author Chris Beall
 * @author Richard Roberts
 * @author Frank Dellaert
 * @author Alex Cunningham
 */
#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../geometry/SimpleCamera.h"
namespace minisam
{
/**
 * Non-linear factor for a constraint derived from a 2D measurement. The calibration is known here.
 * i.e. the main building block for visual SLAM.
 * @addtogroup SLAM
 */
class GenericProjectionFactor: public NoiseModelFactor2
{
protected:

    // Keep a copy of measurement and calibration for I/O
    Eigen::Vector2d measured_;                    ///< 2D measurement
    Cal3_S2* K_;  ///<  pointer to calibration object
    Pose3* body_P_sensor_; ///< The pose of the sensor in the body frame

    // verbosity handling for Cheirality Exceptions
    bool throwCheirality_; ///< If true, rethrows Cheirality exceptions (default: false)
    bool verboseCheirality_; ///< If true, prints text for Cheirality exceptions (default: false)

public:

    /// Default constructor
    GenericProjectionFactor() :
        measured_(0, 0), throwCheirality_(false), verboseCheirality_(false),NoiseModelFactor2(2)
    {
    }

    /**
     * Constructor
     * TODO: Mark argument order standard (keys, measurement, parameters)
     * @param measured is the 2 dimensional location of point in image (the measurement)
     * @param model is the standard deviation
     * @param poseKey is the index of the camera
     * @param pointKey is the index of the landmark
     * @param K  pointer to the constant calibration
     * @param body_P_sensor is the transform from body to sensor frame (default identity)
     */
    GenericProjectionFactor(const Eigen::Vector2d& measured, GaussianNoiseModel* model,
                            int poseKey, int pointKey,  Cal3_S2* K,
                            Pose3* body_P_sensor=NULL) :
        NoiseModelFactor2(model, poseKey, pointKey,2), measured_(measured), K_(K), body_P_sensor_(body_P_sensor),
        throwCheirality_(false), verboseCheirality_(false) {}

    /**
     * Constructor with exception-handling flags
     * TODO: Mark argument order standard (keys, measurement, parameters)
     * @param measured is the 2 dimensional location of point in image (the measurement)
     * @param model is the standard deviation
     * @param poseKey is the index of the camera
     * @param pointKey is the index of the landmark
     * @param K  pointer to the constant calibration
     * @param throwCheirality determines whether Cheirality exceptions are rethrown
     * @param verboseCheirality determines whether exceptions are printed for Cheirality
     * @param body_P_sensor is the transform from body to sensor frame  (default identity)
     */
    GenericProjectionFactor(const  Eigen::Vector2d& measured, GaussianNoiseModel* model,
                            int poseKey, int pointKey,  Cal3_S2* K,
                            bool throwCheirality, bool verboseCheirality,
                            Pose3*  body_P_sensor) :
        NoiseModelFactor2(model, poseKey, pointKey,2), measured_(measured), K_(K), body_P_sensor_(body_P_sensor),
        throwCheirality_(throwCheirality), verboseCheirality_(verboseCheirality) {}

    /** Virtual destructor */
    virtual ~GenericProjectionFactor() {}

    /// Evaluate error h(x)-z
    Eigen::VectorXd evaluateError(const Pose3& pose, const Eigen::Vector3d& point,
                                  Eigen::MatrixXd& H1, Eigen::MatrixXd& H2) const
    {
        try
        {
            Eigen::MatrixXd* fdf=NULL;
            if(body_P_sensor_!=NULL)
            {
                Eigen::MatrixXd H0= (*body_P_sensor_).inverse().AdjointMap();
                PinholeCameraCal3S2 camera(pose*(*body_P_sensor_), *K_);

                Eigen::Vector2d reprojectionError(camera.projectPoint(point, &H1, &H2, fdf) - measured_);
                H1 = H1 * H0;
                return reprojectionError;
            }
            else
            {
                PinholeCameraCal3S2 camera(pose, *K_);
                return camera.projectPoint(point, &H1, &H2, fdf) - measured_;
            }
        }
        catch( exception& e)
        {
            H1 = Eigen::MatrixXd::Zero(2,6);
            H2 = Eigen::MatrixXd::Zero(2,3);
            if (verboseCheirality_)
                std::cout << ": Landmark ";//<< DefaultKeyFormatter(this->key2())
            std::cout<< " moved behind camera "  << std::endl;
            if (throwCheirality_)
                throw e;
        }
        return Eigen::Vector2d::Constant(2.0 * K_->fx());
    }

    /// Evaluate error h(x)-z and optionally derivatives
    Eigen::VectorXd evaluateError(const Pose3& pose, const Eigen::Vector3d& point) const
    {
        try
        {
            if(body_P_sensor_!=NULL)
            {
                Eigen::MatrixXd H0= (*body_P_sensor_).inverse().AdjointMap();
                PinholeCameraCal3S2 camera(pose*(*body_P_sensor_), *K_);

                Eigen::Vector2d reprojectionError(camera.projectPoint(point) - measured_);
                // H1 = H1 * H0;
                return reprojectionError;
            }
            else
            {
                PinholeCameraCal3S2 camera(pose, *K_);
                return camera.projectPoint(point) - measured_;
            }
        }
        catch( exception& e)
        {
            if (verboseCheirality_)
                std::cout << ": Landmark ";//<< DefaultKeyFormatter(this->key2())
            std::cout<< " moved behind camera "  << std::endl;
            if (throwCheirality_)
                throw e;
        }
        return Eigen::Vector2d::Constant(2.0 * K_->fx());
    }


    /** return the measurement */
    const Eigen::Vector2d& measured() const
    {
        return measured_;
    }

    /** return the calibration object */
    inline Cal3_S2* calibration()
    {
        return K_;
    }

    /** return verbosity */
    inline bool verboseCheirality() const
    {
        return verboseCheirality_;
    }

    /** return flag for throwing cheirality exceptions */
    inline bool throwCheirality() const
    {
        return throwCheirality_;
    }

    //nonsense for compling;
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x) const
    {
        Eigen::VectorXd xd;
        xd.setZero(6);
        return xd;
    }
    //nonsense for compling;
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x,std::vector<Eigen::MatrixXd>& H) const
    {
        {
            Eigen::VectorXd xd;
            xd.setZero(6);
            return xd;
        }
    }
    //nonsense for compling;
    virtual Eigen::VectorXd
    evaluateError(const Eigen::VectorXd& X1, const Eigen::VectorXd X2) const
    {
        Eigen::VectorXd xd;
        xd.setZero(6);
        return xd;
    }
//nonsense for compling;
    virtual Eigen::VectorXd
    evaluateError(const Eigen::VectorXd& X1, const Eigen::VectorXd& X2, Eigen::MatrixXd& H1,Eigen::MatrixXd& H2 ) const
    {
        Eigen::VectorXd xd;
        xd.setZero(6);
        return xd;

    }
#ifdef GMF_Using_Pose3
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x1,
                                            const std::map<int,Eigen::VectorXd>& x2,
                                            std::vector<Eigen::MatrixXd>& H) const override
    {
        std::map<int,Pose3>::const_iterator xbegin=x1.find(key1());;
        Pose3 x11=xbegin->second;
        std::map<int,Eigen::VectorXd>::const_iterator  xsecond=x2.find(key2());;
        Eigen::VectorXd x22=xsecond->second;

        // return evaluateError(x1, x2);

        return evaluateError(x11, x22, *(H.begin()), *(H.begin()+1));



    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x1,
                                            const std::map<int,Eigen::VectorXd>& x2) const
    {
        std::map<int,Pose3>::const_iterator xbegin=x1.find(key1());;
        Pose3 x11=xbegin->second;
        std::map<int,Eigen::VectorXd>::const_iterator  xsecond=x2.find(key2());;
        Eigen::VectorXd x22=xsecond->second;

        // return evaluateError(x1, x2);

        return evaluateError(x11, x22);
    }
#endif // GMF_Using_Pose3
    virtual NonlinearFactor* clone()const
    {
        GenericProjectionFactor* newfactor=new GenericProjectionFactor( measured_,noiseModel_,key1(),key2(),K_,
                throwCheirality_, verboseCheirality_,body_P_sensor_);
        return newfactor;
    }

};
};
#endif // PROJECTIONFACTOR_H
