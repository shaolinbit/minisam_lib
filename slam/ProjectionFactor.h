#ifndef PROJECTIONFACTOR_H
#define PROJECTIONFACTOR_H


/**
 * @file ProjectionFactor.h
 * @brief Basic bearing factor from 2D measurement
 * @author
 */

#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../geometry/SimpleCamera.h"
//#include <boost/optional.hpp>


/**
 * Non-linear factor for a constraint derived from a 2D measurement. The calibration is known here.
 * i.e. the main building block for visual SLAM.
 * @addtogroup SLAM
 */
// template<class POSE, class LANDMARK, class CALIBRATION = Cal3_S2>
class GenericProjectionFactor: public NoiseModelFactor2
{
protected:

    // Keep a copy of measurement and calibration for I/O
    Eigen::Vector2d measured_;                    ///< 2D measurement
    Cal3_S2* K_;  ///< shared pointer to calibration object
    Pose3* body_P_sensor_; ///< The pose of the sensor in the body frame

    // verbosity handling for Cheirality Exceptions
    bool throwCheirality_; ///< If true, rethrows Cheirality exceptions (default: false)
    bool verboseCheirality_; ///< If true, prints text for Cheirality exceptions (default: false)

public:

    /// shorthand for base class type
    // typedef NoiseModelFactor2<POSE, LANDMARK> Base;

    /// shorthand for this class
    //typedef GenericProjectionFactor<POSE, LANDMARK, CALIBRATION> This;

    /// shorthand for a smart pointer to a factor
    // typedef boost::shared_ptr<This> shared_ptr;

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
     * @param K shared pointer to the constant calibration
     * @param body_P_sensor is the transform from body to sensor frame (default identity)
     */
    GenericProjectionFactor(const Eigen::Vector2d& measured, SharedNoiseModel* model,
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
     * @param K shared pointer to the constant calibration
     * @param throwCheirality determines whether Cheirality exceptions are rethrown
     * @param verboseCheirality determines whether exceptions are printed for Cheirality
     * @param body_P_sensor is the transform from body to sensor frame  (default identity)
     */
    GenericProjectionFactor(const  Eigen::Vector2d& measured, SharedNoiseModel* model,
                            int poseKey, int pointKey,  Cal3_S2* K,
                            bool throwCheirality, bool verboseCheirality,
                            Pose3*  body_P_sensor) :
        NoiseModelFactor2(model, poseKey, pointKey,2), measured_(measured), K_(K), body_P_sensor_(body_P_sensor),
        throwCheirality_(throwCheirality), verboseCheirality_(verboseCheirality) {}

    /** Virtual destructor */
    virtual ~GenericProjectionFactor() {}

    /// @return a deep copy of this factor
    //virtual NonlinearFactor* clone() const {
    //  return NonlinearFactor*( new GenericProjectionFactor()); }

    /**
     * print
     * @param s optional string naming the factor
     * @param keyFormatter optional formatter useful for printing Symbols
     */
    /// Evaluate error h(x)-z and optionally derivatives
    Eigen::VectorXd evaluateError(const Pose3& pose, const Eigen::Vector3d& point,
                                  Eigen::MatrixXd& H1, Eigen::MatrixXd& H2) const
    {
        try
        {
            Eigen::MatrixXd* fdf=NULL;
            if(body_P_sensor_!=NULL)
            {
                // if(H1!=NULL) {
                Eigen::MatrixXd H0= (*body_P_sensor_).inverse().AdjointMap();
                PinholeCameraCal3S2 camera(pose*(*body_P_sensor_), *K_);

                //Eigen::MatrixXd* fdf=NULL;
                //PinholeCameraCal3S2 camera(pose.compose(*body_P_sensor_, H0), *K_);
                Eigen::Vector2d reprojectionError(camera.projectPoint(point, &H1, &H2, fdf) - measured_);
                H1 = H1 * H0;
                return reprojectionError;
                // } else {
                //    PinholeCameraCal3S2 camera(pose*(*body_P_sensor_), *K_);
                //    return camera.projectPoint(point, H1, H2, NULL) - measured_;
                //  }
            }
            else
            {
                PinholeCameraCal3S2 camera(pose, *K_);
                return camera.projectPoint(point, &H1, &H2, fdf) - measured_;
            }
        }
        catch( exception& e)
        {
            //if (H1!=NULL)
            H1 = Eigen::MatrixXd::Zero(2,6);
            // if (H2!=NULL)
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
            // Eigen::MatrixXd* fdf=NULL;
            if(body_P_sensor_!=NULL)
            {
                // if(H1!=NULL) {
                Eigen::MatrixXd H0= (*body_P_sensor_).inverse().AdjointMap();
                PinholeCameraCal3S2 camera(pose*(*body_P_sensor_), *K_);

                //Eigen::MatrixXd* fdf=NULL;
                //PinholeCameraCal3S2 camera(pose.compose(*body_P_sensor_, H0), *K_);
                Eigen::Vector2d reprojectionError(camera.projectPoint(point) - measured_);
                // H1 = H1 * H0;
                return reprojectionError;
                // } else {
                //    PinholeCameraCal3S2 camera(pose*(*body_P_sensor_), *K_);
                //    return camera.projectPoint(point, H1, H2, NULL) - measured_;
                //  }
            }
            else
            {
                PinholeCameraCal3S2 camera(pose, *K_);
                return camera.projectPoint(point) - measured_;
            }
        }
        catch( exception& e)
        {
            //if (H1!=NULL)
            //  H1 = Eigen::MatrixXd::Zero(2,6);
            // if (H2!=NULL)
            // H2 = Eigen::MatrixXd::Zero(2,3);
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
#endif // PROJECTIONFACTOR_H
