#ifndef PROJECTIONFACTOR_H
#define PROJECTIONFACTOR_H

/**
 * @file ProjectionFactor.h
 * @brief Basic bearing factor from 2D measurement
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
public:

    // Keep a copy of measurement and calibration for I/O
    minivector measured_;                    ///< 2D measurement
    Cal3_S2* K_;  ///<  pointer to calibration object
    Pose3* body_P_sensor_; ///< The pose of the sensor in the body frame

public:

    /// Default constructor
    GenericProjectionFactor() :
        measured_(minivector(2,0.0)), NoiseModelFactor2(),K_(NULL),body_P_sensor_(NULL)//,throwCheirality_(false), verboseCheirality_(false),
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
    GenericProjectionFactor(const minivector& measured, GaussianNoiseModel* model,
                            int poseKey, int pointKey,  Cal3_S2* K,
                            Pose3* body_P_sensor=NULL) :
        NoiseModelFactor2(model, poseKey, pointKey), measured_(measured), K_(K), body_P_sensor_(body_P_sensor)//,throwCheirality_(false), verboseCheirality_(false)
        {

        }


    /** Virtual destructor */
    virtual ~GenericProjectionFactor() {}

    /// Evaluate error h(x)-z
    minivector evaluateError(const minimatrix* pose, const minimatrix* point,
                             minimatrix& H1, minimatrix& H2) const
    {
        try
        {
            if(body_P_sensor_!=NULL)
            {
                minimatrix H0= body_P_sensor_->inverse().AdjointMap();
                Pose3 ppose(pose);
                PinholeCameraCal3S2 camera(ppose.multiply(*body_P_sensor_), *K_);

                minivector reprojectionError=camera.projectPoint(point, &H1, &H2, NULL);
                minivector_sub(&reprojectionError,measured_);

                //H1 = H1 * H0;;
                minimatrix temph(H1.size1,H1.size2);
                minimatrix_memcpy(&temph,H1);
                miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,temph,H0,0.0,&H1);


                return reprojectionError;
            }
            else
            {
                PinholeCameraCal3S2 camera(pose, *K_);
                //return camera.projectPoint(point, &H1, &H2, fdf) - measured_;
                minivector reprojectionError=camera.projectPoint(point, &H1, &H2, NULL);
                minivector_sub(&reprojectionError,measured_);
                return reprojectionError;
            }
        }
        catch( exception& e)
        {
            minimatrix_resize(&H1,2,6);
            minimatrix_set_zero(&H1);
            minimatrix_resize(&H2,2,3);
            minimatrix_set_zero(&H2);
            std::cout << ": Landmark "<< this->key2()<< "  moved behind camera "  << std::endl;
            throw e;
        }
        return minivector_2dim(2.0*K_->fx(),2.0*K_->fx());
    }

    /// Evaluate error h(x)-z and optionally derivatives
    minivector evaluateError(const minimatrix* pose, const minimatrix* point) const
    {
        try
        {

            Pose3 ppose(pose);
            minivector ppoint(point);
            if(body_P_sensor_!=NULL)
            {
                PinholeCameraCal3S2 camera(ppose.multiply(*body_P_sensor_), *K_);
                minivector reprojectionError=camera.projectPoint(ppoint,NULL,NULL,NULL);
                minivector_sub(&reprojectionError,measured_);
                // H1 = H1 * H0;
                return reprojectionError;
            }
            else
            {
                PinholeCameraCal3S2 camera(ppose, *K_);
                minivector reprojectionError=camera.projectPoint(ppoint,NULL,NULL,NULL);
                minivector_sub(&reprojectionError,measured_);
                return reprojectionError;
            }
        }
        catch( exception& e)
        {
                std::cout << ": Landmark ";
            std::cout<< " moved behind camera "  << std::endl;
                throw e;
        }
        return minivector_2dim(2.0*K_->fx(),2.0*K_->fx());
    }


    /** return the measurement */
    const minivector& measured() const
    {
        return measured_;
    }

    /** return the calibration object */
    inline Cal3_S2* calibration()
    {
        return K_;
    }

     virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x,
                                       std::vector<minimatrix> &H) const
    {
        std::map<int,minimatrix*>::const_iterator itb1=x.find(key1());
        std::map<int,minimatrix*>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second,
                             *(H.begin()),*(H.begin()+1));
    }

    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x) const
    {
        std::map<int,minimatrix*>::const_iterator itb1=x.find(key1());
        std::map<int,minimatrix*>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second);
    }


    virtual NoiseModelFactor* clone()const
    {
        GenericProjectionFactor* newfactor=new GenericProjectionFactor( measured_,noiseModel_,key1(),key2(),K_,
               body_P_sensor_ );
        return newfactor;
    }

};
};
#endif // PROJECTIONFACTOR_H
