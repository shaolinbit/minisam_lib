#ifndef PREINTEGRATIONPARAMS_H
#define PREINTEGRATIONPARAMS_H


/**
 *  @file  PreintegrationParams.h
 **/

#pragma once

#include "../navigation/PreintegratedRotation.h"
namespace minisam
{

/// Parameters for pre-integration:
/// Usage: Create just a single Params and pass a  pointer to the constructor
struct PreintegrationParams:public minimatrix// PreintegratedRotationParams minimatrix(4,3)
{
    // minimatrix gyroscopeCovariance;  ///< continuous-time "Covariance" of gyroscope measurements
    //minimatrix accelerometerCovariance; ///< continuous-time "Covariance" of accelerometer
    //minimatrix integrationCovariance; ///< continuous-time "Covariance" describing integration uncertainty
    //bool use2ndOrderCoriolis; ///< Whether Eigen::MatrixXd::Identity(3,3)to use second order Coriolis integration
    //  minivector omegaCoriolis;
    // minivector n_gravity; ///< Gravity vector in nav frame
    PreintegrationParams():minimatrix(11,3)
    {
        minimatrix_set_zero(this);
        data[0]=1.0;
        data[4]=1.0;
        data[8]=1.0;
        data[9]=1.0;
        data[13]=1.0;
        data[17]=1.0;
        data[18]=1.0;
        data[22]=1.0;
        data[26]=1.0;
        data[32]=9.8;
    }
    PreintegrationParams(const minimatrix& copym):minimatrix(11,3)
    {
        minimatrix_memcpy(this,copym);
    }
    PreintegrationParams(const minimatrix* m_memory)
    {
        if(m_memory->size1!=11||m_memory->size2!=3)
        {
            throw std::invalid_argument("size not right");
        }
        size1=11;
        size2=3;
        prd=m_memory->prd;
        data=m_memory->data;
        owner=0;
        dimension=33;
    }


    /// The Params constructor insists on getting the navigation frame gravity vector
    /// For convenience, two commonly used conventions are provided by named constructors below
    PreintegrationParams(const minivector& nn_gravity)
        :minimatrix(11,3)//PreintegratedRotationParams(PreintegratedRotationParams()),use2ndOrderCoriolis(false),accelerometerCovariance(minimatrix(3,3)),integrationCovariance(minimatrix(3,3)),n_gravity(minivector(3))
    {
        minimatrix_set_zero(this);
        data[0]=1.0;
        data[4]=1.0;
        data[8]=1.0;
        data[9]=1.0;
        data[13]=1.0;
        data[17]=1.0;
        data[18]=1.0;
        data[22]=1.0;
        data[26]=1.0;
        data[30]=minivector_get(&nn_gravity,0);
        data[31]=minivector_get(&nn_gravity,1);
        data[32]=minivector_get(&nn_gravity,2);
    }

    // Default Params for a Z-down navigation frame, such as NED: gravity points along positive Z-axis
    static PreintegrationParams MakeSharedD(double g = 9.81)
    {
        return PreintegrationParams(minivector(0, 0, g));
    }

    // Default Params for a Z-up navigation frame, such as ENU: gravity points along negative Z-axis
    static PreintegrationParams MakeSharedU(double g = 9.81)
    {
        return PreintegrationParams(minivector(0, 0, -g));
    }
    ~PreintegrationParams()
    {
    }

    void setGyroscopeCovariance(const minimatrix& cov)
    {
        minimatrix gyroscopeCovariance=minimatrix_blockmatrix_var(this,0,0,3,3);
        minimatrix_memcpy(&gyroscopeCovariance,cov);
    }
    void setAccelerometerCovariance(const minimatrix& cov)
    {
        minimatrix accelerometerCovariance = minimatrix_blockmatrix_var(this,3,0,3,3);
        minimatrix_memcpy(&accelerometerCovariance,cov);
    }
    void setIntegrationCovariance(const minimatrix& cov)
    {
        minimatrix integrationCovariance = minimatrix_blockmatrix_var(this,6,0,3,3);
        minimatrix_memcpy(&integrationCovariance,cov);
    }
    void setOmegaCoriolis(const minivector& omega)
    {
        minivector omegaCoriolis=minimatrix_row(this,9);
        minivector_memcpy(&omegaCoriolis,omega);
    }
    void setGravity(const minivector& Gravity)
    {
        minivector mGravity=minimatrix_row(this,10);
        minivector_memcpy(&mGravity,Gravity);
    }

    /*
    void setUse2ndOrderCoriolis(bool flag)
    {
        use2ndOrderCoriolis = flag;
    }
    */
     minimatrix getGyroscopeCovariance()     const
    {
        minimatrix gyroscopeCovariance=minimatrix_blockmatrix(*this,0,0,3,3);
        return gyroscopeCovariance;
    }
     minimatrix getAccelerometerCovariance() const
    {
        minimatrix accelerometerCovariance = minimatrix_blockmatrix(*this,3,0,3,3);
        return accelerometerCovariance;
    }
     minimatrix getIntegrationCovariance()   const
    {
        minimatrix integrationCovariance = minimatrix_blockmatrix(*this,6,0,3,3);
        return integrationCovariance;
    }
     minivector getOmegaCoriolis() const
    {
        minivector omegaCoriolis=minimatrix_row(*this,9);
        return omegaCoriolis;
    }

     minivector getGravity() const
    {
        minivector mGravity=minimatrix_row(*this,10);
        return mGravity;
    }

    /*
    bool getUse2ndOrderCoriolis()     const
    {
       return use2ndOrderCoriolis;
    }


    void setBodyPSensor(const Pose3& pose)
    {
       *body_P_sensor=pose;
    }
    */




    bool operator==(const PreintegrationParams& other)const
    {
        bool equp=minimatrix_equal(*this,other);
        return equp;
    }



};
};
#endif
