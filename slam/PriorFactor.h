#ifndef PRIORFACTOR_H
#define PRIORFACTOR_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  PriorFactor.h
 *  @author Frank Dellaert
 **/
#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../linear/NoiseModel.h"
namespace minisam
{

/**
 * A class for a soft prior on Vector Value type
 * @addtogroup SLAM
 */

class PriorFactor: public NoiseModelFactor1
{

private:
    Eigen::VectorXd prior_; /** The measurement */
public:
    /** default constructor - only use for serialization */
     PriorFactor() {}

    virtual ~PriorFactor() {}

    /** Constructor */
     PriorFactor(int key, const Eigen::VectorXd& prior, GaussianNoiseModel* model) :
        NoiseModelFactor1(model, key), prior_(prior)
    {
    }

    /** Convenience constructor that takes a full covariance argument */
     PriorFactor(int key, const Eigen::VectorXd& prior, const Eigen::MatrixXd& covariance) :
        NoiseModelFactor1(GaussianNoiseModel_Covariance(covariance), key), prior_(prior)
    {}


    /** implement functions needed for Testable */

    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Eigen::VectorXd>& x)const
    {
        std::map<int,Eigen::VectorXd>::const_iterator itb=x.find(key());
        return -(prior_-itb->second);//evaluateError(x);(prior_-x);
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Eigen::VectorXd>& x,std::vector<Eigen::MatrixXd>& H) const
    {
        std::map<int,Eigen::VectorXd>::const_iterator itb=x.find(key());
        return evaluateError(itb->second,*(H.begin()));
    }

    Eigen::VectorXd evaluateError(const Eigen::VectorXd& x) const
    {
        return -(prior_-x);
    }
    Eigen::VectorXd evaluateError(const Eigen::VectorXd& x, Eigen::MatrixXd& H) const
    {
        H=Eigen::MatrixXd::Identity(x.rows(),x.rows());
        Eigen::VectorXd h(x.rows());
        h=prior_-x;
        return -h;
    }




    const Eigen::VectorXd prior() const
    {
        return prior_;
    }


    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x1,
                                            const std::map<int,Eigen::VectorXd>& x2,
                                            std::vector<Eigen::MatrixXd>& H) const
    {
        Eigen::VectorXd uw(6);
        uw.setZero();
        return uw;
    }

#ifdef GMF_Using_Pose3
    //nonsense for virtual;
    virtual Eigen::VectorXd evaluateError(const Pose3 x) const
    {
        Eigen::VectorXd pp(6);
        pp.setZero();
        return pp;
    }
    //nonsense for virtual;
    virtual Eigen::VectorXd evaluateError(const Pose3 x, Eigen::MatrixXd& H) const
    {
        Eigen::VectorXd pp(6);
        pp.setZero();
        return pp;

    }
#else
    virtual Eigen::VectorXd evaluateError(const Pose2& x) const
    {
        Eigen::VectorXd xb;
        return xb;
    }
    virtual Eigen::VectorXd evaluateError(const Pose2& x, Eigen::MatrixXd& H) const
    {

        Eigen::VectorXd xb;
        return xb;

    }
#endif // GMF_Using_Pose3
    virtual NonlinearFactor* clone()const
    {
        PriorFactor* newfactor=new PriorFactor(key(),prior(),noiseModel_);
        return newfactor;
    }
};
};

#endif
