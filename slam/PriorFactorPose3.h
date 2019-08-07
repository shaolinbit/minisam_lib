#ifndef PRIORFACTORPOSE3_H
#define PRIORFACTORPOSE3_H


#pragma once
#include "../nonlinear/NonlinearFactor.h"
/**
 * A class for a  prior on Pose3
 * @addtogroup SLAM
 */
namespace minisam
{
class PriorFactorPose3: public NoiseModelFactor1
{
private:
    Pose3 prior_; /** The measurement */

public:

    /** default constructor - only use for serialization */
    PriorFactorPose3():NoiseModelFactor1(1) {}

    virtual ~PriorFactorPose3() {}

    /** Constructor */
    PriorFactorPose3(int key, const Pose3& prior,
                     GaussianNoiseModel* model) :
        NoiseModelFactor1(model, key,1), prior_(prior)
    {
    }

    /** Convenience constructor Pose3 takes a full covariance argument */
    PriorFactorPose3(int key, const Pose3& prior, const Eigen::MatrixXd& covariance) :
        NoiseModelFactor1(GaussianNoiseModel_Covariance(covariance), key,1), prior_(prior)
    {
    }

    /** implement functions needed to derive from Factor */

    /** vector of errors */

    virtual Eigen::VectorXd evaluateError(const Pose3& x) const
    {
        // Pose3 h=x.inverse()*prior_;
        Eigen::VectorXd v(6);
        Pose3 fb(x);
        v=-fb.LocalCoordinates(prior_);
        return v;
    }

    virtual Eigen::VectorXd evaluateError(const Pose3& x, Eigen::MatrixXd& H)const
    {
        H=Eigen::MatrixXd::Identity(6,6);
        Eigen::VectorXd v(6);
        v=-x.LocalCoordinates(prior_);

        return v;
    }


    const Pose3 & prior() const
    {
        return prior_;
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
//nonsense for virtual;
    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd& x) const
    {

        Eigen::VectorXd uw(6);
        uw.setZero();
        return uw;

    }
    //nonsense for virtual;
    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd& x, Eigen::MatrixXd& H)const
    {
        Eigen::VectorXd uw(6);
        uw.setZero();
        return uw;
    }
    //nonsense for virtual;
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x1,
                                            const std::map<int,Eigen::VectorXd>& x2,
                                            std::vector<Eigen::MatrixXd>& H) const
    {
        Eigen::VectorXd uw(6);
        uw.setZero();
        return uw;
    }
    virtual NonlinearFactor* clone()const
    {
        PriorFactorPose3* newfactor=new PriorFactorPose3(key(),prior(),noiseModel_);
        return newfactor;
    }

};
};
#endif // PRIORFACTORPOSE3_H
