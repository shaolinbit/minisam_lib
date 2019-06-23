#ifndef BETWEENFACTORPOSE2_H_INCLUDED
#define BETWEENFACTORPOSE2_H_INCLUDED


#pragma once
#include "../nonlinear/NonlinearFactor.h"
namespace minisam
{
/**
 * A class for a measurement predicted by "between(config[key1],config[key2])"
 * @tparam VALUE the Pose2 Value type
 * @addtogroup SLAM
 */
class BetweenFactorPose2: public NoiseModelFactor2
{

private:
    Pose2 measured_; /** The measurement */

public:


    /** default constructor - only use for serialization */
    BetweenFactorPose2() {}

    /** Constructor */
    BetweenFactorPose2(int key1, int key2, const Pose2& measured,
                       SharedNoiseModel* model) :
        NoiseModelFactor2(model, key1, key2), measured_(measured)
    {
    }

    ~BetweenFactorPose2() {}


    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2) const
    {
        Eigen::VectorXd hx = p2-p1; // h(x)
        return hx;
    }

    virtual  Eigen::VectorXd evaluateError(const Eigen::VectorXd& p1,
                                           const Eigen::VectorXd& p2, Eigen::MatrixXd& H1,Eigen::MatrixXd& H2) const
    {
        Eigen::VectorXd hx = p2-p1; // h(x)
        int rankn=p1.rows();
        H1=Eigen::MatrixXd(rankn,rankn);
        H1.setIdentity();
        H1=-(H1);
        H2=Eigen::MatrixXd(rankn,rankn);
        H2.setIdentity();
        return hx;
    }
    /** return the measured */
    const Pose2& measured() const
    {
        return measured_;
    }

    /** number of variables attached to this factor */
    int size() const
    {
        return 2;
    }

    Eigen::VectorXd unwhitenedError(const std::map<int,Pose2>& x) const
    {
        std::map<int,Pose2>::const_iterator itb1=x.find(key1());
        std::map<int,Pose2>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second);
    }

    Eigen::VectorXd evaluateError(const Pose2& p2) const
    {
        Pose2 h=p2->inverse()*measured();
        Eigen::VectorXd v(3);
        v=-Pose2::ChartAtOrigin::Local(h);
        return v;

    }

    Eigen::VectorXd unwhitenedError(const std::map<int,Pose2>& x,std::vector<Eigen::MatrixXd>& H) const
    {
        std::map<int,Pose2>::const_iterator itb1=x.find(key1());
        std::map<int,Pose2>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second,*(H.begin()),*(H.begin()+1));
    }
    Eigen::VectorXd unwhitenedError(const std::map<int,Pose2>& x1,const std::map<int,Eigen::VectorXd>& x2) const
    {
        std::map<int,Pose2>::const_iterator itb1=x1.find(key1());
        std::map<int,Pose2>::const_iterator itb2=x1.find(key2());

        return evaluateError(itb1->second,itb2->second);
    }

    Eigen::VectorXd unwhitenedError(const std::map<int,Pose2>& x1,const std::map<int,Eigen::VectorXd>& x2,
                                    std::vector<Eigen::MatrixXd>& H) const
    {
        std::map<int,Pose2>::const_iterator itb1=x1.find(key1());
        std::map<int,Pose2>::const_iterator itb2=x1.find(key2());

        return evaluateError(itb1->second,itb2->second,*(H.begin()),*(H.begin()+1));
    }
    Eigen::VectorXd evaluateError(const Pose2& x, Eigen::MatrixXd& H)const
    {
        H=Eigen::MatrixXd::Identity(3,3);
        Pose2 h=x.inverse()*measured();
        Eigen::VectorXd v(3);
        v=-Pose2::ChartAtOrigin::Local(h);
        return v;
    }
    Eigen::VectorXd
    evaluateError(const Pose2& X1, const Pose2& X2, Eigen::MatrixXd& H1,Eigen::MatrixXd& H2 ) const
    {
        Pose2 hx = X1->between(X2, H1, H2); // h(x)
        tx1=measured();
        return X1->LocalCoordinates(hx);
    }

    Eigen::VectorXd
    evaluateError(const Pose2& X1, const Pose2& X2) const
    {
        Pose2 hx = X1.between(X2); // h(x)
        tx1=measured();
        return tx1.LocalCoordinates(hx);
    }

    virtual NonlinearFactor* clone()const
    {
        BetweenFactorPose2* newfactor=new BetweenFactorPose2(key1(),key2(),measured(),noiseModel_);
        return newfactor;
    }

};
};

#endif // BETWEENFACTORPOSE2_H_INCLUDED
