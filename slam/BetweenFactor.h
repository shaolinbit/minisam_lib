#ifndef BETWEENFACTOR_H
#define BETWEENFACTOR_H

/**
 *  @file  BetweenFactor.h
 *  @author
 **/
#pragma once
#include "../nonlinear/NonlinearFactor.h"
namespace minisam
{
/**
 * A class for a measurement predicted by "between(config[key1],config[key2])"
 * @tparam VALUE the Value type
 * @addtogroup SLAM
 */
class BetweenFactor: public NoiseModelFactor2
{

private:
    Eigen::VectorXd measured_; /** The measurement */

public:

    /** default constructor - only use for serialization */
     BetweenFactor() {}

    /** Constructor */
     BetweenFactor(int key1, int key2, const Eigen::VectorXd measured,
                  SharedNoiseModel* model) :
        NoiseModelFactor2(model, key1, key2), measured_(measured)
    {
    }

    ~BetweenFactor() {}


    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2) const
    {
        Eigen::VectorXd hx = p2-p1; // h(x)
        return (hx-measured_);
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x1,
                                            const std::map<int, Eigen::VectorXd>& x2,
                                            std::vector<Eigen::MatrixXd> &H) const
    {
        std::map<int,Eigen::VectorXd>::const_iterator itb1=x2.find(key1());
        std::map<int,Eigen::VectorXd>::const_iterator itb2=x2.find(key2());

        return evaluateError(itb1->second,itb2->second,
                             *(H.begin()),*(H.begin()+1));
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int, Pose3>& x1,
                                            const std::map<int, Eigen::VectorXd>& x2) const
    {
        std::map<int,Eigen::VectorXd>::const_iterator itb1=x2.find(key1());
        std::map<int,Eigen::VectorXd>::const_iterator itb2=x2.find(key2());

        return evaluateError(itb1->second,itb2->second);
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

        return (hx-measured_);
    }
    /** return the measured */
    const Eigen::VectorXd measured() const
    {
        return measured_;
    }

    /** number of variables attached to this factor */
    int size() const
    {
        return 2;
    }

}; // \class BetweenFactor
};

#endif // BETWEENFACTOR_H
