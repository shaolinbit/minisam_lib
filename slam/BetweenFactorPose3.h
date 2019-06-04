#ifndef BETWEENFACTORPOSE3_H_INCLUDED
#define BETWEENFACTORPOSE3_H_INCLUDED


/**
 *  @file  BetweenFactor.h
 *  @author Frank Dellaert, Viorela Ila
 **/
#pragma once
#include "../nonlinear/NonlinearFactor.h"

/**
 * A class for a measurement predicted by "between(config[key1],config[key2])"
 * @tparam VALUE the Value type
 * @addtogroup SLAM
 */
class BetweenFactorPose3: public NoiseModelFactor2
{

private:
    Pose3 measured_; /** The measurement */

public:


    /** default constructor - only use for serialization */
    BetweenFactorPose3() {}

    /** Constructor */
    BetweenFactorPose3(int key1, int key2, const Pose3& measured,
                       SharedNoiseModel* model) :
        NoiseModelFactor2(model, key1, key2), measured_(measured)
    {
    }

    ~BetweenFactorPose3() {}


    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2) const
    {
        Eigen::VectorXd hx = p2-p1; // h(x)
        // manifold equivalent of h(x)-z -> log(z,h(x))
        // cout<<p2<<endl;
        //  cout<<p1<<endl;
        // cout<<measured_<<endl;
        return hx;
        //return (hx-measured_);
    }

    virtual  Eigen::VectorXd evaluateError(const Eigen::VectorXd& p1,
                                           const Eigen::VectorXd& p2, Eigen::MatrixXd& H1,Eigen::MatrixXd& H2) const
    {
        Eigen::VectorXd hx = p2-p1; // h(x)
        //  cout<<p2<<endl;
        //  cout<<p1<<endl;
        // cout<<measured_<<endl;
        int rankn=p1.rows();
        H1=Eigen::MatrixXd(rankn,rankn);
        H1.setIdentity();
        H1=-(H1);
        H2=Eigen::MatrixXd(rankn,rankn);
        H2.setIdentity();
        //cout<<H1<<endl;
        return hx;
        // return (hx-measured_);
    }
    /** return the measured */
    const Pose3& measured() const
    {
        return measured_;
    }

    /** number of variables attached to this factor */
    int size() const
    {
        return 2;
    }
//#ifdef GMF_Using_Pose3
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x,std::vector<Eigen::MatrixXd>& H) const
    {
        std::map<int,Pose3>::const_iterator itb1=x.find(key1());
        std::map<int,Pose3>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second,*(H.begin()),*(H.begin()+1));
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x1,
                                            const std::map<int,Eigen::VectorXd>& x2,
                                            std::vector<Eigen::MatrixXd>& H) const
    {
        Eigen::VectorXd uw(6);
        uw.setZero();
        return uw;
    }
    Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x) const
    {
        std::map<int,Pose3>::const_iterator itb1=x.find(key1());
        std::map<int,Pose3>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second);
    }

    Eigen::VectorXd
    evaluateError(const Pose3& X1, const Pose3& X2, Eigen::MatrixXd& H1,Eigen::MatrixXd& H2 ) const
    {
       // Pose3 tx1=X1;
       // Pose3 hx = tx1.between(*X2, H1, H2); // h(x)
        Pose3 hx = X1.between(X2, H1, H2); // h(x)

      Pose3  tx1=measured();
       // cout<<"measured_"<<tx1<<endl;
        return tx1.LocalCoordinates(hx);
    }

    Eigen::VectorXd
    evaluateError(const Pose3& X1, const Pose3& X2) const
    {
       // Pose3 tx1=X1;
        //Pose3 hx = tx1.between(X2); // h(x)
        Pose3 hx = X1.between(X2); // h(x)
       Pose3 tx1=measured();
       // cout<<measured()<<endl;
        return tx1.LocalCoordinates(hx);
        // manifold equivalent of h(x)-z -> log(z,h(x))
        // return traits<T>::Local(measured_, hx);
    }
    virtual NonlinearFactor* clone()const
    {
        BetweenFactorPose3* newfactor=new BetweenFactorPose3(key1(),key2(),measured(),noiseModel_);
        return newfactor;
    }
};


#endif // BETWEENFACTORPOSE3_H_INCLUDED
