#ifndef BETWEENFACTORPOSE2_H_INCLUDED
#define BETWEENFACTORPOSE2_H_INCLUDED


/**
 *  @file  BetweenFactor.h
 *  @author
 **/
#pragma once
#include "../nonlinear/NonlinearFactor.h"

/**
 * A class for a measurement predicted by "between(config[key1],config[key2])"
 * @tparam VALUE the Value type
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

    /** implement functions needed to derive from Factor */

    /** vector of errors
    Eigen::VectorXd evaluateError(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2, boost::optional<Matrix&> H1 =
      boost::none, boost::optional<Matrix&> H2 = boost::none) const {
      T hx = traits<T>::Between(p1, p2, H1, H2); // h(x)
      // manifold equivalent of h(x)-z -> log(z,h(x))
    #ifdef SLOW_BUT_CORRECT_BETWEENFACTOR
      typename traits<T>::ChartJacobian::Jacobian Hlocal;
      Vector rval = traits<T>::Local(measured_, hx, boost::none, (H1 || H2) ? &Hlocal : 0);
      if (H1) *H1 = Hlocal * (*H1);
      if (H2) *H2 = Hlocal * (*H2);
      return rval;
    #else
      return traits<T>::Local(measured_, hx);
    #endif
    } */


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
        // cout<<"keys"<<endl;
        //cout<<key1()<<endl;
        //cout<<key2()<<endl;
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
        //Pose2 tx1=*X1;
        //Pose2 hx = tx1.between(*X2, H1, H2); // h(x)
        Pose2 hx = X1->between(X2, H1, H2); // h(x)
        tx1=measured();
        //  cout<<"measured_"<<tx1<<endl;
        //return tx1.LocalCoordinates(hx);
        return X1->LocalCoordinates(hx);
    }

    Eigen::VectorXd
    evaluateError(const Pose2& X1, const Pose2& X2) const
    {
        //Pose2 tx1=X1;
       // Pose2 hx = tx1.between(*X2); // h(x)
        Pose2 hx = X1.between(X2); // h(x)
        tx1=measured();
        //return tx1.LocalCoordinates(hx);
        return tx1.LocalCoordinates(hx);
    }

    virtual NonlinearFactor* clone()const
    {
        BetweenFactorPose2* newfactor=new BetweenFactorPose2(key1(),key2(),measured(),noiseModel_);
        return newfactor;
    }

};


#endif // BETWEENFACTORPOSE2_H_INCLUDED
