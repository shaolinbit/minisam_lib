/**
 *  @file   PseudorangeMaxMix.cpp
 *  @author Ryan
 *  @brief  Implementation file for pseudorange max-mix factor
 **/

#include "PseudorangeMaxMix.h"

using namespace std;

namespace minisam {



  Eigen::VectorXd PseudorangeMaxMix::evaluateError(const Eigen::VectorXd& q) const
  {

  Eigen::VectorXd hv(h_);
  double error=hv.dot(q)-measured_;
  Eigen::MatrixXd qg1(1,1);
  qg1<<1/hyp_;
  GaussianNoiseModel* g1=new GaussianNoiseModel(qg1);
    Eigen::MatrixXd qg2(1,1);
  qg2<<w_/hyp_;
  GaussianNoiseModel* g2=new GaussianNoiseModel(qg2);


      Eigen::VectorXd errorv(1);
      errorv<<error;

    double m1 = this->noiseModel()->distance(errorv);
    Eigen::MatrixXd info1(g1->information());

    double nu1 = 1.0/sqrt(info1.inverse().determinant());
    double l1 = nu1 * exp(-0.5*m1);

    double m2 = nullHypothesisModel_->distance(errorv);
    Eigen::MatrixXd info2(g2->information());
    double nu2 = 1.0/sqrt(info2.inverse().determinant());
    double l2 = nu2 * exp(-0.5*m2);



    if (l2>l1) {

      error *= sqrt(w_);
    }

    Eigen::VectorXd result(1);
    result<<error;
    return result;

  }
  Eigen::VectorXd PseudorangeMaxMix::evaluateError(const Eigen::VectorXd& q,Eigen::MatrixXd& H1 ) const
  {

  Eigen::VectorXd hv(h_);
  double error=hv.dot(q)-measured_;
  Eigen::MatrixXd qg1(1,1);
  qg1<<1/hyp_;
  GaussianNoiseModel* g1=new GaussianNoiseModel(qg1);
    Eigen::MatrixXd qg2(1,1);
  qg2<<w_/hyp_;
  GaussianNoiseModel* g2=new GaussianNoiseModel(qg2);


    Eigen::VectorXd errorv(1);
      errorv<<error;

    double m1 = this->noiseModel()->distance(errorv);
    Eigen::MatrixXd info1(g1->information());

    double nu1 = 1.0/sqrt(info1.inverse().determinant());
    double l1 = nu1 * exp(-0.5*m1);

    double m2 = nullHypothesisModel_->distance(errorv);
    Eigen::MatrixXd info2(g2->information());
    double nu2 = 1.0/sqrt(info2.inverse().determinant());
    double l2 = nu2 * exp(-0.5*m2);


   H1.resize(1,5);
   H1<<h_.transpose();


    if (l2>l1) {
      H1*=w_;
      error *= sqrt(w_);
    }

    Eigen::VectorXd result(1);
    result<<error;
    return result;;

  }

};
