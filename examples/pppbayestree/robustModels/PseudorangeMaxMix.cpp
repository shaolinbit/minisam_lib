/**
 *  @file   PseudorangeMaxMix.cpp
 *  @author Ryan
 *  @brief  Implementation file for pseudorange max-mix factor
 **/

#include "PseudorangeMaxMix.h"

using namespace std;

namespace minisam
{



minivector PseudorangeMaxMix::evaluateError(const minimatrix* q) const
{

    double error;
    miniblas_ddot(*h_,minivector(q),&error);
    error-=measured_;


    minimatrix qg1(1,1);
    qg1.data[0]=1/hyp_;
    GaussianNoiseModel* g1=new GaussianNoiseModel(qg1);

    //Eigen::MatrixXd qg2(1,1);
    //qg2<<w_/hyp_;
    minimatrix qg2(1,1);
    qg2.data[0]=w_/hyp_;
    GaussianNoiseModel* g2=new GaussianNoiseModel(qg2);


    minivector errorv(1,error);
   // errorv.data[0]=error;
    //errorv<<error;

    double m1 = this->noiseModel()->distance(errorv);
    minimatrix info1=g1->information();
    //mini_int_vector p1(info1.size1);
    int p1[info1.size1];
    int signum;
    nr_ludcmp(&info1, p1, &signum);
    double info1_inverse_determinant = minilinalg_LU_det(&info1, signum);
    double nu1=1.0/sqrt(info1_inverse_determinant);
    // double nu1 = 1.0/sqrt(info1.inverse().determinant());
    double l1 = nu1 * exp(-0.5*m1);

    double m2 = nullHypothesisModel_->distance(errorv);
    minimatrix info2=g2->information();
    // mini_int_vector p2 (info2.size1);
     int p2[info2.size1];
    nr_ludcmp(&info2, p2, &signum);
    double info2_inverse_determinant = minilinalg_LU_det(&info2, signum);

    double nu2=1.0/sqrt(info2_inverse_determinant);
    // double nu2 = 1.0/sqrt(info2.inverse().determinant());
    double l2 = nu2 * exp(-0.5*m2);


    if (l2>l1)
    {

        error *= sqrt(w_);
    }

    minivector result(1,error);
    //result<<error;
    //result.data[0]=error;

    return result;

}
minivector PseudorangeMaxMix::evaluateError(const minimatrix* q,minimatrix& H1 ) const
{

// minivector hv(h_);
// double error=hv.dot(q)-measured_;

    double error;
    miniblas_ddot(*h_,minivector(q),&error);
    error-=measured_;

// Eigen::MatrixXd qg1(1,1);
// qg1<<1/hyp_;
    minimatrix qg1(1,1);
    qg1.data[0]=1/hyp_;
    GaussianNoiseModel g1(qg1);

    //Eigen::MatrixXd qg2(1,1);
    //qg2<<w_/hyp_;
    minimatrix qg2(1,1);
    qg2.data[0]=w_/hyp_;
    GaussianNoiseModel g2(qg2);


    minivector errorv(1);
    errorv.data[0]=error;
    //errorv<<error;

    double m1 = this->noiseModel()->distance(errorv);
    minimatrix info1=g1.information();

    //mini_int_vector p1(info1.size1);
    int p1[info1.size1];
    int signum;
    nr_ludcmp(&info1, p1, &signum);
    double info1_inverse_determinant = minilinalg_LU_det(&info1, signum);

    double nu1=1.0/sqrt(info1_inverse_determinant);
    // double nu1 = 1.0/sqrt(info1.inverse().determinant());
    double l1 = nu1 * exp(-0.5*m1);

    double m2 = nullHypothesisModel_->distance(errorv);
    minimatrix info2=g2.information();

    //mini_int_vector p2(info2.size1);
    int p2[info2.size1];
    nr_ludcmp(&info2, p2, &signum);
    double info2_inverse_determinant = minilinalg_LU_det(&info2, signum);

    double nu2=1.0/sqrt(info2_inverse_determinant);
    // double nu2 = 1.0/sqrt(info2.inverse().determinant());
    double l2 = nu2 * exp(-0.5*m2);


    //H1.resize(1,5);
    //H1<<h_.transpose();
    minimatrix_resize(&H1,1,5);
    H1.data[0]=h_->data[0];
    H1.data[1]=h_->data[1];
    H1.data[2]=h_->data[2];
    H1.data[3]=h_->data[3];
    H1.data[4]=h_->data[4];



    if (l2>l1)
    {
        // H1*=w_;
        minimatrix_scale(&H1,w_);
        error *= sqrt(w_);
    }

    minivector result(1,error);
    //result<<error;
   // result.data[0]=error;


    return result;

}

};
