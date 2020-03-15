#include <iostream>


#include "linear/NoiseModel.h"
#include "linear/KalmanFilter.h"

using namespace std;
using namespace minisam;
int main()
{

    //Declare A Kalman filter using CHolesky and Hessian factor.
    KalmanFilter  kf(2,QR);

    //Declare the initial two dimentional state
    minivector x0(2,0.0);
    minimatrix PMatrix(2,2);
    minimatrix_set_identity(&PMatrix);
    minimatrix_scale(&PMatrix,0.01);


    cout<<"PMatrix:"<<endl;
    minimatrix_print(PMatrix);
    cout<<endl;

    GaussianDensity* gdp0(kf.init(x0,PMatrix));






    cout<<"gdp0.size():"<<gdp0->size()<<endl;
    cout<<"gdp0.keys().size():"<<gdp0->keys().size()<<endl;
    cout<<"gdp0.mean()";
    minivector_print(gdp0->mean());
    cout<<endl;
    cout<<"gdp0.information()";
    minimatrix_print(gdp0->information());
    cout<<endl;
    cout<<"gdp0.covariance()"<<endl;
    minimatrix_print(gdp0->covariance());
    cout<<endl;







    //Initialize State transition matrix;
    minimatrix FMatrix(2,2);
    minimatrix_set_identity(&FMatrix);

    //Initialize process noise
    minivector Q0(2,0.1);

    GaussianNoiseModel* DQ0=new GaussianNoiseModel(Q0);

    minimatrix QMatrix(2,2);
    minimatrix_set_identity(&QMatrix);
    minimatrix_scale(&QMatrix,0.01);

    minivector u(2);
    u.data[0]=1.0;u.data[1]=0.0;


    minimatrix B(2,2);
    minimatrix_set_identity(&B);


    cout<<"+++++++++++++++++++++++++++++++++++"<<"begin predict"<<endl;

    //Predict Process
    GaussianDensity* gdp1=kf.predictQ(gdp0,FMatrix,B,u,QMatrix);
     cout<<"gdp1.size():"<<gdp1->size()<<endl;
    cout<<"gdp1.keys().size():"<<gdp1->keys().size()<<endl;
    cout<<"gdp1.mean()";
    minivector_print(gdp1->mean());
    cout<<endl;
    cout<<"gdp1.information()";
    minimatrix_print(gdp1->information());
    cout<<endl;
    cout<<"gdp1.covariance()"<<endl;
    minimatrix_print(gdp1->covariance());
    cout<<endl;



    //Initialize H
    minimatrix HMatrix(2,2);
    minimatrix_set_identity(&HMatrix);

    //Initialize Measurement noise
    minivector R0(2,0.1);
    GaussianNoiseModel* DR0=new GaussianNoiseModel(R0);

    minivector Z1(2);
    Z1.data[0]=1.0;Z1.data[1]=0.0;

    minivector Z2(2);
    Z2.data[0]=2.0;Z2.data[1]=0.0;

    minivector Z3(2);
    Z3.data[0]=3.0;Z3.data[1]=0.0;

     if(gdp0->model_!=NULL)
     {
       delete gdp0->model_;
     }

    delete gdp0;
    gdp0=NULL;
    cout<<"+++++++++++++++++++++++++++++++++++"<<"begin update"<<endl;
    gdp0=kf.update(gdp1,HMatrix,Z1,DR0);




     cout<<"gdp0.size():"<<gdp0->size()<<endl;
    cout<<"gdp0.keys().size():"<<gdp0->keys().size()<<endl;
    cout<<"gdp0.mean()";
    minivector_print(gdp0->mean());
    cout<<endl;
    cout<<"gdp0.information()";
    minimatrix_print(gdp0->information());
    cout<<endl;
    cout<<"gdp0.covariance()"<<endl;
    minimatrix_print(gdp0->covariance());
    cout<<endl;
     //Predict Process

      if(gdp1->model_!=NULL)
     {
       delete gdp1->model_;
     }
    delete gdp1;
    gdp1=NULL;


    gdp1=kf.predictQ(gdp0,FMatrix,B,u,QMatrix);
     cout<<"gdp1.size():"<<gdp1->size()<<endl;
    cout<<"gdp1.keys().size():"<<gdp1->keys().size()<<endl;
    cout<<"gdp1.mean()";
    minivector_print(gdp1->mean());
    cout<<endl;
    cout<<"gdp1.information()";
    minimatrix_print(gdp1->information());
    cout<<endl;
    cout<<"gdp1.covariance()"<<endl;
    minimatrix_print(gdp1->covariance());
    cout<<endl;
    cout<<"+++++++++++++++++++++++++++++++++++"<<"begin update"<<endl;

     if(gdp0->model_!=NULL)
     {
       delete gdp0->model_;
     }

    delete gdp0;
    gdp0=NULL;


    gdp0=kf.update(gdp1,HMatrix,Z2,DR0);



    cout<<"gdp0.size():"<<gdp0->size()<<endl;
    cout<<"gdp0.keys().size():"<<gdp0->keys().size()<<endl;
    cout<<"gdp0.mean()";
    minivector_print(gdp0->mean());
    cout<<endl;
    cout<<"gdp0.information()";
    minimatrix_print(gdp0->information());
    cout<<endl;
    cout<<"gdp0.covariance()"<<endl;
    minimatrix_print(gdp0->covariance());
    cout<<endl;
    cout<<"+++++++++++++++++++++++++++++++++++"<<"begin predict"<<endl;

     if(gdp1->model_!=NULL)
     {
       delete gdp1->model_;
     }

    delete gdp1;
    gdp1=NULL;




       //Predict Process
     gdp1=kf.predictQ(gdp0,FMatrix,B,u,QMatrix);
     cout<<"gdp1.size():"<<gdp1->size()<<endl;
    cout<<"gdp1.keys().size():"<<gdp1->keys().size()<<endl;
    cout<<"gdp1.mean()";
    minivector_print(gdp1->mean());
    cout<<endl;
    cout<<"gdp1.information()";
    minimatrix_print(gdp1->information());
    cout<<endl;
    cout<<"gdp1.covariance()"<<endl;
    minimatrix_print(gdp1->covariance());
    cout<<endl;
    cout<<"+++++++++++++++++++++++++++++++++++"<<"begin update"<<endl;

     if(gdp0->model_!=NULL)
     {
       delete gdp0->model_;
     }

    delete gdp0;
    gdp0=NULL;





     gdp0=kf.update(gdp1,HMatrix,Z3,DR0);

     cout<<"gdp0.size():"<<gdp0->size()<<endl;
    cout<<"gdp0.keys().size():"<<gdp0->keys().size()<<endl;
    cout<<"gdp0.mean()";
    minivector_print(gdp0->mean());
    cout<<endl;
    cout<<"gdp0.information()";
    minimatrix_print(gdp0->information());
    cout<<endl;
    cout<<"gdp0.covariance()"<<endl;
    minimatrix_print(gdp0->covariance());
    cout<<endl;

    if(gdp0->model_!=NULL)
     {
       delete gdp0->model_;
     }

    delete gdp0;
    delete DQ0;
     if(gdp1->model_!=NULL)
     {
       delete gdp1->model_;
     }
    delete gdp1;
    delete DR0;


    return 0;
}

