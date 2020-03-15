#include <iostream>
#include "nonlinear/ExtendedKalmanFilter.h"
#include "slam/PriorFactor.h"
#include "slam/BetweenFactor.h"

using namespace std;
using namespace minisam;

int main()
{

  minivector* x_initial=new minivector(2,0.0);

  minivector sigmanoise(2,0.1);

  //DiagonalNoiseModel P_initial = DiagonalNoiseModel::Sigmas(sigmanoise);
  GaussianNoiseModel* P_initial=new GaussianNoiseModel(sigmanoise);
 // cout<<P_initial->thisR()<<endl;
  // Create Key for initial pose
  //Symbol x0('x',0);

  // Create an ExtendedKalmanFilter object
  ExtendedKalmanFilter ekf(0, x_initial, P_initial);



  // Now predict the state at t=1, i.e. argmax_{x1} P(x1) = P(x1|x0) P(x0)
  // In Kalman Filter notation, this is x_{t+1|t} and P_{t+1|t}
  // For the Kalman Filter, this requires a motion model, f(x_{t}) = x_{t+1|t)
  // Assuming the system is linear, this will be of the form f(x_{t}) = F*x_{t} + B*u_{t} + w
  // where F is the state transition model/matrix, B is the control input model,
  // and w is zero-mean, Gaussian white noise with covariance Q
  // Note, in some models, Q is actually derived as G*w*G^T where w models uncertainty of some
  // physical property, such as velocity or acceleration, and G is derived from physics
  //
  // For the purposes of this example, let us assume we are using a constant-position model and
  // the controls are driving the point to the right at 1 m/s. Then, F = [1 0 ; 0 1], B = [1 0 ; 0 1]
  // and u = [1 ; 0]. Let us also assume that the process noise Q = [0.1 0 ; 0 0.1].
  minivector u(2);
  u.data[0]=1.0;u.data[1]=0.0;

  minivector sigmanoiseq(2,0.1);

  //DiagonalNoiseModel Q = DiagonalNoiseModel::Sigmas(sigmanoiseq, false);
  GaussianNoiseModel* Q =new GaussianNoiseModel(sigmanoiseq);
  //cout<<Q->thisR()<<endl;

  //cout<<Q<<endl;


  // This simple motion can be modeled with a BetweenFactor
  // Create Key for next pose
  //Symbol x1('x',1);
  // Predict delta based on controls
  //Point2 difference(1,0);

  minivector* difference=new minivector(2);
  difference->data[0]=1.0;difference->data[1]=0.0;
  //difference<<1.0, 0.0;


  // Create Factor
  BetweenFactor* factor1=new BetweenFactor(0, 1, difference, Q);

  //cout<<Q->R()<<endl;
  //cout<<Q->thisR()<<endl;

GaussianFactorGraph lngraph;
lngraph.reserve(2);
  // Predict the new value with the EKF class
   ekf.predict(factor1,lngraph);

  lngraph.clearall();


  //traits<Point2>::Print(x1_predict, "X1 Predict");
  cout<<"x1_predict:"<<endl;
  minimatrix_print(ekf.x_);
  cout<<endl;




  // Now, a measurement, z1, has been received, and the Kalman Filter should be "Updated"/"Corrected"
  // This is equivalent to saying P(x1|z1) ~ P(z1|x1)*P(x1)
  // For the Kalman Filter, this requires a measurement model h(x_{t}) = \hat{z}_{t}
  // Assuming the system is linear, this will be of the form h(x_{t}) = H*x_{t} + v
  // where H is the observation model/matrix, and v is zero-mean, Gaussian white noise with covariance R
  //
  // For the purposes of this example, let us assume we have something like a GPS that returns
  // the current position of the robot. Then H = [1 0 ; 0 1]. Let us also assume that the measurement noise
  // R = [0.25 0 ; 0 0.25].

  minivector SigmaR(2,0.25);
  GaussianNoiseModel* R=new GaussianNoiseModel(SigmaR);



  // This simple measurement can be modeled with a PriorFactor
  minivector* z1=new minivector(2);
  z1->data[0]=1.0;
  z1->data[1]=0.0;
  PriorFactor* factor2=new PriorFactor(1, z1, R);


  // Update the Kalman Filter with the measurement
  ekf.update(factor2,lngraph);
   lngraph.clearall();
 // traits<Point2>::Print(x1_update, "X1 Update");
  cout<<"x1_update:"<<endl;
  minimatrix_print(ekf.x_);
  cout<<endl;
 // cout<<"ekf.x_:"<<endl<<ekf.x_<<endl;





  // Do the same thing two more times...
  // Predict
 // Symbol x2('x',2);
  minivector* difference3=new minivector(2);
  difference3->data[0]=1.0;
  difference3->data[1]=0.0;

  BetweenFactor* factor3=new BetweenFactor(1, 2, difference3, Q);
   ekf.predict(factor3,lngraph);
  lngraph.clearall();
 // traits<Point2>::Print(x2_predict, "X2 Predict");
  cout<<"x2_predict:"<<endl;
  minimatrix_print(ekf.x_);
  cout<<endl;

  // Update
  minivector* z2=new minivector(2);
  z2->data[0]=2.0;
  z2->data[1]=0.0;

  PriorFactor* factor4=new PriorFactor(2, z2, R);
  ekf.update(factor4,lngraph);
  lngraph.clearall();
  //traits<Point2>::Print(x2_update, "X2 Update");
  cout<<"x2_update:"<<endl;
  minimatrix_print(ekf.x_);
  cout<<endl;



  // Do the same thing one more time...
  // Predict
  //Symbol x3('x',3);
  minivector* difference4=new minivector(2);
  difference4->data[0]=1.0;
  difference4->data[1]=0.0;
  BetweenFactor* factor5=new BetweenFactor(2, 3, difference4, Q);
   ekf.predict(factor5,lngraph);
  lngraph.clearall();
  cout<<"x3_predict:"<<endl;
  minimatrix_print(ekf.x_);
  cout<<endl;

  // Update
  minivector* z3=new minivector(2);
  z3->data[0]=3.0;
  z3->data[1]=0.0;
  PriorFactor* factor6=new PriorFactor(3, z3, R);
   ekf.update(factor6,lngraph);
  lngraph.clearall();
   cout<<"x3_update:"<<endl;
   minimatrix_print(ekf.x_);
   cout<<endl;

    delete P_initial;
  delete factor1;
  delete Q;
  delete factor2;
  delete R;
  delete factor3;
  delete factor5;
  delete factor4;
  delete factor6;



  return 0;
}

