/**
 * @file Pose2SLAMExample.cpp
 * @brief A 2D Pose SLAM example
 * @date
 * @author
 */

/**
 * A simple 2D pose slam example
 *  - The robot moves in a 2 meter square
 *  - The robot moves 2 meters each step, turning 90 degrees after each step
 *  - The robot initially faces along the X axis (horizontal, to the right in 2D)
 *  - We have full odometry between pose
 *  - We have a loop closure constraint when the robot returns to the first position
 */

// In planar SLAM example we use Pose2 variables (x, y, theta) to represent the robot poses
#include "geometry/Pose2.h"

// We will use simple integer Keys to refer to the robot poses.
//#include <gtsam/inference/Key.h>

// In GTSAM, measurement functions are represented as 'factors'. Several common factors
// have been provided with the library for solving robotics/SLAM/Bundle Adjustment problems.
// Here we will use Between factors for the relative motion described by odometry measurements.
// We will also use a Between Factor to encode the loop closure constraint
// Also, we will initialize the robot at the origin using a Prior factor.
#include "linear/NoiseModel.h"
#include "nonlinear/NonlinearFactorGraph.h"
#include "slam/BetweenFactor.h"
#include "slam/PriorFactor.h"

// When the factors are created, we will add them to a Factor Graph. As the factors we are using
// are nonlinear factors, we will need a Nonlinear Factor Graph.

// Finally, once all of the factors have been added to our factor graph, we will want to
// solve/optimize to graph to find the best (Maximum A Posteriori) set of variable values.
// GTSAM includes several nonlinear optimizers to perform this step. Here we will use the
// a Gauss-Newton solver
#include "nonlinear/GaussNewtonOptimizer.h"

// Once the optimized values have been calculated, we can also calculate the marginal covariance
// of desired variables
#include "nonlinear/Marginals.h"

// The nonlinear solvers within GTSAM are iterative solvers, meaning they linearize the
// nonlinear functions around an initial linearization point, then solve the linear system
// to update the linearization point. This happens repeatedly until the solver converges
// to a consistent set of variable values. This requires us to specify an initial guess
// for each variable, held in a Values container.
using namespace std;
using namespace minisam;
int main(int argc, char** argv)
{
  // 1. Create a factor graph container and add factors to it
  NonlinearFactorGraph graph;

  // 2a. Add a prior on the first pose, setting it to the origin
  // A prior factor consists of a mean and a noise model (covariance matrix)
  //noiseModel::Diagonal::shared_ptr priorNoise = noiseModel::Diagonal::Sigmas(Vector3(0.3, 0.3, 0.1));

  minivector sigmanoiseq(3);
  sigmanoiseq.data[0]=0.3;
  sigmanoiseq.data[1]=0.3;
  sigmanoiseq.data[2]=0.1;
  GaussianNoiseModel* priorNoise=new GaussianNoiseModel(sigmanoiseq);
  Pose2* priorMean=new Pose2(0.0, 0.0, 0.0); // prior at origin
  int priorKey = 1;//This key must be in the poses key list.
  PriorFactor* npfp=new PriorFactor(priorKey, priorMean, priorNoise);
  graph.push_back(npfp);
  //graph.emplace_shared<PriorFactor<Pose2> >(1, Pose2(0, 0, 0), priorNoise);

  // For simplicity, we will use the same noise model for odometry and loop closures
 // noiseModel::Diagonal::shared_ptr model = noiseModel::Diagonal::Sigmas(Vector3(0.2, 0.2, 0.1));

  minivector sigmanoiseodm(3);
  sigmanoiseodm.data[0]=0.2;
  sigmanoiseodm.data[1]=0.2;
  sigmanoiseodm.data[2]=0.1;
  GaussianNoiseModel* model=new GaussianNoiseModel(sigmanoiseodm);

  // 2b. Add odometry factors
  // Create odometry (Between) factors between consecutive poses
  Pose2* odometryMeasurement1=new Pose2(2.0, 0.0, 0.0);
  BetweenFactor *nf1=new BetweenFactor(1, 2, odometryMeasurement1, model);
  graph.push_back(nf1);
  //graph.emplace_shared<BetweenFactor<Pose2> >(1, 2, Pose2(2, 0, 0     ), model);

  Pose2* odometryMeasurement2=new Pose2(2.0, 0.0, M_PI_2);
  BetweenFactor *nf2=new BetweenFactor(2, 3, odometryMeasurement2, model);
  graph.push_back(nf2);
  //graph.emplace_shared<BetweenFactor<Pose2> >(2, 3, Pose2(2, 0, M_PI_2), model);
  BetweenFactor *nf3=new BetweenFactor(3, 4, new Pose2(odometryMeasurement2), model);
  graph.push_back(nf3);
  //graph.emplace_shared<BetweenFactor<Pose2> >(3, 4, Pose2(2, 0, M_PI_2), model);
  BetweenFactor *nf4=new BetweenFactor(4, 5, new Pose2(odometryMeasurement2), model);
  graph.push_back(nf4);
 // graph.emplace_shared<BetweenFactor<Pose2> >(4, 5, Pose2(2, 0, M_PI_2), model);

  // 2c. Add the loop closure constraint
  // This factor encodes the fact that we have returned to the same pose. In real systems,
  // these constraints may be identified in many ways, such as appearance-based techniques
  // with camera images. We will use another Between Factor to enforce this constraint:
  BetweenFactor *nf5=new BetweenFactor(5, 2, new Pose2(odometryMeasurement2), model);
  graph.push_back(nf5);
 // graph.emplace_shared<BetweenFactor<Pose2> >(5, 2, Pose2(2, 0, M_PI_2), model);
 //graph.print("\nFactor Graph:\n"); // print

  // 3. Create the data structure to hold the initialEstimate estimate to the solution
  // For illustrative purposes, these have been deliberately set to incorrect values
  std::map<int,minimatrix*> initialEstimate;

  initialEstimate.insert(std::make_pair(1, new Pose2(0.5, 0.0,  0.2)));
  initialEstimate.insert(std::make_pair(2, new Pose2(2.3, 0.1, -0.2)));
  initialEstimate.insert(std::make_pair(3, new Pose2(4.1, 0.1,  M_PI_2)));
  initialEstimate.insert(std::make_pair(4, new Pose2(4.0, 2.0,  M_PI  )));
  initialEstimate.insert(std::make_pair(5, new Pose2(2.1, 2.1, -M_PI_2)));


  cout<<endl;
  cout<<"initialEstimate"<<endl;
  for(auto& kpairp:initialEstimate)
  {
      cout<<kpairp.first<<endl;
      minimatrix_print(kpairp.second);
      cout<<endl;
  }


  // 4. Optimize the initial values using a Gauss-Newton nonlinear optimizer
  // The optimizer accepts an optional set of configuration parameters,
  // controlling things like convergence criteria, the type of linear
  // system solver to use, and the amount of information displayed during
  // optimization. We will set a few parameters as a demonstration.
  GaussNewtonParams parameters;
  // Stop iterating once the change in error between steps is less than this value
  parameters.relativeErrorTol = 1e-5;
  // Do not perform more than N iteration steps
  parameters.maxIterations = 100;

  parameters.linearSolverType=NonlinearOptimizerParams::MULTIFRONTAL_QR;
  // Create the optimizer ...
  GaussNewtonOptimizer optimizer(graph, initialEstimate, parameters);
  // ... and optimize

  std::map<int,minimatrix*> result = optimizer.optimize();


  cout<<endl;
  cout<<"resultp"<<endl;
  for(auto& kpair:result)
  {
      cout<<kpair.first<<endl;
      minimatrix_print(kpair.second);
      cout<<endl;
  }

  // 5. Calculate and print marginal covariances for all variables
  cout.precision(3);
  Marginals marginals(graph, result);
  cout << "x1 covariance:\n" <<endl;
  minimatrix m1=marginals.marginalCovariance(1);
   minimatrix_print(m1);
  cout << "x2 covariance:\n" ;
  cout<< endl;
  minimatrix m2=marginals.marginalCovariance(2);
  minimatrix_print(m2);
  cout << "x3 covariance:\n" <<endl;
  minimatrix m3=marginals.marginalCovariance(3);
   minimatrix_print(m3);
   cout << endl;
  cout << "x4 covariance:\n" << endl;
    minimatrix m4=marginals.marginalCovariance(4);
   minimatrix_print(m4);
  cout<<endl;
  cout << "x5 covariance:\n" << endl;
      minimatrix m5=marginals.marginalCovariance(5);
   minimatrix_print(m5);
  cout << endl;


  graph.clearall();
  delete priorNoise;
  delete model;



  return 0;
}

