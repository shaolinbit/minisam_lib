/**
 * @file BearingFactor.cpp
 * @brief A 2D Range SLAM example
 */

// Both relative poses and recovered trajectory poses will be stored as Pose2 objects
#include "geometry/Pose2.h"

// Each variable in the system (poses and landmarks) must be identified with a unique key.
// We can either use simple integer keys (1, 2, 3, ...) or symbols (X1, X2, L1).
// Here we will use Symbols
#include "inference/Symbol.h"

// We want to use iSAM2 to solve the range-SLAM problem incrementally
#include "nonlinear/ISAM2.h"

// iSAM2 requires as input a set set of new factors to be added stored in a factor graph,
// and initial guesses for any new variables used in the added factors
#include "nonlinear/NonlinearFactorGraph.h"

// We will use a non-liear solver to batch-inituialize from the first 150 frames
#include "nonlinear/LevenbergMarquardtOptimizer.h"

// Measurement functions are represented as 'factors'. Several common factors
// have been provided with the library for solving robotics SLAM problems.
#include "slam/PriorFactor.h"
#include "slam/BetweenFactor.h"
#include "slam/BearingFactor.h"
#include "slam/dataset.h"
#include "slam/BearingRangeFactor.h"
#include "linear/mEstimator.h"
#include "nonlinear/Marginals.h"
// Standard headers, added last, so we know headers above work on their own
#include <fstream>
#include <iostream>

using namespace std;
using namespace minisam;


// main
int main (int argc, char** argv)
{

    // Create the same factor graph as in PlanarSLAMExample
    int i1 = symbol('x',1).key();
    int i2 = symbol('x',2).key();
    int i3 = symbol('x',3).key();
    NonlinearFactorGraph  graph;
    Pose2* priorMean =new  Pose2(0.0, 0.0, 0.0);
    //% prior at origin
    GaussianNoiseModel* priorNoise=new GaussianNoiseModel(minivector(0.3,0.3,0.1));

    graph.push_back(new PriorFactor(i1, priorMean, priorNoise));

    //% add directly to graph
    Pose2* odometry =new  Pose2(2.0, 0.0, 0.0);
    GaussianNoiseModel* odometryNoise=new GaussianNoiseModel(minivector(0.2,0.2,0.1));
    graph.push_back(new BetweenFactor(i1, i2, odometry, odometryNoise));
    graph.push_back(new BetweenFactor(i2, i3, new Pose2(odometry), odometryNoise));

    //%% Except, for measurements we offer a choice
    int j1 = symbol('l',1).key();
    int j2 = symbol('l',2).key();
    double degrees = M_PI/180;

    GaussianNoiseModel* bearingModel = new GaussianNoiseModel(minivector(1,0.1));

    minivector bfv(2);
    bfv.data[0]=0.1;bfv.data[1]=0.2;
    GaussianNoiseModel* brNoise = new GaussianNoiseModel(bfv);


    graph.push_back(new BearingFactor2D(i1, j1, new Rot2(45*degrees), bearingModel));
    graph.push_back(new BearingFactor2D(i2, j1, new Rot2(90*degrees), bearingModel));

    graph.push_back(new BearingRangeFactor2d(i3, j2, new Rot2(90*degrees), 2, brNoise));

    //%% Initialize MCMC sampler with ground truth
    std::map<int,minimatrix*> sample;
    sample.insert(std::make_pair(i1, new Pose2(0,0,0)));
    sample.insert(std::make_pair(i2, new Pose2(2,0,0)));
    sample.insert(std::make_pair(i3, new Pose2(4,0,0)));
    minivector* vj1=new minivector(2);
    vj1->data[0]=2;vj1->data[1]=2;
    sample.insert(std::make_pair(j1,vj1));
    minivector* vj2=new minivector(2);
    vj2->data[0]=4;vj2->data[1]=2;
    sample.insert(std::make_pair(j2, vj2));
   LevenbergMarquardtOptimizer optimizer(graph,sample);

   std::map<int,minimatrix*> result = optimizer.optimize();

    //%% Calculate and plot Covariance Ellipses
    Marginals marginals(graph, result);


    cout<<endl;
    cout<<"result"<<endl;
    for(auto& kpair:result)
    {
        cout<<symbolChr(kpair.first)<<endl;
        cout<<symbolIndex(kpair.first)<<endl;
        minimatrix_print(kpair.second);
        cout<<endl;
    }

    minimatrix m1=marginals.marginalCovariance(i1);
    cout << "i1 covariance:\n" <<endl;
    minimatrix_print(m1);
    cout<< endl;
    cout << "i2 covariance:\n" <<endl;
    minimatrix m2=marginals.marginalCovariance(i2);
    minimatrix_print(m2);
    cout<< endl;

    graph.clearall();
    marginals.clearlineargraph();
    delete priorNoise;
    delete brNoise;
    delete odometryNoise;
    delete bearingModel;


  for(auto& kst:optimizer.state_->values)
     {
       if(kst.second!=NULL)
       {
         delete kst.second;
         kst.second=NULL;
       }
     }
     optimizer.state_->values.clear();



    return 0;
}

