#include <iostream>
#include "linear/NoiseModel.h"
#include "nonlinear/NonlinearFactorGraph.h"
#include "slam/BetweenFactor.h"
#include "slam/PriorFactor.h"
#include "nonlinear/ISAM2Params.h"
#include "gmfconfig.h"
// Finally, once all of the factors have been added to our factor graph, we will want to
// solve/optimize to graph to find the best (Maximum A Posteriori) set of variable values.
// GTSAM includes several nonlinear optimizers to perform this step. Here we will use the
// standard Levenberg-Marquardt solver
#include "nonlinear/LevenbergMarquardtOptimizer.h"

// Once the optimized values have been calculated, we can also calculate the marginal covariance
// of desired variables
#include "nonlinear/Marginals.h"
// Before we begin the example, we must create a custom unary factor to implement a
// "GPS-like" functionality. Because standard GPS measurements provide information
// only on the position, and not on the orientation, we cannot use a simple prior to
// properly model this measurement.
//
// The factor will be a unary factor, affect only a single system variable. It will
// also use a standard Gaussian noise model. Hence, we will derive our new factor from
// the NoiseModelFactor1.
#include "nonlinear/NonlinearFactor.h"
#include "slam/UnaryFactor.h"


using namespace minisam;

int main()
{
    // 1. Create a factor graph container and add factors to it
    NonlinearFactorGraph graph;

// 2a. Add odometry factors
    // For simplicity, we will use the same noise model for each odometry factor
    minivector sigmanoiseodm(3);
    sigmanoiseodm.data[0]=0.2;
    sigmanoiseodm.data[1]= 0.2;
    sigmanoiseodm.data[2]= 0.1;
    GaussianNoiseModel* odometryNoise1=new GaussianNoiseModel(sigmanoiseodm);
    GaussianNoiseModel* odometryNoise2=new GaussianNoiseModel(sigmanoiseodm);
    Pose2* odometryMeasurement1=new Pose2(2.0, 0.0, 0.0);
    // Create odometry (Between) factors between consecutive poses
    BetweenFactor* nf1=new BetweenFactor(1, 2, odometryMeasurement1, odometryNoise1);
    graph.push_back(nf1);
    BetweenFactor* nf2=new BetweenFactor(2, 3, new Pose2(*odometryMeasurement1), odometryNoise2);
    graph.push_back(nf2);
    // 2b. Add "GPS-like" measurements
    // We will use our custom UnaryFactor for this.
// noiseModel::Diagonal::shared_ptr unaryNoise = NoiseModel::Diagonal::Sigmas(Vector2(0.1, 0.1)); // 10cm std on x,y
    minivector sigmanoiseodm2(2);
    sigmanoiseodm2.data[0]=0.1;
    sigmanoiseodm2.data[1]=0.1;
    GaussianNoiseModel* unaryNoise1=new GaussianNoiseModel(sigmanoiseodm2);
    GaussianNoiseModel* unaryNoise2=new GaussianNoiseModel(sigmanoiseodm2);
    GaussianNoiseModel* unaryNoise3=new GaussianNoiseModel(sigmanoiseodm2);
    UnaryFactor* u1=new UnaryFactor(1,0.0,0.0,unaryNoise1);
    UnaryFactor* u2=new UnaryFactor(2,2.0,0.0,unaryNoise2);
    UnaryFactor* u3=new UnaryFactor(3,4.0,0.0,unaryNoise3);

    graph.push_back(u1);
    graph.push_back(u2);
    graph.push_back(u3);

    // 3. Create the data structure to hold the initialEstimate estimate to the solution
    // For illustrative purposes, these have been deliberately set to incorrect values
    std::map<int,minimatrix*> initialEstimate;
    initialEstimate.insert(std::make_pair(1,new  Pose2(0.5, 0.0, 0.2)));
    initialEstimate.insert(std::make_pair(2,new Pose2(2.3, 0.1, -0.2)));
    initialEstimate.insert(std::make_pair(3,new Pose2(4.1, 0.1, 0.1)));

    //std::map<int,minimatrix*> initialEstimate;
    // 4. Optimize using Levenberg-Marquardt optimization. The optimizer
    // accepts an optional set of configuration parameters, controlling
    // things like convergence criteria, the type of linear system solver
    // to use, and the amount of information displayed during optimization.
    // Here we will use the default set of parameters.  See the
    // documentation for the full set of parameters.
    LevenbergMarquardtOptimizer optimizer(graph,initialEstimate);


    std::map<int,minimatrix*> result = optimizer.optimize();
    cout<<endl;
    cout<<"result"<<endl;
    for(auto& kpair:result)
    {
        cout<<kpair.first<<endl;
        minimatrix_print(kpair.second);
        cout<<endl;
    }
    //result.print("Final Result:\n");

    // 5. Calculate and print marginal covariances for all variables
    Marginals marginals(graph, result);
    for(auto& btc:*(marginals.bayesTree_->nodesbtc))
    {
        if(btc.second->cachedSeparatorMarginal_!=NULL)
        {
            cout<<"btc.cachedSeparatorMarginal_.size(): "<<btc.second->cachedSeparatorMarginal_->size()<<endl;
        }
        else
        {
            cout<<"btc.cachedSeparatorMarginal_=NULL"<<endl;
        }
    }
    cout << "x1 covariance:\n" <<endl;
    minimatrix m1=marginals.marginalCovariance(1);
    minimatrix_print(m1);
    cout<< endl;
    for(auto& btc:*(marginals.bayesTree_->nodesbtc))
    {
        if(btc.second->cachedSeparatorMarginal_!=NULL)
        {
            cout<<"btc.cachedSeparatorMarginal_.size(): "<<btc.second->cachedSeparatorMarginal_->size()<<endl;
        }
        else
        {
            cout<<"btc.cachedSeparatorMarginal_=NULL"<<endl;
        }
    }
    cout << "x2 covariance:\n" <<endl;
    minimatrix m2=marginals.marginalCovariance(2);
    minimatrix_print(m2);
    cout<< endl;
    cout << "x3 covariance:\n" <<endl;
    minimatrix m3=marginals.marginalCovariance(3);
    minimatrix_print(m3);
    cout << endl;


    graph.clearall();
    marginals.clearlineargraph();
    delete odometryNoise1;
    delete odometryNoise2;
    delete unaryNoise1;
    delete unaryNoise2;
    delete unaryNoise3;


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


