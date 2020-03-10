/**
 * @file Pose3SLAMExample_initializePose3.cpp
 * @brief A 3D Pose SLAM example that reads input from g2o, and initializes the Pose3 using InitializePose3
 * Syntax for the script is ./Pose3SLAMExample_initializePose3 input.g2o output.g2o
 */

#include "slam/dataset.h"
#include "slam/BetweenFactor.h"
#include "slam/PriorFactor.h"
#include "nonlinear/GaussNewtonOptimizer.h"
#include <fstream>

using namespace std;
using namespace minisam;

int main(const int argc, const char *argv[]) {

  // Read graph from file
  string g2oFile;
  if (argc < 2)
    g2oFile = "examples_tuning/data/pose3example.txt";
  else
    g2oFile = argv[1];

  bool is3D = true;
  GraphAndValues readresults = readG2o(g2oFile, is3D);

  // Add prior on the first key
  NonlinearFactorGraph* graphWithPrior = readresults.first;
  minivector bn(6);
  bn.data[0]=1e-6;bn.data[1]=1e-6;bn.data[2]=1e-6;
  bn.data[3]=1e-4;bn.data[4]=1e-4;bn.data[5]=1e-4;

  GaussianNoiseModel* priorModel =GaussianNoiseModel::Variances(bn);


  int firstKey = 0;
    for(auto& key_value: *(readresults.second)) {
    std::cout << "Adding prior to g2o file " << std::endl;
    firstKey = key_value.first;
    graphWithPrior->push_back(new PriorFactor(firstKey, new Pose3(), priorModel));
    break;
  }

  std::cout << "Optimizing the factor graph" << std::endl;
  GaussNewtonParams params;
  std::cout << "initial error=" <<graphWithPrior->error(*readresults.second)<< std::endl;

  params.setVerbosity("TERMINATION"); // this will show info about stopping conditions
  GaussNewtonOptimizer optimizer(*graphWithPrior, *readresults.second, params);
  std::map<int,minimatrix*> result = optimizer.optimize();
  std::cout << "Optimization complete" << std::endl;

  std::cout << "final error=" <<graphWithPrior->error(result)<< std::endl;

  if (argc < 3) {
     printf("Final Result:\n");
  for(auto& bf:result)
  {
      cout<<bf.first<<endl;
      minimatrix_print(bf.second);
      cout<<endl;
  }
    const string outputFile = "examples_tuning/data/pose3exampleresult.txt";
    std::cout << "Writing results to file: " << outputFile << std::endl;
    writeG2o(*graphWithPrior, result, outputFile);
    std::cout << "done! " << std::endl;
  } else {
    const string outputFile = argv[2];
    std::cout << "Writing results to file: " << outputFile << std::endl;
    writeG2o(*graphWithPrior, result, outputFile);
    std::cout << "done! " << std::endl;
  }

  for(NoiseModelFactor* nf:*graphWithPrior)
  {
    delete nf->noiseModel_;
    nf->noiseModel_=NULL;
  }
  graphWithPrior->clearall();
  delete readresults.first;
  delete readresults.second;
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

