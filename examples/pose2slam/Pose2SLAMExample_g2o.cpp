
/**
 * @file Pose2SLAMExample_g2o.cpp
 * @brief A 2D Pose SLAM example that reads input from g2o, converts it to a factor graph and does the
 * optimization. Output is written on a file, in g2o format
 * Syntax for the script is ./Pose2SLAMExample_g2o input.g2o output.g2o
 */

#include "slam/dataset.h"
#include "slam/PriorFactor.h"
#include "nonlinear/GaussNewtonOptimizer.h"
#include <fstream>

using namespace std;
using namespace minisam;

// HOWTO: ./Pose2SLAMExample_g2o inputFile outputFile (maxIterations) (tukey/huber)
int main(const int argc, const char *argv[]) {

  string kernelType = "none";
  int maxIterations = 100; // default
  string g2oFile = "examples_tuning/data/noisyToyGraph.txt"; // default

  // Parse user's inputs
  if (argc > 1){
    g2oFile = argv[1]; // input dataset filename
    // outputFile = g2oFile = argv[2]; // done later
  }
  if (argc > 3){
    maxIterations = atoi(argv[3]); // user can specify either tukey or huber
  }
  if (argc > 4){
    kernelType = argv[4]; // user can specify either tukey or huber
  }

  // reading file and creating factor graph
  GraphAndValues readresult;
  bool is3D = false;
  if(kernelType.compare("none") == 0){
    readresult = readG2o(g2oFile,is3D);
  }
  if(kernelType.compare("huber") == 0){
    std::cout << "Using robust kernel: huber " << std::endl;
    readresult = readG2o(g2oFile,is3D, KernelFunctionTypeHUBER);
  }

  if(kernelType.compare("tukey") == 0){
    std::cout << "Using robust kernel: tukey " << std::endl;
    readresult = readG2o(g2oFile,is3D, KernelFunctionTypeTUKEY);
  }

  // Add prior on the pose having index (key) = 0
  NonlinearFactorGraph graphWithPrior = *(readresult.first);
  std::map<int,minimatrix*> initial=*(readresult.second);


  GaussianNoiseModel* priorModel = //
      GaussianNoiseModel::Variances(minivector(1e-6, 1e-6, 1e-8));
  graphWithPrior.push_back(new PriorFactor(0,new Pose2(),priorModel));

  //graphWithPrior.add(PriorFactor<Pose2>(0, Pose2(), priorModel));
  std::cout << "Adding prior on pose 0 " << std::endl;
  printf("initial Result:\n");
  for(auto& bf:initial)
  {
      cout<<bf.first<<endl;
      minimatrix_print(bf.second);
      cout<<endl;
  }
  GaussNewtonParams params;
  params.setVerbosity("TERMINATION");
  if (argc > 3) {
    params.maxIterations = maxIterations;
    std::cout << "User required to perform maximum  " << params.maxIterations << " iterations "<< std::endl;
  }

  std::cout << "Optimizing the factor graph" << std::endl;
  GaussNewtonOptimizer optimizer(graphWithPrior, initial, params);
  std::cout << "initial error=" <<graphWithPrior.error(initial)<< std::endl;

  std::map<int,minimatrix*> result = optimizer.optimize();
  std::cout << "Optimization complete" << std::endl;


  std::cout << "final error=" <<graphWithPrior.error(result)<< std::endl;


  if (argc < 3) {
    printf("Final Result:\n");
  for(auto& bf:result)
  {
      cout<<bf.first<<endl;
      minimatrix_print(bf.second);
      cout<<endl;
  }
    const string outputFile = "examples_tuning/data/noisyToyGraphresult.txt";
    std::cout << "Writing results to file: " << outputFile << std::endl;
    writeG2o(*readresult.first, result, outputFile);
    std::cout << "done! " << std::endl;
  } else {
    const string outputFile = argv[2];
    std::cout << "Writing results to file: " << outputFile << std::endl;
    //GraphAndValues readresult2 = readG2o(g2oFile);
    writeG2o(*readresult.first, result, outputFile);
    std::cout << "done! " << std::endl;
  }
  delete priorModel;
  delete readresult.first->at(0)->noiseModel_;

  graphWithPrior.clearall();
  delete readresult.first;
  delete readresult.second;
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

