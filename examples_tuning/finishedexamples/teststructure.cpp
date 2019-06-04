#include <iostream>
//#include "geometry/SimpleCamera.h"
#include "navigation/ImuBias.h"
#include "navigation/ImuFactor.h"
#include "nonlinear/ISAM2.h"
#include "nonlinear/NonlinearFactorGraph.h"
#include "slam/BetweenFactor.h"
#include "slam/PriorFactor.h"
#include "slam/PriorFactorPose3.h"
#include "gmfconfig.h"
//#include "slam/dataset.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include "inference/ClusterTree.h"

using namespace std;

struct isam2cliquedata
{

    std::vector<unsigned long>* children_;
    std::map<int,int>* assignedkeys;
    std::map<int,std::map<int,Eigen::VectorXd>::iterator>* solnPointers_;

};

class ISAM2Clique
{
public:
    GaussianConditional*  conditional_;
    //int parentindex_;
    ISAM2CliquePointer* parent_;
    GaussianFactorGraph*  cachedSeparatorMarginal_;
    //int currentindex_;
    std::vector<ISAM2CliquePointer*>* children_;
    //std::vector<ISAM2CliquePointer*> children_;
    std::map<int,int>* aaignedkeys;
    //std::map<int,unsigned long> solnPointers_;
  //sam2cliquedata dbb;
    int problemSize_;

   RealGaussianFactor*  cachedFactor_;
   // Eigen::Vector3d* gradientContribution_;
    std::map<int, std::map<int,Eigen::VectorXd>::iterator>* solnPointers_;
    //std::map<int, unsigned long> solnPointers_;

    bool setErased;

public:
ISAM2Clique()
{
conditional_=NULL;
 cachedSeparatorMarginal_=NULL;
    cachedFactor_=NULL;
    parent_=NULL;
    //solnPointers_.;
//    gradientContribution_=NULL;

}
~ISAM2Clique()
{
 if(conditional_!=NULL)
   {
      delete conditional_;
        conditional_=NULL;
    }
    if(cachedSeparatorMarginal_!=NULL)
    {
        delete cachedSeparatorMarginal_;
        cachedSeparatorMarginal_=NULL;
    }
    if(cachedFactor_!=NULL)
    {
        delete cachedFactor_;
        cachedFactor_=NULL;
    }
}
};
struct EliminationData
{
//EliminationData* parentdata;
ISAM2Clique* bnode_;
//std::vector<RealGaussianFactor*> childFactors(3,new RealGaussianFactor());
std::vector<int*> childFactors;
EliminationData(int* rf,ISAM2Clique* bnode):bnode_(bnode){
childFactors.push_back(rf);
//childFactors.reserve(5);
}
};

/* ************************************************************************* */
int main(int argc, char* argv[])
{
  //std::vector<unsigned long> datalist;
  //datalist.reserve(10000);
   int* rf0=new int(5);
   int* rf1=new int(6);
   ISAM2Clique* bnode0=new ISAM2Clique();
   ISAM2Clique* bnode1=new ISAM2Clique();
   EliminationData edata0(rf0,bnode0);
   EliminationData edata1(rf1,bnode1);
   delete edata1.bnode_;
   edata1.bnode_=NULL;
return 0;
}
