/**
 * @file RangeISAMExample_plaza1.cpp
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
#include "slam/RangeFactor.h"
#include "slam/dataset.h"
#include "linear/mEstimator.h"

// Standard headers, added last, so we know headers above work on their own
#include <fstream>
#include <iostream>

using namespace std;
using namespace minisam;
//namespace NM = noiseModel;


// load the odometry
// DR: Odometry Input (delta distance traveled and delta heading change)
//    Time (sec)  Delta Dist. Trav. (m) Delta Heading (rad)
typedef pair<double, Pose2*> TimedOdometry;
list<TimedOdometry> readOdometry() {
  list<TimedOdometry> odometryList;
  string data_file = "examples_tuning/data/Plaza2_DR.txt";
  ifstream is(data_file.c_str());

  while (is) {
    double t, distance_traveled, delta_heading;
    is >> t >> distance_traveled >> delta_heading;
    odometryList.push_back(
        TimedOdometry(t, new Pose2(distance_traveled, 0, delta_heading)));
  }
  is.clear(); /* clears the end-of-file and error flags */
  return odometryList;
}

// load the ranges from TD
//    Time (sec)  Sender / Antenna ID Receiver Node ID  Range (m)
typedef std::tuple<double, size_t, double> RangeTriple;
vector<RangeTriple> readTriples() {
  vector<RangeTriple> triples;
  string data_file = "data/Plaza2_TD.txt";
  ifstream is(data_file.c_str());

  while (is) {
    double t, sender, receiver, range;
    is >> t >> sender >> receiver >> range;
    triples.push_back(RangeTriple(t, receiver, range));
  }
  is.clear(); /* clears the end-of-file and error flags */
  return triples;
}

// main
int main (int argc, char** argv) {

  // load Plaza2 data
  list<TimedOdometry> odometryvalues = readOdometry();
//  size_t M = odometry.size();

  vector<RangeTriple> triples = readTriples();
  size_t K = triples.size();

  // parameters
  size_t minK = 150; // minimum number of range measurements to process initially
  size_t incK = 25; // minimum number of range measurements to process after
  bool groundTruth = false;
  bool robust = true;

  // Set Noise parameters
  //Vector priorSigmas = Vector3(1,1,M_PI);
  minivector priorSigmas(1,1,M_PI);
  //Vector odoSigmas = Vector3(0.05, 0.01, 0.2);
  minivector odoSigmas(0.05, 0.01, 0.2);
  double sigmaR = 100; // range standard deviation
  GaussianNoiseModel* // all same type
  priorNoise = new GaussianNoiseModel(priorSigmas); //prior
  GaussianNoiseModel* odoNoise = new GaussianNoiseModel(odoSigmas); // odometry
  GaussianNoiseModel* gaussian =new IsotropicNoiseModel(1, sigmaR); // non-robust
  GaussianNoiseModel* tukey = RobustNoiseModel::Create(Tukey_mEstimator::Create(15), gaussian); //robust
  GaussianNoiseModel* rangeNoise = robust ? tukey : gaussian;

  // Initialize iSAM
  ISAM2 isam;

  ISAM2Data isam2data;

  // Add prior on first pose
  Pose2* pose0 =new Pose2(-34.2086489999201, 45.3007639991120,
      M_PI - 2.02108900000000);
  NonlinearFactorGraph newFactors;
  //newFactors.push_back(PriorFactor<Pose2>(0, pose0, priorNoise));
  newFactors.push_back(new PriorFactor(Symbol('P',0).key(), new Pose2(*pose0), priorNoise));
  std::map<int,minimatrix*> initial;
  initial.insert(std::make_pair(Symbol('P',0).key(), pose0));

  //  initialize points
  if (groundTruth) { // from TL file
    minivector* v1=new minivector(2);
    v1->data[0]=-68.9265;v1->data[1]= 18.3778;
    initial.insert(std::make_pair(Symbol('L', 1).key(), v1));
    //initial.insert(symbol('L', 1), Point2(-68.9265, 18.3778));
     minivector* v6=new minivector(2);
    v6->data[0]=-37.5805;v6->data[1]= 69.2278;
    initial.insert(std::make_pair(Symbol('L', 6).key(), v6));
    //initial.insert(symbol('L', 6), Point2(-37.5805, 69.2278));
    minivector* v0=new minivector(2);
    v0->data[0]=-33.6205;v0->data[1]= 26.9678;
    initial.insert(std::make_pair(Symbol('L', 0).key(), v0));
    //initial.insert(symbol('L', 0), Point2(-33.6205, 26.9678));
    minivector* v5=new minivector(2);
    v5->data[0]=1.7095;v5->data[1]= -5.8122;
    initial.insert(std::make_pair(Symbol('L', 5).key(), v5));
    //initial.insert(symbol('L', 5), Point2(1.7095, -5.8122));
  } else { // drawn from sigma=1 Gaussian in matlab version
    minivector* v1=new minivector(2);
    v1->data[0]=3.5784;v1->data[1]= 2.76944;
    initial.insert(std::make_pair(Symbol('L', 1).key(), v1));
    //initial.insert(symbol('L', 1), Point2(3.5784, 2.76944));
    minivector* v6=new minivector(2);
    v6->data[0]=-1.34989;v6->data[1]= 3.03492;
    initial.insert(std::make_pair(Symbol('L', 6).key(), v6));
    //initial.insert(symbol('L', 6), Point2(-1.34989, 3.03492));
    minivector* v0=new minivector(2);
    v0->data[0]=0.725404;v0->data[1]= -0.0630549;
    initial.insert(std::make_pair(Symbol('L', 0).key(), v0));
    //initial.insert(symbol('L', 0), Point2(0.725404, -0.0630549));
    minivector* v5=new minivector(2);
    v5->data[0]=0.714743;v5->data[1]=-0.204966;
    initial.insert(std::make_pair(Symbol('L', 5).key(), v5));
    //initial.insert(symbol('L', 5), Point2(0.714743, -0.204966));
  }

  // set some loop variables
  size_t i = 1; // step counter
  size_t k = 0; // range measurement counter
  bool initialized = false;
  Pose2* lastPose = new Pose2(*pose0);
  size_t countK = 0;
  cout<<*pose0<<endl;

  std::map<int,minimatrix*> isaminit;

  // Loop over odometry
  for(const TimedOdometry& timedOdometry:odometryvalues)
    {
    //--------------------------------- odometry loop -----------------------------------------
   // boost::tie(t, odometry) = timedOdometry;
   double t=timedOdometry.first;
   Pose2* odometry=timedOdometry.second;

  // cout<<*odometry<<endl;
    // add odometry factor
    newFactors.push_back(new BetweenFactor(Symbol('P',i-1).key(), Symbol('P',i).key(), odometry,new GaussianNoiseModel(odoSigmas)));

    // predict pose and add as initial estimate
    Pose2* predictedPose = lastPose->compose(*odometry);
   // lastPose = predictedPose;
    minimatrix_memcpy(lastPose,predictedPose);
    initial.insert(std::make_pair(Symbol('P',i).key(), predictedPose));

    //cout<<*predictedPose<<endl;

    // Check if there are range factors to be added
    while (k < K && t >= std::get<0>(triples[k])) {
      size_t j = std::get<1>(triples[k]);
      double range = std::get<2>(triples[k]);
      newFactors.push_back(new RangeFactor2D(Symbol('P',i).key(), Symbol('L', j).key(), range,rangeNoise));
      k = k + 1;
      countK = countK + 1;
    }


    /*
    for(auto& bid:initial)
    {
        cout<<symbolChr(bid.first)<<endl;
        cout<<symbolIndex(bid.first)<<endl;
        minimatrix_print(bid.second);
    }
    cout<<"Finish printding"<<endl;
    cout<<"----------------------------------------------"<<endl;
    */

    // Check whether to update iSAM 2
   // if ((k > minK) && (countK > incK))
    if(k>5)
        {


  /*
    cout<<"initial.size()"<<endl;
    cout<<initial.size()<<endl;
    for(auto& bid:initial)
    {
        cout<<symbolChr(bid.first)<<endl;
        cout<<symbolIndex(bid.first)<<endl;
        minimatrix_print(bid.second);
    }
    cout<<"Finish printding"<<endl;
    cout<<"----------------------------------------------"<<endl;
    */



      if (!initialized) { // Do a full optimize for first minK ranges
        LevenbergMarquardtOptimizer batchOptimizer(newFactors, initial);
        cout<<batchOptimizer.error()<<endl;
        batchOptimizer.optimize();
        initialized = true;



      }
      isam.update(newFactors, initial,isam2data);


      isam.calculateEstimate(isam2data);
      //lastPose = result.at<Pose2>(i);
      minimatrix_memcpy(lastPose,isam2data.resulttheta_.at(Symbol('P',i).key()));
      cout<<Symbol('P',i).chr()<<" "<<Symbol('P',i).index()<<" th value."<<endl;
      minimatrix_print(lastPose);


      newFactors.clear();
      initial.clear();

     // newFactors = NonlinearFactorGraph();
     // initial = Values();
      countK = 0;
    }
    i += 1;
    //--------------------------------- odometry loop -----------------------------------------
  }
        for(auto& kst:initial)
     {
       if(kst.second!=NULL)
       {
         delete kst.second;
         kst.second=NULL;
       }
     }

  exit(0);
}

