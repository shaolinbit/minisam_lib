#include <iostream>
#include "linear/NoiseModel.h"
#include "nonlinear/NonlinearFactorGraph.h"
#include "nonlinear/ISAM2.h"
#include "geometry/Cal3_S2Stereo.h"
#include "slam/PriorFactor.h"
#include "slam/ProjectionFactor.h"
#include "SFMdata.h"
#include "gmfconfig.h"

using namespace minisam;

int main()
{

   // Define the camera calibration parameters
  Cal3_S2* K(new Cal3_S2(50.0, 50.0, 0.0, 50.0, 50.0));

  // Define the camera observation noise model

  minivector posenoisesigma(6);
  posenoisesigma.data[0]=0.3;posenoisesigma.data[1]=0.3;posenoisesigma.data[2]=0.3;
  posenoisesigma.data[3]=0.1;posenoisesigma.data[4]=0.1;posenoisesigma.data[5]=0.1;

   GaussianNoiseModel* poseNoise =new GaussianNoiseModel(posenoisesigma);
 // Add a prior on landmark l0
      IsotropicNoiseModel* pointnoise = new IsotropicNoiseModel(3, 0.1);

  // Create the set of ground-truth landmarks
  std::vector<minivector*> points = createPoints();

  // Create the set of ground-truth poses
  std::vector<Pose3> poses = createPoses();

 // Create an iSAM2 object. Unlike iSAM1, which performs periodic batch steps to maintain proper linearization
  // and efficient variable ordering, iSAM2 performs partial relinearization/reordering at each step. A parameter
  // structure is available that allows the user to set various properties, such as the relinearization threshold
  // and type of linear solver. For this example, we we set the relinearization threshold small so the iSAM2 result
  // will approach the batch result.
  ISAM2Params parameters;
  parameters.optimizationParamsGaussNewton=new ISAM2GaussNewtonParams;
  parameters.relinearizeThresholdDouble = 0.01;
  parameters.relinearizeSkip = 1;
  ISAM2 isam(parameters);
  parameters.setFactorization("QR");
  ISAM2Data isam2data;

  // Create a Factor Graph and Values to hold the new data
  NonlinearFactorGraph graph;
  NonlinearFactorGraph nullgraph;
  std::map<int,minimatrix*>nullEstimate;

  std::map<int,minimatrix*> initialEstimate;
   IsotropicNoiseModel* isnnoise = //
      new IsotropicNoiseModel(2, 1.0); // one pixel in u and v
  // Loop over the different poses, adding the observations to iSAM incrementally
  for (int i = 0; i < poses.size(); ++i)
    {

    // Add factors for each landmark observation
    for (int j = 0; j < points.size(); ++j)
        {
      // Create ground truth measurement
      SimpleCamera camera(poses[i], *K);
      minivector measurement = camera.projectPoint(*(points[j]),NULL,NULL,NULL);


      // Add measurement
      GenericProjectionFactor* nf=new GenericProjectionFactor(measurement, isnnoise,Symbol('p',i).key(),Symbol('O',j).key(),K);
              // i+PoseConst, j+PointConst, K);
      graph.push_back(nf);
    }

    cout<<"graph.size(): "<<graph.size()<<endl;

    // Intentionally initialize the variables off from the ground truth
    Pose3 noise(Rot3::Rodrigues(-0.1, 0.2, 0.25), minivector(0.05, -0.10, 0.20));
    //Pose3 initial_xi = poses[i].compose(noise);
    Pose3* initial_xi =new Pose3(poses[i].multiply(noise));

    // Add an initial guess for the current pose
    initialEstimate.insert(std::make_pair(Symbol('p',i).key(),initial_xi));//i+PoseConst, initial_xi));

    // If this is the first iteration, add a prior on the first pose to set the coordinate frame
    // and a prior on the first landmark to set the scale
    // Also, as iSAM solves incrementally, we must wait until each is observed at least twice before
    // adding it to iSAM.
    if (i == 0)
        {
      // Add a prior on pose x0, with 30cm std on x,y,z 0.1 rad on roll,pitch,yaw

      cout<<"poseNoise.thisR():"<<endl;
      minimatrix_print(poseNoise->thisR());
      cout<<endl;
      PriorFactor* npfp=new PriorFactor(Symbol('p',0).key(), new Pose3(poses[0]), poseNoise);//(0+PoseConst, poses[0], poseNoise);


       minimatrix_print(poses[0].matrix());
       cout<<endl;
      graph.push_back(npfp);

       cout<<"pointnoise.thisR():"<<endl;
       minimatrix_print(pointnoise->thisR());
       cout<<endl;
      PriorFactor* np=new PriorFactor(Symbol('O',0).key(), new minivector(points[0]), pointnoise);//(0+PointConst, points[0], pointnoise);
      graph.push_back(np);

      // Add initial guesses to all observed landmarks
      minivector noise(3);
      noise.data[0]=-0.25;noise.data[1]=0.2;noise.data[2]=0.15;
      //(-0.25, 0.20, 0.15);
      for (int j = 0; j < points.size(); ++j)
        {
        // Intentionally initialize the variables off from the ground truth
        minivector* initial_lj=new minivector(3);
        minivector_add(initial_lj,*points[j],noise);
        initialEstimate.insert(std::make_pair(Symbol('O',j).key(), initial_lj));//(j+PointConst, initial_lj));
      }

    }
    else
     {
      // Update iSAM with the new factors

      isam.update(graph, initialEstimate,isam2data);

      isam.update(nullgraph, nullEstimate,isam2data);

     isam.calculateEstimate(isam2data);

   cout<<"result  map"<<endl;
    for(auto& pv:isam2data.resulttheta_)
   {
       cout<<pv.first<<endl;
       minimatrix_print(pv.second);
       cout<<endl;
   }


      cout << "****************************************************" << endl;
      cout << "Frame " << i << ": " << endl;
      if(i==2)
      {
        cout<<endl;
      }
      // Clear the factor graph and values for the next iteration


      graph.clear();
      initialEstimate.clear();


    }
  }
   delete pointnoise;
      delete poseNoise;
      delete parameters.optimizationParamsGaussNewton;
      isam.clearall();
  isam2data.clearvalues();
  isam2data.clearfactors();
  SFMData_ClearPoints(points);
  delete isnnoise;
  delete K;
    return 0;
}


