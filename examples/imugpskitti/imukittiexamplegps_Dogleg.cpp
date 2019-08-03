#include <iostream>
#include "minisam/navigation/ImuBias.h"
#include "minisam/navigation/ImuFactor.h"
#include "minisam/nonlinear/ISAM2.h"
#include "minisam/nonlinear/NonlinearFactorGraph.h"
#include "minisam/slam/BetweenFactor.h"
#include "minisam/slam/PriorFactor.h"
#include "minisam/slam/PriorFactorPose3.h"
#include "minisam/gmfconfig.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

using namespace std;
using namespace minisam;

/* ************************************************************************* */
int main(int argc, char* argv[])
{
    double AccelerometerSigma,GyroscopeSigma,IntegrationSigma,AccelerometerBiasSigma,GyroscopeBiasSigma, AverageDeltaT;
    int BodyPtx, BodyPty, BodyPtz, BodyPrx, BodyPry, BodyPrz;
    FILE *kittimetadatafile = fopen("examples/imugpskitti/data/KittiEquivBiasedImu_metadata.txt", "r");
    FILE *kittiIMU=fopen("examples/imugpskitti/data/KittiEquivBiasedImu.txt", "r");
    FILE *KittiGps=fopen("examples/imugpskitti/data/KittiGps_converted.txt", "r");
    FILE *fpstate=fopen("examples/imugpskitti/data/isam2Wholeresult.txt","w+");
    FILE *fprealtime=fopen("examples/imugpskitti/data/isam2realtimeb.txt","w+");
    double GPSTime,GPSX,GPSY,GPSZ;
    double IMUTime, IMUdt, IMUaccelX,IMUaccelY,IMUaccelZ,IMUomegaX,IMUomegaY,IMUomegaZ;
    int kittiindex=0;
    char kittimetadatabuf[400];
    int GPSskip=2;
  //   clock_t isamtime1,isamtime2;

    if (!kittimetadatafile)
    {
        return -1;
    }
    else
    {
        while(!feof(kittimetadatafile))
        {
            if(kittiindex==0)
            {
                fscanf(kittimetadatafile,"%[^\n]",kittimetadatabuf);
                cout<<kittimetadatabuf<<endl;
            }
            else
            {
                fscanf(kittimetadatafile, "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf\n",&BodyPtx,&BodyPty, &BodyPtz, &BodyPrx, &BodyPry, &BodyPrz,
                       &AccelerometerSigma,&GyroscopeSigma,&IntegrationSigma,&AccelerometerBiasSigma,&GyroscopeBiasSigma,&AverageDeltaT);
            }

            kittiindex++;
        }
    }

    int  kittidataindex=0;
    if (!kittiIMU)
    {
        return -1;
    }
    else
    {
        fscanf(kittiIMU,"%[^\n]",kittimetadatabuf);
        cout<<kittimetadatabuf<<endl;
    }
    if (!KittiGps)
    {
        return -1;
    }
    else
    {
        fscanf(KittiGps,"%[^\n]",kittimetadatabuf);
    }

    //noiseModelGPS = noiseModel.Diagonal.Precisions([ [0;0;0]; 1.0/0.07 * [1;1;1] ]);
    Eigen::VectorXd GPSPrecisions(6);
    GPSPrecisions<<0.001,0.001,0.001,1.0/0.07,1.0/0.07,1.0/0.07;
    int firstGPSPose = 2;
    Pose3 currentPoseGlobal;//(Rot3(),firstGPSPosition);
    // Get initial conditions for the estimated trajectory
    //currentPoseGlobal = Pose3(Rot3, GPS_data(firstGPSPose).Position); // initial pose is the reference frame (navigation frame)
    Eigen::VectorXd currentVelocityGlobal(3);
    currentVelocityGlobal=Eigen::VectorXd::Zero(3);
    Eigen::VectorXd currentBias(6);
    currentBias=Eigen::VectorXd::Zero(6);
    Eigen::VectorXd init_xPrecisions(6);
    init_xPrecisions<<1000.0,1000.0,1000.0, 1,1,1;
    DiagonalNoiseModel*  sigma_init_x=new DiagonalNoiseModel(init_xPrecisions);
    DiagonalNoiseModel*  sigma_init_v=new IsotropicNoiseModel(3,1000);
    Eigen::VectorXd sigma_init_b_sigmas(6);
    sigma_init_b_sigmas<<0.100,0.100,0.100,5.00e-05,5.00e-05,5.00e-05;
    DiagonalNoiseModel*  sigma_init_b=new DiagonalNoiseModel(sigma_init_b_sigmas);
    Eigen::VectorXd sigma_between_b(6);
    sigma_between_b<<10*AccelerometerBiasSigma,10*AccelerometerBiasSigma,10*AccelerometerBiasSigma,
                    10*GyroscopeBiasSigma, 10*GyroscopeBiasSigma, 10*GyroscopeBiasSigma;
    Eigen::Vector3d g(0,0,-9.8);
    Eigen::VectorXd w_coriolis(3);
    w_coriolis<<0,0,0;
    PreintegrationParams IMU_params(g);
    IMU_params.setAccelerometerCovariance(AccelerometerSigma*AccelerometerSigma *Eigen::Matrix3d::Identity());
    IMU_params.setGyroscopeCovariance(GyroscopeSigma*GyroscopeSigma * Eigen::Matrix3d::Identity());
    IMU_params.setIntegrationCovariance(IntegrationSigma*IntegrationSigma * Eigen::Matrix3d::Identity());
    IMU_params.setOmegaCoriolis(w_coriolis);

    ISAM2Params parameters;
    parameters.optimizationParamsDogleg=new ISAM2DoglegParams();
    parameters.setFactorization("CHOLESKY");
    //parameters.relinearizeThresholdDouble = 0.01;
    parameters.relinearizeSkip = 1;
    ISAM2 isam(parameters);
    ISAM2Data isam2data;

    int kittigpsindex=firstGPSPose;
    int currentPoseKey,currentVelKey,currentBiasKey;
    double nowtime,previoustime;

    std::map<int,Eigen::VectorXd> newValuesV,resultV;
    std::map<int,Pose3> newValuesP,resultP;
    NonlinearFactorGraph newFactors;
   // Pose3* GPSPose;
    Eigen::Vector3d accMeas, omegaMeas;
    int firstimupre=0;
    int updatecount=0;

    while((!feof(KittiGps))&&(!feof(kittiIMU)))
    {
        fscanf(KittiGps, "%lf,%lf,%lf,%lf\n",&GPSTime,&GPSX, &GPSY, &GPSZ);
        //currentPoseKey=kittigpsindex+PoseConst;
        currentPoseKey=Symbol('p',kittigpsindex).key();
        //currentVelKey=kittigpsindex+VelConst;
        currentVelKey=Symbol('v',kittigpsindex).key();
        //currentBiasKey=kittigpsindex+BiasConst;
        currentBiasKey=Symbol('b',kittigpsindex).key();

        nowtime=GPSTime;
        if(kittigpsindex==firstGPSPose)//??????
        {
            currentPoseGlobal=Pose3(Rot3(),Eigen::Vector3d(GPSX,GPSY,GPSZ));
            newValuesP.insert(std::make_pair(currentPoseKey, currentPoseGlobal));
            newValuesV.insert(std::make_pair(currentVelKey, currentVelocityGlobal));
            newValuesV.insert(std::make_pair(currentBiasKey, currentBias));
            PriorFactorPose3* npfp=new PriorFactorPose3(currentPoseKey, currentPoseGlobal, sigma_init_x);
            newFactors.push_back(npfp);
            PriorFactor* npfv=new PriorFactor(currentVelKey, currentVelocityGlobal, sigma_init_v);
            newFactors.push_back(npfv);
            PriorFactor* npfcb=new PriorFactor(currentBiasKey, currentBias, sigma_init_b);
            newFactors.push_back(npfcb);
        }
        else
        {
            PreintegratedImuMeasurements currentSummarizedMeasurement(&IMU_params,ConstantBias(currentBias));
            if(firstimupre>0)
                currentSummarizedMeasurement.integrateMeasurement(accMeas, omegaMeas, IMUdt);
            for(;;)
            {
                firstimupre=1;
                fscanf(kittiIMU, "%lf %lf %lf %lf %lf %lf %lf %lf\n",&IMUTime, &IMUdt,  &IMUaccelX, &IMUaccelY, &IMUaccelZ, &IMUomegaX, &IMUomegaY, &IMUomegaZ);
                accMeas(0)=IMUaccelX;
                accMeas(1)=IMUaccelY;
                accMeas(2)=IMUaccelZ;
                omegaMeas(0)=IMUomegaX;
                omegaMeas(1)=IMUomegaY;
                omegaMeas(2)=IMUomegaZ;
                if(IMUTime>nowtime)
                {
                    break;
                }
                else
                {
                    currentSummarizedMeasurement.integrateMeasurement(accMeas, omegaMeas, IMUdt);
                }
            }
            ImuFactor* InteIMUFactor=new ImuFactor(currentPoseKey-1, currentVelKey-1, currentPoseKey, currentVelKey,
                                                   currentBiasKey, currentSummarizedMeasurement);
            newFactors.push_back(InteIMUFactor);
            //Bias evolution as given in the IMU metadata
            DiagonalNoiseModel*  noise_between_b=new DiagonalNoiseModel(sigma_between_b);
            BetweenFactor* ncbias=new BetweenFactor(currentBiasKey-1,currentBiasKey,ConstantBias().vector(),noise_between_b);
            newFactors.push_back(ncbias);
            Eigen::Vector3d gpsposition(GPSX,GPSY,GPSZ);
            Pose3 GPSPose(currentPoseGlobal.rotation(),gpsposition);
            if(kittigpsindex%GPSskip==0)
            {
                DiagonalNoiseModel* noiseModelGPS=NDiagonalNoiseModelPrecision(GPSPrecisions);
                PriorFactorPose3* npfp3=new  PriorFactorPose3(currentPoseKey,GPSPose, noiseModelGPS);
                newFactors.push_back(npfp3);
            }
            newValuesP.insert(std::make_pair(currentPoseKey, currentPoseGlobal));
            newValuesV.insert(std::make_pair(currentVelKey, currentVelocityGlobal));
            newValuesV.insert(std::make_pair(currentBiasKey, currentBias));
            if (kittigpsindex > (firstGPSPose + 2*GPSskip))
            {
                if(DEBUGSTATE)
                    cout<<newFactors.size()<<endl;
                updatecount++;
                if(updatecount>105)
                {
                    cout<<endl;
                }
                if(kittigpsindex>=30)
                {


                /*
                isam.clearall();
                isam2data.clearpose();
                isam2data.clearfactors();
                delete parameters.optimizationParamsGaussNewton;
                return 0;*/
                //for(auto& ipair:*isam.nodesbtc)
               // {
               //   cout<<ipair.first<<endl;
               // }
               // cout<<endl;
                }

                isam.updatetimefortune=kittigpsindex;
                isam.update(newFactors, newValuesV,newValuesP,isam2data);


                newFactors.clear();
                newValuesV.clear();
                newValuesP.clear();

              //  isamtime1=clock();

                resultV = isam.calculateEstimate(isam2data,&resultP);
               // isamtime2=clock();
              //  cout<<"calculateEstimate costs "<<(double)(isamtime2-isamtime1)/CLOCKS_PER_SEC<<" seconds"<<endl;

                currentPoseGlobal = resultP.at(currentPoseKey);
                currentVelocityGlobal = resultV.at(currentVelKey);
                currentBias = resultV.at(currentBiasKey);

                //resultP.clear();
                //resultV.clear();
                cout<<"isam.lastBacksubVariableCount: "<<isam.lastBacksubVariableCount<<endl;;
                fprintf(fprealtime,"%d %.15f %.15f %.15f %.15f %.15f %.15f %d\n",
                        kittigpsindex, currentPoseGlobal.translation()(0), currentPoseGlobal.translation()(1),
                        currentPoseGlobal.translation()(2),
                        currentVelocityGlobal(0),currentVelocityGlobal(1),currentVelocityGlobal(2),isam.lastBacksubVariableCount);
                cout<<" "<<kittigpsindex<<" seconds has passed."<<endl;

                cout<<"GPSTime is "<<GPSTime<<endl;

                cout<<currentPoseGlobal<<endl;
                cout<<currentVelocityGlobal<<endl;
                cout<<endl;

            }
            if(kittigpsindex%10==0)
            {
                //    cout<<kittigpsindex<<"seconds has passed."<<endl;
                //    cout<<"GPSTime is "<<GPSTime<<endl;

                //    cout<<currentPoseGlobal<<endl;
                //    cout<<currentVelocityGlobal<<endl;
                //    cout<<endl;
                if(kittigpsindex==50)
                {
                 //   isam.clearall();
            //   delete parameters.optimizationParamsGaussNewton;
              //      return 0;
                }
                if(kittigpsindex==470)
                {
                    int PoseKeyindex;
                    int VelKeyindex;
                    int poseindex=firstGPSPose;
                    Pose3 Posetrjectory;
                    Eigen::VectorXd Veltrjectory(3);

                    for(poseindex=firstGPSPose; poseindex<kittigpsindex; poseindex++)
                    {
                       // PoseKeyindex=poseindex+PoseConst;
                        PoseKeyindex=Symbol('p',poseindex).key();
                        //VelKeyindex=poseindex+VelConst;
                        VelKeyindex=Symbol('v',poseindex).key();
                        Posetrjectory=resultP.at(PoseKeyindex);
                        Veltrjectory=resultV.at(VelKeyindex);

                        fprintf(fpstate,"%.15f %.15f %.15f %.15f %.15f %.15f\n",
                                Posetrjectory.translation()(0), Posetrjectory.translation()(1),Posetrjectory.translation()(2),
                                Veltrjectory(0),Veltrjectory(1),Veltrjectory(2));
                    }

                }

            }
        }
        previoustime=nowtime;
        kittigpsindex++;
        resultP.clear();
        resultV.clear();
    }

    isam.clearall();
    isam2data.clearpose();
    isam2data.clearfactors();
    delete parameters.optimizationParamsDogleg;
    return 0;
}

