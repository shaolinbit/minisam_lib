#include <iostream>
#include "navigation/ImuBias.h"
#include "navigation/ImuFactor.h"
#include "nonlinear/ISAM2.h"
#include "nonlinear/NonlinearFactorGraph.h"
#include "slam/BetweenFactor.h"
#include "slam/PriorFactor.h"
#include "gmfconfig.h"
#include "inference/Symbol.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

using namespace std;
using namespace minisam;

/* ************************************************************************* */
int main(int argc, char* argv[])
{

    GaussianNoiseModel*  sigma_init_v=new IsotropicNoiseModel(3,1000);



    double AccelerometerSigma,GyroscopeSigma,IntegrationSigma,AccelerometerBiasSigma,GyroscopeBiasSigma, AverageDeltaT;
    int BodyPtx, BodyPty, BodyPtz, BodyPrx, BodyPry, BodyPrz;
    FILE *kittimetadatafile = fopen("data/KittiEquivBiasedImu_metadata.txt", "r");
    FILE *kittiIMU=fopen("data/KittiEquivBiasedImu.txt", "r");
    FILE *KittiGps=fopen("data/KittiGps_converted.txt", "r");
    FILE *fpstate=fopen("data/isam2Wholeresult.txt","w+");
    FILE *fprealtime=fopen("data/isam2realtimeb.txt","w+");
    double GPSTime,GPSX,GPSY,GPSZ;
    double IMUTime, IMUdt, IMUaccelX,IMUaccelY,IMUaccelZ,IMUomegaX,IMUomegaY,IMUomegaZ;
    int kittiindex=0;
    char kittimetadatabuf[400];
    int GPSskip=2;

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

    minivector GPSPrecisions(6);
    GPSPrecisions.data[0]=0.001;
    GPSPrecisions.data[1]=0.001;
    GPSPrecisions.data[2]=0.001;
    GPSPrecisions.data[3]=1.0/0.07;
    GPSPrecisions.data[4]=1.0/0.07;
    GPSPrecisions.data[5]=1.0/0.07;


    int firstGPSPose = 2;
    Pose3* currentPoseGlobal=new Pose3();//(Rot3(),firstGPSPosition);
    // Get initial conditions for the estimated trajectory

    minivector* currentVelocityGlobal=new minivector(3,0.0);

    ConstantBias* currentBias=new ConstantBias(0.0,0.0,0.0,0.0,0.0,0.0);

    minivector init_xPrecisions(6);
    init_xPrecisions.data[0]=1000.0;
    init_xPrecisions.data[1]=1000.0;
    init_xPrecisions.data[2]=1000.0;
    init_xPrecisions.data[3]=1.0;
    init_xPrecisions.data[4]=1.0;
    init_xPrecisions.data[5]=1.0;

    GaussianNoiseModel*  sigma_init_x=new GaussianNoiseModel(init_xPrecisions);

    minivector sigma_init_b_sigmas(6);
    sigma_init_b_sigmas.data[0]=0.10000;
    sigma_init_b_sigmas.data[1]=0.10000;
    sigma_init_b_sigmas.data[2]=0.10000;
    sigma_init_b_sigmas.data[3]=5.00e-05;
    sigma_init_b_sigmas.data[4]=5.00e-05;
    sigma_init_b_sigmas.data[5]=5.00e-05;

    GaussianNoiseModel*  sigma_init_b=new GaussianNoiseModel(sigma_init_b_sigmas);

    minivector sigma_between_b(6);
    sigma_between_b.data[0]=10*AccelerometerBiasSigma;
    sigma_between_b.data[1]=10*AccelerometerBiasSigma;
    sigma_between_b.data[2]=10*AccelerometerBiasSigma;
    sigma_between_b.data[3]=10*GyroscopeBiasSigma;
    sigma_between_b.data[4]=10*GyroscopeBiasSigma;
    sigma_between_b.data[5]=10*GyroscopeBiasSigma;

    GaussianNoiseModel*  noise_between_b=new GaussianNoiseModel(sigma_between_b);

    minivector g(3);
    g.data[0]=0.0;
    g.data[1]=0.0;
    g.data[2]=-9.8;

    minivector w_coriolis(3);
    minivector_set_zero(&w_coriolis);
    PreintegrationParams IMU_params(g);

    minivector iden3(AccelerometerSigma*AccelerometerSigma,
                                     AccelerometerSigma*AccelerometerSigma,AccelerometerSigma*AccelerometerSigma);

    minimatrix AccelerometerCovariance= minimatrix_vector_asDiagonal(iden3);
    IMU_params.setAccelerometerCovariance(AccelerometerCovariance);
    minivector_set_all(&iden3,GyroscopeSigma*GyroscopeSigma);
    minimatrix GyroscopeCovariance= minimatrix_vector_asDiagonal(iden3);
    IMU_params.setGyroscopeCovariance(GyroscopeCovariance);

    minivector_set_all(&iden3,IntegrationSigma*IntegrationSigma);
    minimatrix IntegrationCovariance= minimatrix_vector_asDiagonal(iden3);
    IMU_params.setIntegrationCovariance(IntegrationCovariance);
    IMU_params.setOmegaCoriolis(w_coriolis);

    ISAM2Params parameters;
    parameters.optimizationParamsGaussNewton=new ISAM2GaussNewtonParams;
    parameters.setFactorization("QR");
    parameters.relinearizeSkip = 1;
    ISAM2 isam(parameters);
    ISAM2Data isam2data;

    int kittigpsindex=firstGPSPose;
    int currentPoseKey,currentVelKey,currentBiasKey;
    double nowtime,previoustime;

    std::map<int,minimatrix*> newValues;

    NonlinearFactorGraph newFactors;
    minivector accMeas(3);
    minivector omegaMeas(3);
    int firstimupre=0;
    int updatecount=0;

    minivector temppose3dimv(3);

    GaussianNoiseModel* noiseModelGPS=NDiagonalNoiseModelPrecision(GPSPrecisions);

    while((!feof(KittiGps))&&(!feof(kittiIMU)))
    {
        fscanf(KittiGps, "%lf,%lf,%lf,%lf\n",&GPSTime,&GPSX, &GPSY, &GPSZ);
        currentPoseKey=Symbol('p',kittigpsindex).key();
        currentVelKey=Symbol('v',kittigpsindex).key();
        currentBiasKey=Symbol('b',kittigpsindex).key();

        nowtime=GPSTime;
        if(kittigpsindex==firstGPSPose)//??????
        {

            temppose3dimv.data[0]=GPSX;
            temppose3dimv.data[1]=GPSY;
            temppose3dimv.data[2]=GPSZ;

            Pose3* cupose=new Pose3(Rot3(),temppose3dimv);
            minivector* cvelocity=new minivector(*currentVelocityGlobal);
            ConstantBias* cbias=new  ConstantBias(*currentBias);

            newValues.insert(std::make_pair(currentPoseKey,cupose));
            minimatrix_memcpy(currentPoseGlobal,*cupose);
            newValues.insert(std::make_pair(currentVelKey, cvelocity));
            newValues.insert(std::make_pair(currentBiasKey, cbias));
            PriorFactor* npfp=new PriorFactor(currentPoseKey,new Pose3(*cupose), sigma_init_x);
            newFactors.push_back(npfp);
            PriorFactor* npfv=new PriorFactor(currentVelKey,new minivector(*cvelocity), sigma_init_v);
            newFactors.push_back(npfv);
            PriorFactor* npfcb=new PriorFactor(currentBiasKey,new ConstantBias(*cbias), sigma_init_b);
            newFactors.push_back(npfcb);
        }
        else
        {
            int imucount=0;

            PreintegratedImuMeasurements currentSummarizedMeasurement(IMU_params,*currentBias);
            if(firstimupre>0)
            {
                currentSummarizedMeasurement.integrateMeasurement(accMeas, omegaMeas, IMUdt);
            }
            for(;;)
            {
                firstimupre=1;
                fscanf(kittiIMU, "%lf %lf %lf %lf %lf %lf %lf %lf\n",&IMUTime, &IMUdt,  &IMUaccelX, &IMUaccelY, &IMUaccelZ, &IMUomegaX, &IMUomegaY, &IMUomegaZ);
                accMeas.data[0]=IMUaccelX;
                accMeas.data[1]=IMUaccelY;
                accMeas.data[2]=IMUaccelZ;
                omegaMeas.data[0]=IMUomegaX;
                omegaMeas.data[1]=IMUomegaY;
                omegaMeas.data[2]=IMUomegaZ;
                if(IMUTime>nowtime)
                {
                    break;
                }
                else
                {
                    imucount++;
                    currentSummarizedMeasurement.integrateMeasurement(accMeas, omegaMeas, IMUdt);
                }
            }

            ImuFactor* InteIMUFactor=new ImuFactor(currentPoseKey-1, currentVelKey-1, currentPoseKey, currentVelKey,
                                                   currentBiasKey, currentSummarizedMeasurement);

            newFactors.push_back(InteIMUFactor);

            BetweenFactor* ncbias=new BetweenFactor(currentBiasKey-1,currentBiasKey,new ConstantBias(),noise_between_b);
            newFactors.push_back(ncbias);


            minivector gpsposition(GPSX,GPSY,GPSZ);
            Pose3* GPSPose=new Pose3(currentPoseGlobal->rotation(),gpsposition);

            if(kittigpsindex%GPSskip==0)
            {
                PriorFactor* npfp3=new  PriorFactor(currentPoseKey,GPSPose, noiseModelGPS);
                newFactors.push_back(npfp3);
            }
            else
                delete GPSPose;
            newValues.insert(std::make_pair(currentPoseKey, new Pose3(*currentPoseGlobal)));
            newValues.insert(std::make_pair(currentVelKey,new minivector(*currentVelocityGlobal)));
            newValues.insert(std::make_pair(currentBiasKey, new ConstantBias(*currentBias)));

            if (kittigpsindex > (firstGPSPose + 2*GPSskip))
            {
                updatecount++;
                if(updatecount>=1000)
                {

                    for(std::vector<NoiseModelFactor*>::iterator bi=newFactors.begin();
                            bi!=newFactors.end(); bi++)
                    {
                        delete *bi;
                        *bi=NULL;
                    }

                    for(auto& dlt:newValues)
                    {
                        delete dlt.second;
                    }

                    isam.clearall();
                    isam2data.clearvalues();
                    isam2data.clearfactors();
                    delete parameters.optimizationParamsGaussNewton;


                    newFactors.clearall();
                    delete currentVelocityGlobal;
                    delete currentPoseGlobal;
                    delete noiseModelGPS;
                    delete  noise_between_b;
                    delete sigma_init_b;
                    delete  sigma_init_x;
                    delete currentBias;
                    delete  sigma_init_v;
                    return  0;
                }


                isam.update(newFactors, newValues,isam2data);


                newFactors.clear();
                newValues.clear();


                isam.calculateEstimate(isam2data);
                minimatrix_memcpy(currentPoseGlobal,isam2data.resulttheta_.at(currentPoseKey));
                minimatrix_memcpy(currentVelocityGlobal,isam2data.resulttheta_.at(currentVelKey));
                minimatrix_memcpy(currentBias,isam2data.resulttheta_.at(currentBiasKey));

                cout<<"isam.lastBacksubVariableCount: "<<isam.lastBacksubVariableCount<<endl;;
                fprintf(fprealtime,"%d %.15f %.15f %.15f %.15f %.15f %.15f %d\n",
                        kittigpsindex, currentPoseGlobal->translation().data[0], currentPoseGlobal->translation().data[1],
                        currentPoseGlobal->translation().data[2],
                        currentVelocityGlobal->data[0],currentVelocityGlobal->data[1],currentVelocityGlobal->data[2],
                        isam.lastBacksubVariableCount);

                   cout<<" "<<kittigpsindex<<" seconds has passed."<<endl;



                cout<<*currentPoseGlobal<<endl;
                minivector_print(*currentVelocityGlobal);
                cout<<endl;
                minimatrix_print(*currentBias);
                cout<<endl;


            }
        }
        previoustime=nowtime;
        kittigpsindex++;
    }

    for(std::vector<NoiseModelFactor*>::iterator bi=newFactors.begin();
            bi!=newFactors.end(); bi++)
    {
        delete *bi;
        *bi=NULL;
    }


    isam.clearall();
    isam2data.clearvalues();
    isam2data.clearfactors();
    delete parameters.optimizationParamsGaussNewton;


    newFactors.clearall();
    delete currentVelocityGlobal;
    delete currentPoseGlobal;
    delete noiseModelGPS;
    delete  noise_between_b;
    delete sigma_init_b;
    delete  sigma_init_x;
    delete currentBias;
    delete  sigma_init_v;
    return  0;
}

