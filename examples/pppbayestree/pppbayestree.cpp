/*
 * @file pppBayesTree.cpp
 * @brief Iterative GPS Range/Phase Estimator with collected data
 * @author Ryan Watson & Jason Gross
 */


#include "nonlinear/ISAM2.h"
#include "nonlinear/NonlinearFactorGraph.h"
#include "inference/Symbol.h"

#include "slam/PriorFactor.h"
#include "slam/BetweenFactor.h"
#include "pppbayestree/gnssNavigation/GnssData.h"
#include "pppbayestree/gnssNavigation/GnssTools.h"
#include "pppbayestree/gnssNavigation/PhaseFactor.h"
#include "pppbayestree/gnssNavigation/nonBiasStates.h"
#include "pppbayestree/configReader/ConfDataReader.hpp"

#include "pppbayestree/gnssNavigation/PseudorangeFactor.h"


// GPSTK
#include "pppbayestree/gpstk/MJD.hpp"
#include "pppbayestree/gpstk/PowerSum.hpp"
#include "pppbayestree/gpstk/Decimate.hpp"
#include "pppbayestree/gpstk/SolidTides.hpp"
#include "pppbayestree/gpstk/PoleTides.hpp"
#include "pppbayestree/gpstk/TropModel.hpp"
#include "pppbayestree/gpstk/BasicModel.hpp"
#include "pppbayestree/gpstk/CommonTime.hpp"
#include "pppbayestree/gpstk/PCSmoother.hpp"
#include "pppbayestree/gpstk/OceanLoading.hpp"
#include "pppbayestree/gpstk/CodeSmoother.hpp"
#include "pppbayestree/gpstk/SimpleFilter.hpp"
#include "pppbayestree/gpstk/MWCSDetector.hpp"
#include "pppbayestree/gpstk/SatArcMarker.hpp"
#include "pppbayestree/gpstk/DCBDataReader.hpp"
#include "pppbayestree/gpstk/ComputeWindUp.hpp"
#include "pppbayestree/gpstk/Rinex3NavData.hpp"
#include "pppbayestree/gpstk/GNSSconstants.hpp"
#include "pppbayestree/gpstk/ComputeLinear.hpp"
#include "pppbayestree/gpstk/GPSWeekSecond.hpp"
#include "pppbayestree/gpstk/LICSDetector2.hpp"
#include "pppbayestree/gpstk/DataStructures.hpp"
#include "pppbayestree/gpstk/RinexObsStream.hpp"
#include "pppbayestree/gpstk/Rinex3ObsStream.hpp"
#include "pppbayestree/gpstk/Rinex3NavStream.hpp"
#include "pppbayestree/gpstk/ComputeTropModel.hpp"
#include "pppbayestree/gpstk/SP3EphemerisStore.hpp"
#include "pppbayestree/gpstk/ComputeSatPCenter.hpp"
#include "pppbayestree/gpstk/EclipsedSatFilter.hpp"
#include "pppbayestree/gpstk/GPSEphemerisStore.hpp"
#include "pppbayestree/gpstk/CorrectCodeBiases.hpp"
#include "pppbayestree/gpstk/ComputeSatPCenter.hpp"
#include "pppbayestree/gpstk/RequireObservables.hpp"
#include "pppbayestree/gpstk/CorrectObservables.hpp"
#include "pppbayestree/gpstk/LinearCombinations.hpp"
#include "pppbayestree/gpstk/GravitationalDelay.hpp"
#include "pppbayestree/gpstk/PhaseCodeAlignment.hpp"
#include "pppbayestree/slam/dataset.h"

// BOOST
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

// STD
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <Eigen/Core>

using namespace std;
using namespace gpstk;
using namespace boost;
using namespace std::chrono;
using namespace minisam;
namespace po = boost::program_options;


int main(int argc, char* argv[])
{
    // define std out print color
    vector<int> prn_vec;
    vector<rnxData> data;
    const string red("\033[0;31m");
    const string green("\033[0;32m");
    string confFile, gnssFile, station, sp3File, p1p2File, p1c1File, antennaModel;
    string rnx_file, nav_file, sp3_file, out_file, antexFile;
    double xn, yn, zn, xp, yp, range, phase, rho, minElev, weightFactor;
    double antennaOffSetH, antennaOffSetE, antennaOffSetN;
    int startKey(0), currKey, startEpoch(0), svn, doy;
    int nThreads(-1), phase_break, break_count(0), nextKey, dec_int, itsBelowThree=0, count=0;
    bool printECEF, printENU, printAmb, printUpdateRate, first_ob(true), usingP1(false);

    FILE *fprealtime=fopen("examples_tuning/gpsdata/pppgpsposbackcount.txt","w+");

    cout.precision(12);

    /*po::options_description desc("Available options");
    desc.add_options()
            ("help,h", "Print help message")
            ("confFile,c", po::value<string>(&confFile)->default_value(""),
            "Input config file" )
            ("out", po::value<string>(&out_file)->default_value(""),
            "output file.")
            ("usingP1", "Are you using P1 instead of C1?");
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);*/

    ConfDataReader confReader;
    confReader.open("examples_tuning/gpsdata/phastball.conf");
    ISAM2Data isam2data;

    /*if (confFile.empty() ) {
            cout << red << "\n\n Currently, you need to provide a conf file \n"
                 << "\n\n"  << green << desc << endl;
    }*/

    while ( (station = confReader.getEachSection()) != "" )
    {
        // Fetch nominal station location [m]
        xn = confReader.fetchListValueAsDouble("nominalECEF",station);
        yn = confReader.fetchListValueAsDouble("nominalECEF",station);
        zn = confReader.fetchListValueAsDouble("nominalECEF",station);
        // day of year ( used for Niell Trop model)
        doy = confReader.getValueAsInt("DOY", station);
        // Elevation cut-off
        minElev = confReader.getValueAsDouble("minElev", station);
        // Code/carrier ratio
        weightFactor = confReader.getValueAsDouble("weightFactor", station);
        // Data file names
        gnssFile = confReader("rnxFile", station);
        sp3File = confReader("sp3File", station);
        p1p2File = confReader("p1p2", station);
        p1c1File = confReader("p1c1", station);
        // Print statements
        printENU = confReader.getValueAsBoolean("printENU", station);
        printAmb = confReader.getValueAsBoolean("printAmb", station);
        printECEF = confReader.getValueAsBoolean("printECEF", station);
        printUpdateRate = confReader.getValueAsBoolean("printUpdateRate", station);
    }

    //usingP1 = (vm.count("usingP1")>0);

    minivector nomXYZ(xn, yn, zn);
    minivector prop_xyz = nomXYZ;


    /*#ifdef USE_TBB
    std::auto_ptr<tbb::task_scheduler_init> init;
    if(nThreads > 0) {
            init.reset(new tbb::task_scheduler_init(nThreads));
    }
    else
            cout << green << " \n\n Using threads for all processors" << endl;
    #else
    */
    if(nThreads > 0)
    {
        cout << red <<" \n\n not compiled with TBB, so threading is"
             << " disabled and the --threads option cannot be used."
             << endl;
        exit(1);
    }
    //#endif

    //ISAM2DoglegParams doglegParams;
    ISAM2Params parameters;
    //parameters.optimizationParamsDogleg=new ISAM2DoglegParams;
    parameters.optimizationParamsGaussNewton=new ISAM2GaussNewtonParams;
    parameters.relinearizeThresholdDouble = 0.1;
    parameters.relinearizeSkip = 10;
    ISAM2 isam(parameters);

    double output_time = 0.0;
    double rw = 2.5;
    double rangeWeight = pow(rw,2);
    double phaseWeight = pow(rw*1/weightFactor,2);

    string value;

    minivector prior_nonBias(5,0.0);
    //prior_nonBias<<0.0, 0.0, 0.0, 0.0, 0.0;

    double bias_state(0.0);
    minivector phase_arc(34,0.0);
    //phase_arc.setZero();
    minivector bias_counter(34,0.0);
   // bias_counter.setZero();

    for (int i=1; i<34; i++)
    {
        bias_counter.data[i] = bias_counter.data[i-1] + 10000;
    }

    minivector* initEst5=new minivector(5);
    minivector_set_zero(initEst5);
   // initEst5<<0.0, 0.0, 0.0, 0.0, 0.0;

    minivector*  between_nonBias_State=new minivector(5);
    minivector_set_zero(between_nonBias_State);
    //between_nonBias_State<<0.0, 0.0, 0.0, 0.0, 0.0;

    //Values initial_values;
    std::map<int,minimatrix*> initial_values;
   // std::map<int,Pose3> initial_valuesP;
    //Values result;
   // std::map<int,Eigen::VectorXd> result;

    minivector BiasinitNoise(5);
    BiasinitNoise.data[0]=10.0;BiasinitNoise.data[1]=10.0;BiasinitNoise.data[2]=10.0;
    BiasinitNoise.data[3]=3e8;BiasinitNoise.data[4]=1e-1;
    //BiasinitNoise<<10.0, 10.0, 10.0, 3e8, 1e-1;
    minivector_cwisesqrt(&BiasinitNoise);
   // GaussianNoiseModel* nonBias_InitNoise=new GaussianNoiseModel(BiasinitNoise.cwiseSqrt());
    GaussianNoiseModel* nonBias_InitNoise=new GaussianNoiseModel(BiasinitNoise);

    minivector nonBiasProcessNoise(5);
    nonBiasProcessNoise.data[0]=0.1;nonBiasProcessNoise.data[1]=0.1;nonBiasProcessNoise.data[2]=0.1;
    nonBiasProcessNoise.data[3]=3e6;nonBiasProcessNoise.data[4]=3e-5;
    minivector_cwisesqrt(&nonBiasProcessNoise);

     GaussianNoiseModel* nonBias_ProcessNoise=new GaussianNoiseModel(nonBiasProcessNoise);



    //0.1, 0.1, 0.1, 3e6, 3e-5;
    minivector init_Noise(1,100);//init_Noise[0]=100;
    //init_Noise<<100;
    minivector_cwisesqrt(&init_Noise);

     GaussianNoiseModel* initNoise=new GaussianNoiseModel(init_Noise);

    string obs_path = findExampleDataFile(gnssFile,"examples_tuning/gpsdata");
    string p1p2_path = findExampleDataFile(p1p2File,"examples_tuning/gpsdata");
    string p1c1_path = findExampleDataFile(p1c1File,"examples_tuning/gpsdata");



    if ( obs_path.empty() )
    {
        // cout << " Must pass in obs file !!! " << desc << endl;
        exit(1);
    }
    if ( sp3File.empty() )
    {
        //cout << " Must pass in ephemeris file !!! " << desc << endl;
        exit(1);
    }
    // Declare a "SP3EphemerisStore" object to handle precise ephemeris
    SP3EphemerisStore SP3EphList;

    size_t pos = 0;
    string path, token;
    string delimiter = " ";
    while ((pos = sp3File.find(delimiter)) != string::npos)
    {
        path = findExampleDataFile(sp3File.substr(0,pos),"examples_tuning/gpsdata");
        SP3EphList.loadFile(path);
        sp3File.erase(0, pos + delimiter.length());
    }
    path = findExampleDataFile(sp3File,"examples_tuning/gpsdata");
    SP3EphList.loadFile(path);

    // Set flags to reject satellites with bad or absent positional
    // values or clocks
    SP3EphList.rejectBadPositions(true);
    SP3EphList.rejectBadClocks(true);

    // Create the input observation file stream
    Rinex3ObsStream rin(obs_path);

    // station nominal position
    Position nominalPos(xn, yn, zn);

    CorrectCodeBiases corrCode;
    corrCode.setDCBFile(p1p2_path, p1c1_path);

    if (!usingP1)
    {
        corrCode.setUsingC1(true);
    }


    // This is the GNSS data structure that will hold all the
    // GNSS-related information
    gnssRinex gRin;

    RequireObservables requireObs;
    requireObs.addRequiredType(TypeID::L1);
    requireObs.addRequiredType(TypeID::L2);

    SimpleFilter pObsFilter;
    pObsFilter.setFilteredType(TypeID::C1);

    if ( usingP1 )
    {
        requireObs.addRequiredType(TypeID::P1);
        pObsFilter.addFilteredType(TypeID::P1);
        requireObs.addRequiredType(TypeID::P2);
        pObsFilter.addFilteredType(TypeID::P2);
    }
    else
    {
        requireObs.addRequiredType(TypeID::C1);
        pObsFilter.addFilteredType(TypeID::C1);
        requireObs.addRequiredType(TypeID::P2);
        pObsFilter.addFilteredType(TypeID::P2);
    }

    // Declare a couple of basic modelers
    BasicModel basic(nominalPos, SP3EphList);
    basic.setMinElev(minElev);

    // Object to correct for SP3 Sat Phase-center offset
    ComputeSatPCenter svPcenter(SP3EphList, nominalPos);

    // Objects to mark cycle slips
    MWCSDetector markCSMW;  // Checks Merbourne-Wubbena cycle slip

    // object def several linear combinations
    LinearCombinations comb;

    // Object to compute linear combinations for cycle slip detection
    ComputeLinear linear1;
    if ( usingP1 )
    {
        linear1.addLinear(comb.pdeltaCombination);
        linear1.addLinear(comb.mwubbenaCombination);
    }
    else
    {
        linear1.addLinear(comb.pdeltaCombWithC1);
        linear1.addLinear(comb.mwubbenaCombWithC1);
    }
    linear1.addLinear(comb.ldeltaCombination);
    linear1.addLinear(comb.liCombination);

    ComputeLinear linear2;

    // Read if we should use C1 instead of P1
    if ( usingP1 )
    {
        linear2.addLinear(comb.pcCombination);
    }
    else
    {
        linear2.addLinear(comb.pcCombWithC1);
    }
    linear2.addLinear(comb.lcCombination);

    LICSDetector2 markCSLI2;       // Checks LI cycle slips

    // Object to keep track of satellite arcs
    SatArcMarker markArc;
    markArc.setDeleteUnstableSats(true);

    // Objects to compute gravitational delay effects
    GravitationalDelay grDelay(nominalPos);

    // Object to remove eclipsed satellites
    EclipsedSatFilter eclipsedSV;

    //Object to compute wind-up effect
    ComputeWindUp windup( SP3EphList, nominalPos );

    // Object to compute prefit-residuals
    ComputeLinear linear3(comb.pcPrefit);
    linear3.addLinear(comb.lcPrefit);

    TypeIDSet tset;
    tset.insert(TypeID::prefitC);
    tset.insert(TypeID::prefitL);

    // Declare a NeillTropModel object, setting the defaults
    NeillTropModel neillTM( nominalPos.getAltitude(),
                            nominalPos.getGeodeticLatitude(),
                            doy);

    // Objects to compute the tropospheric data
    ComputeTropModel computeTropo(neillTM);

    // initialize factor graph
    //NonlinearFactorGraph *graph = new NonlinearFactorGraph();
    NonlinearFactorGraph graph;

    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    // Loop over all data epochs
    while(rin >> gRin)
    {
        TimeSystem sys;
        sys.fromString("GPS");
        CommonTime time(gRin.header.epoch);
        time.setTimeSystem(sys);
        GPSWeekSecond gpstime( time );

        // update nominal ECEF with propogated pos.
        NeillTropModel neillTM( nominalPos.getAltitude(),
                                nominalPos.getGeodeticLatitude(),
                                doy);
        try
        {
            gRin >> requireObs // Check if required observations are present
                 >> pObsFilter // Filter out spurious data
                 >> linear1 // Compute linear combinations to detect CS
                 >> markCSLI2 // Mark cycle slips
                 >> markArc // Keep track of satellite arcs
                 >> basic // Compute the basic components of model
                 >> eclipsedSV // Remove satellites in eclipse
                 >> grDelay // Compute gravitational delay
                 >> svPcenter // Computer delta for sat. phase center
                 >> corrCode // Correct for differential code biases
                 >> windup // phase windup correction
                 >> computeTropo // neill trop function
                 >> linear2  // Compute ionosphere-free combinations
                 >> linear3;   // Compute prefit residuals
        }
        catch(Exception& e)
        {
            continue;
        }
        catch(...)
        {
            cerr << "Unknown exception at epoch: " << time << endl;
            continue;
        }

        // Iterate through the GNSS Data Structure
        satTypeValueMap::const_iterator it;
        typeValueMap::const_iterator itObs;
        if ( itsBelowThree > 0 )
        {
            itsBelowThree = 0;
            continue;
        }

        if (gRin.body.size() == 0)
        {
            continue;
        }

        // Loop over all observed sats at current epoch
        for (it = gRin.body.begin(); it!= gRin.body.end(); it++)
        {

            start = std::chrono::steady_clock::now();
            svn = ((*it).first).id;
            double satX, satY, satZ;
            satX = (*it).second.getValue(TypeID::satX);
            satY = (*it).second.getValue(TypeID::satY);
            satZ = (*it).second.getValue(TypeID::satZ);
            minivector satXYZ(satX,satY,satZ);
            double range, rangeRes;
            range = (*it).second.getValue(TypeID::PC);
            rangeRes = (*it).second.getValue(TypeID::prefitC);
            double phase, phaseRes;
            phase = (*it).second.getValue(TypeID::LC);
            phaseRes = (*it).second.getValue(TypeID::prefitL);
            int phase_break;
            phase_break = (*it).second.getValue(TypeID::satArc);

            if (first_ob)
            {
                startKey = count;
                first_ob=false;
                PriorFactor* npn=new PriorFactor(Symbol('X',count).key(), initEst5,  nonBias_InitNoise);

                graph.push_back(npn);
                initial_values.insert(std::make_pair(Symbol('X',count).key(), new minivector(initEst5)));

            }

            if (phase_arc.data[svn]!=phase_break)
            {
                bias_state = phase - range;
                if (count > startKey)
                {
                    bias_counter.data[svn] = bias_counter.data[svn] +1;
                }
                minivector *biasb=new minivector(1);
                biasb->data[0]=bias_state;
               // biasb<<bias_state;

                //GaussianNoiseModel* initNoise=new GaussianNoiseModel(init_Noise.cwiseSqrt());


                PriorFactor* nphb=new PriorFactor(Symbol('B',bias_counter.data[svn]).key(),
                biasb,initNoise);
                graph.push_back(nphb);
                initial_values.insert(std::make_pair(Symbol('B',bias_counter.data[svn]).key(), new minivector(biasb)));
                phase_arc.data[svn] = phase_break;
            }
            // Generate pseudorange factor
            minivector gpsRangeFactorvec(1,sqrt(elDepWeight(satXYZ, nomXYZ, rangeWeight)));

            //gpsRangeFactorvec<<elDepWeight(satXYZ, nomXYZ, rangeWeight);
            //GaussianNoiseModel* ngpsrfn=new GaussianNoiseModel(gpsRangeFactorvec.cwiseSqrt());
            GaussianNoiseModel* ngpsrfn=new GaussianNoiseModel(gpsRangeFactorvec);
            PseudorangeFactor* ngpsrf=new PseudorangeFactor(Symbol('X',count).key(), rangeRes,
             new minivector(satXYZ), new minivector(nomXYZ),ngpsrfn);

            graph.push_back(ngpsrf);
            gpsRangeFactorvec.data[0]=sqrt(elDepWeight(satXYZ, nomXYZ, phaseWeight));
            GaussianNoiseModel* ngpspfn2=new GaussianNoiseModel(gpsRangeFactorvec);
            //GaussianNoiseModel* ngpspfn2=new GaussianNoiseModel(gpsRangeFactorvec.cwiseSqrt());
            PhaseFactor* ngpspf=new PhaseFactor(Symbol('X',count).key(),
             Symbol('B',bias_counter.data[svn]).key(),phaseRes, new minivector(satXYZ),
             new minivector( nomXYZ),ngpspfn2);
            graph.push_back(ngpspf);

            prn_vec.push_back(svn);
        }
        if (count > startKey )
        {
           // GaussianNoiseModel* nonBias_ProcessNoise=
           // new GaussianNoiseModel(nonBiasProcessNoise.cwiseSqrt());

            BetweenFactor* nbfnbs=new BetweenFactor(Symbol('X',count).key(),
            Symbol('X',count-1).key(), new minivector(between_nonBias_State), nonBias_ProcessNoise);
            graph.push_back(nbfnbs);
        }

        isam.update(graph, initial_values,isam2data);
        //result = isam.calculateEstimate(isam2data);
        isam.calculateEstimate(isam2data);

        end = std::chrono::steady_clock::now();


        //prior_nonBias =result.at(Symbol('X',count).key());
        //prior_nonBias =isam2data.result.at(Symbol('X',count).key());
        minimatrix_memcpy(&prior_nonBias,isam2data.resulttheta_.at(Symbol('X',count).key()));

        minivector delta_xyz(prior_nonBias.data[0], prior_nonBias.data[1], prior_nonBias.data[2]);
        Position deltaPos(prior_nonBias.data[0], prior_nonBias.data[1], prior_nonBias.data[2]);
        //prop_xyz = nomXYZ - delta_xyz;
        minivector_sub(&prop_xyz,nomXYZ,delta_xyz);
        nominalPos -= deltaPos;

        if (printECEF)
        {
            cout << "xyz " << gpstime.week << " " << gpstime.sow << " " << prop_xyz.x() <<
            " " << prop_xyz.y() << " " << prop_xyz.data[2] << endl;
        }

        if (printENU)
        {
            minivector enu = xyz2enu(prop_xyz, nomXYZ);
            cout << "enu " << gpstime.week << " " << gpstime.sow << " " <<
            enu.x() << " " << enu.y() << " " << enu.data[2] << endl;
            fprintf(fprealtime,"%d %.15f %.15f %.15f %.15f %d\n",
                    gpstime.week,gpstime.sow,enu.x(),enu.y(),enu.data[2],
                    isam.lastBacksubVariableCount);
        }

        if (printAmb)
        {
            for (int k=0; k<prn_vec.size(); k++)
            {
                cout << "amb. " << gpstime.week << " " << gpstime.sow << " ";
                cout << prn_vec[k] << " ";
                //cout << result.at(Symbol('B',bias_counter[prn_vec[k]]).key()) << endl;
                minimatrix_print(isam2data.resulttheta_.at(Symbol('B',
                bias_counter.data[prn_vec[k]]).key()));
                cout << endl;
            }
        }

        if (printUpdateRate)
        {
            cout << "Elapsed time "
                 << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
                 << " µs" << endl;
        }

        output_time = output_time +1;
        graph.resize(0);
        initial_values.clear();
        prn_vec.clear();
        count++;
        initial_values.insert(std::make_pair(Symbol('X',count).key(),
        new minivector(prior_nonBias)));
       if(count>20000)
        {
         isam.clearall();
    isam2data.clearvalues();
    isam2data.clearfactors();
    delete parameters.optimizationParamsGaussNewton;
    delete between_nonBias_State;
    delete nonBias_InitNoise;
    delete initNoise;
     delete nonBias_ProcessNoise;
     for(auto& dlt:initial_values)
    {
       if(dlt.second!=NULL)
       delete dlt.second;
    }
    return 0;
        }
    }
   // ofstream bs("examples_tuning/data/gnsstree.dot");
    //isam.saveGraph(bs);

    isam.clearall();
    isam2data.clearvalues();
    isam2data.clearfactors();
    delete parameters.optimizationParamsGaussNewton;
    delete between_nonBias_State;
    delete nonBias_InitNoise;
    delete initNoise;
    delete nonBias_ProcessNoise;
    for(auto& dlt:initial_values)
    {
       if(dlt.second!=NULL)
       delete dlt.second;
    }

    return 0;
}
/* ************************************************************************* */


