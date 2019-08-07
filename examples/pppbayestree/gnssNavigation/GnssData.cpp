/**
 * @file   GnssData.cpp
 * @brief  Tools required to read/write GNSS data
 * @author Ryan Watson
 */

#include <iomanip>      // std::setprecision
#include "GnssData.h"
#include "../slam/dataset.h"

using namespace std;


vector<rnxData> readGNSS(const std::string &fileLoc) {
        /*
           inputs ::
           fileLoc ---> path to data file
           output ::
           data ---> gnss data in gtsam format
                           { epoch, svn, satXYZ, computed_range, rangeLC, phaseLC }
         */
        vector<rnxData> data;
        string data_file = findExampleDataFile(fileLoc);
        ifstream is(data_file.c_str());

        while (is) {
                int svn, count;
                double break_flag, grav_delay, windup, satPC, trop_slant, c1Del, c2Del;
                double week, sow, satX, satY, satZ, rho, cb, rel, rangeLC, phaseLC;
                Eigen::Vector3d satXYZ, computed_range;
                string constellation;
                is >> week >> sow  >> count >> constellation
                >> svn >> rangeLC >> phaseLC
                >> rho >> cb >> rel >> grav_delay >> trop_slant >>  windup >> satPC >> satX >> satY >> satZ >> break_flag >> c1Del >> c2Del;
                data.push_back(rnxData(sow, count, svn,Eigen::Vector3d(satX,satY,satZ),
                                       (rho - cb  + rel  + grav_delay + trop_slant - satPC), (rangeLC - c1Del + c2Del), (phaseLC - windup*0.017), break_flag));
                // 0.01702215881 == LC wavelength/2*pi
        }
        is.clear();         /* clears the end-of-file and error flags */
        return data;
}

void writeStates(std::map<int,Eigen::VectorXd> &results, string outputFile){
        /*
           inputs ::
           results -->
           outputFile --> name of file to write state est. to. [string]
         */
        ofstream outFile(outputFile.c_str());
        int epoch = 0;
        //  Values::ConstFiltered<nonBiasStates> result_poses = results.filter<nonBiasStates>();
        foreach (const auto& keypair, results)
        {
                Eigen::VectorXd p = keypair.second;
                outFile << "stateVec " << epoch++
                        << " "  << p(0) << " " << p(1)
                        << " "  << p(2) << " " << p(3)
                        << " " << p(4) << endl;
        }
}

void writeNavFrame(std::map<int,Eigen::VectorXd> &results, Eigen::Vector3d &nom, string outputFile){
        ofstream outFile(outputFile.c_str());
        int epoch = 0;
         foreach (const auto& keypair, results)
        {
                Eigen::VectorXd p = keypair.second;
                //Eigen::Vector3d delta(p.x(),p.y(),p.z());
                Eigen::Vector3d delta(p(0),p(1),p(2));
                Eigen::Vector3d ecef = (nom - delta);
                Eigen::Vector3d enu = xyz2enu(ecef,nom);
                outFile << epoch++ << " " << enu(0)
                        << " " << enu(1) << " " << enu(2) << endl;

        }

}

void writeEarthFrame(std::map<int,Eigen::VectorXd> &results, Eigen::Vector3d &nom, string outputFile){
        ofstream outFile(outputFile.c_str());
        int epoch = 0;
        foreach (const auto& keypair, results)
        {
                Eigen::VectorXd p = keypair.second;
                Eigen::Vector3d delta(p(0),p(1),p(2));
                Eigen::Vector3d ecef = (nom - delta);
                outFile << epoch++ << " " << std::setprecision(10) << ecef.x()
                        << " " << ecef(1) << " " << ecef(2) << endl;

        }
}

void writeSwitches( std::map<int,Eigen::VectorXd> &results, string outputFile, vector<string> switchIndex){
        /*
           inputs ::
           results --> optimizer output
           outputFile --> name of file to write switch states to [string]
           Optional ::
           switchIndex --> index by epoch and visible satellite (i.e. obs 4 would be Switch_0_4) [vector]
         */
        ofstream outFile(outputFile.c_str());
        int epoch = 0;
        foreach (const auto& keypair, results) {
                int index = epoch++;
                outFile << switchIndex[index] << " "
                        << index << " " <<  keypair.second << endl;
        }


}

void writeAmbiguity(std::map<int,Eigen::VectorXd> &results, string outputFile, vector<string> satIndex){
        ofstream outFile(outputFile.c_str());
        int epoch = 0;
        foreach (const auto& keypair, results)
        {
                int index = epoch++;
                //phaseBias p = key_value.value;
                outFile << satIndex[index] <<  " " << keypair.second<< endl;

        }
}

