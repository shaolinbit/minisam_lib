/**
 * @file   GnssData.cpp
 * @brief  Tools required to read/write GNSS data
 * @author Ryan Watson
 */

#include <iomanip>      // std::setprecision
#include "GnssData.h"
#include "../slam/dataset.h"

using namespace std;


vector<rnxData> readGNSS(const std::string &fileLoc,const std::string& filepath)
{
    /*
       inputs ::
       fileLoc ---> path to data file
       output ::
       data ---> gnss data in gtsam format
                       { epoch, svn, satXYZ, computed_range, rangeLC, phaseLC }
     */
    vector<rnxData> data;
    string data_file = findExampleDataFile(fileLoc,filepath);
    ifstream is(data_file.c_str());

    while (is)
    {
        int svn, count;
        double break_flag, grav_delay, windup, satPC, trop_slant, c1Del, c2Del;
        double week, sow, satX, satY, satZ, rho, cb, rel, rangeLC, phaseLC;
        //Eigen::Vector3d satXYZ, computed_range;
        string constellation;
        is >> week >> sow  >> count >> constellation
           >> svn >> rangeLC >> phaseLC
           >> rho >> cb >> rel >> grav_delay >> trop_slant >>  windup >> satPC >> satX >> satY >> satZ >> break_flag >> c1Del >> c2Del;
        data.push_back(rnxData(sow, count, svn,minivector(satX,satY,satZ),
                               (rho - cb  + rel  + grav_delay + trop_slant - satPC), (rangeLC - c1Del + c2Del), (phaseLC - windup*0.017), break_flag));
        // 0.01702215881 == LC wavelength/2*pi
    }
    is.clear();         /* clears the end-of-file and error flags */
    return data;
}

void writeStates(std::map<int,minivector> &results, string outputFile)
{
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
        minivector p = keypair.second;
        outFile << "stateVec " << epoch++
                << " "  << p.data[0] << " " << p.data[1]
                << " "  <<p.data[2]  << " " << p.data[3]
                << " " << p.data[4] << endl;
    }
}

void writeNavFrame(std::map<int,minivector> &results, minivector &nom, string outputFile)
{
    ofstream outFile(outputFile.c_str());
    int epoch = 0;
    foreach (const auto& keypair, results)
    {
        minivector p = keypair.second;
        //Eigen::Vector3d delta(p.x(),p.y(),p.z());
        // Eigen::Vector3d delta(p(0),p(1),p(2));
        minivector delta=minivector_subvector(p,0,3);

        minivector ecef=minivector(nom.data[0] - delta.data[0],nom.data[1] - delta.data[1],
                                        nom.data[2] - delta.data[2]);
        minivector enu = xyz2enu(ecef,nom);
        outFile << epoch++ << " " << enu.data[0]
                << " " <<enu.data[1] << " " << enu.data[2]<< endl;

    }

}

void writeEarthFrame(std::map<int,minivector> &results, minivector &nom, string outputFile)
{
    ofstream outFile(outputFile.c_str());
    int epoch = 0;
    foreach (const auto& keypair, results)
    {
        minivector p = keypair.second;
        minivector delta=minivector_subvector(p,0,3);
        minivector ecef(nom.data[0] - delta.data[0],nom.data[1] - delta.data[1],
                                        nom.data[2] - delta.data[2]);
        outFile << epoch++ << " " << std::setprecision(10) << ecef.data[0]
                << " " << ecef.data[1] << " " << ecef.data[2] << endl;
    }
}

void writeSwitches( std::map<int,minivector> &results, string outputFile, vector<string> switchIndex)
{
    /*
       inputs ::
       results --> optimizer output
       outputFile --> name of file to write switch states to [string]
       Optional ::
       switchIndex --> index by epoch and visible satellite (i.e. obs 4 would be Switch_0_4) [vector]
     */
    ofstream outFile(outputFile.c_str());
    int epoch = 0;
    foreach (const auto& keypair, results)
    {
        int index = epoch++;
        outFile << switchIndex[index] << " "
                << index << " " <<endl;
        minivector_print(keypair.second);
        outFile << endl;
    }


}

void writeAmbiguity(std::map<int,minivector> &results, string outputFile, vector<string> satIndex)
{
    ofstream outFile(outputFile.c_str());
    int epoch = 0;
    foreach (const auto& keypair, results)
    {
        int index = epoch++;
        //phaseBias p = key_value.value;
        outFile << satIndex[index] <<  " " ;
        //minivector_ofstream(outFile,keypair.second);
        for (int i = 0; i < keypair.second.size1; i++)
    {
        outFile<<keypair.second.data[i*keypair.second.prd]<<std::endl;
    }
        outFile<< endl;

    }
}
