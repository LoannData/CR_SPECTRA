#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;

int writeXE(std::string filename, int index, vector<vector<double>> data, double NX, double NE)
{
    ofstream data_file;
    index += 1;
    std::string full_name = "./data_out/";
    std::string ext=".dat";
    full_name.append(filename); 
    full_name.append("_");
    std::string idcomplement;
    if (index < 10) {idcomplement = "0000";}
    if (index >= 10 && index < 100) {idcomplement = "000";}
    if (index >= 100 && index < 1000) {idcomplement = "00";}
    if (index >= 1000 && index < 10000) {idcomplement = "0";}
    if (index >= 10000 && index < 100000) {idcomplement = "";}

    std::ostringstream istring;
    istring<<index;
    full_name.append(idcomplement);
    full_name.append(istring.str());
    full_name.append(ext); 

    data_file.open(full_name.c_str());
    data_file<<",";
    for (int ei = 0; ei < NE; ei++){data_file<<ei<<",";}
    data_file<<"Index"<<endl;

    for (int xi = 0; xi < NX; xi++)
    {
        data_file<<xi<<",";
        for (int ei = 0; ei < NE; ei++)
        {
            data_file<<data[xi][ei]<<",";
        }
        data_file<<xi<<endl;
    }
    data_file.close();

    //index += 1;
    return index;
}  