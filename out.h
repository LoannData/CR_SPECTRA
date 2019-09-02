#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;

// Output data function for a regular and large amount of output data
vector<double> outputData(double tmin, double tmax, int nvalues, string modulation, double eps)  
{
    vector<double> data;
    if (modulation == "linear") 
    {
        double DT = tmax - tmin;
        double dt = DT/nvalues;
        double loc_value;
        for (int vi=0; vi < nvalues; vi++)
        {
            loc_value = tmin + (vi+1)*dt;
            data.push_back(loc_value);
        }
    }
    if (modulation == "log10")
    {
        double log_tmax = log10(tmax);
        double log_tmin;
        if (tmin == 0.) {log_tmin = log10(eps);}
        else            {log_tmin = log10(tmin);}
        double loc_value;
        vector<double> log_data = outputData(log_tmin, log_tmax, nvalues, "linear", eps);
        for (int vi=0; vi < nvalues; vi++)
        {   
            loc_value = pow(10., log_data[vi]);
            data.push_back(loc_value);
        }
    }
    return data;
}


vector<double> specificOutputData()
{
    vector<double> data; 
    double loc_data[] = {1, 2, 3, 4, 5, 6};
    for (int vi=0; vi < sizeof(loc_data)/sizeof(*loc_data); vi++)
    {
        data.push_back(loc_data[vi]);
    }
    return data;
}






//==========================================================================//
// OUTPUT PARAMETERS                                                        // 
//==========================================================================//


// Main output file function. This one will be used in the main file
vector<double> getOutput()
{
    // Choose the way you want to output your data 
    vector<double> data = outputData(0.*kyr, 100.*kyr, 1000, "linear", 0.001);
    return data; 
}

// Main output log file function. This one will be used in the main file 
int getLogOutput()
{
    // Choose the interval of time step between two log outputs 
    int step = 1000;
    return step;
}

// Define the limit time of your simulation 
double setTmax()
{
    double Tmax = 100.*kyr;
    return Tmax;
}





// Exemple output vector 
/*int main()
{
    //vector<double> data = outputData(0., 10., 10, "linear", 0.01);
    vector<double> data = specificOutputData();

    for (int line=0; line<data.size(); line++)
    {
        cout<<data[line]<<"\n";
    }
    return 0;
}*/
