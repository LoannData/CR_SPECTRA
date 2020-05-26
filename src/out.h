#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;



/// This function allows you to specify your own data output time array 
/// You just need to fill the list loc_data,
/// example: double loc_data[] = {t1, t2, ..., tn}; where t_i is in seconds 
vector<double> specificOutputData()
{
    vector<double> data; 
    double loc_data[] = {0.1*kyr, 0.5*kyr, 1*kyr, 5*kyr, 10*kyr, 50*kyr, 100*kyr, 200*kyr};
    for (int vi=0; vi < sizeof(loc_data)/sizeof(*loc_data); vi++)
    {
        data.push_back(loc_data[vi]);
    }
    return data;
}




/// Output data function for a regular and large amount of output data
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
        double log_tmin = log10(eps);
        //if (tmin == 0.) {log_tmin = log10(eps);}
        //else            {log_tmin = log10(tmin);}
        double loc_value;
        vector<double> log_data = outputData(log_tmin, log_tmax, nvalues, "linear", eps);
        for (int vi=0; vi < nvalues; vi++)
        {   
            loc_value = pow(10., log_data[vi]);
            data.push_back(loc_value);
        }
    }
    if (modulation == "custom")
    {
        data = specificOutputData();
    }
    return data;
}









//==========================================================================//
// OUTPUT PARAMETERS                                                        // 
//==========================================================================//


/// Main output file function. This one will be used in the main file
vector<double> getOutput()
{
    std::string distribution;
    if (time_distrib_of_data == 0){distribution = "linear";}
    if (time_distrib_of_data == 1){distribution = "log10";}
    if (time_distrib_of_data == 2){distribution = "custom";}

    // Choose the way you want to output your data 
    vector<double> data = outputData(t_data_out_min, t_data_out_max, number_out_data, distribution, log_first_data);
    return data; 
}

/// Main output log file function. This one will be used in the main file 
int getLogOutput()
{
    // Choose the interval of time step between two log outputs 
    int step = delta_log_output;
    return step;
}

/// Define the limit time of your simulation 
double setTmax()
{
    //double Tmax = 100.*kyr;
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
