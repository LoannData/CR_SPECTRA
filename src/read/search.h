#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;


std::string search(string file_name, string variable) 
{
    ifstream file; 
    std::string value;
    int vect_size = 10000; 
    std::vector<string> values(vect_size);
    file.open(file_name);
    if (!file) cout<<"ERROR: parameter_file \n";
    int kk=0;
    while (!file.eof())
    {
        file >> values[kk]; 
        kk++;
    }
    file.close();
    
    for (int kk=0; kk<vect_size; kk++)
    {
        
        if (variable == values[kk])
        {
            value = values[kk+1];
        };
    };
    return value;
}