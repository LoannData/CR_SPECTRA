#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;



double maxElement1D(vector<double> data)
{
    int size = data.size();
    double max_element = 0.; 
    for (int ii=0; ii<size; ii++)
    {
        if (data[ii] >= max_element) {max_element = data[ii];}
    }
    return max_element; 
}

double minElement1D(vector<double> data)
{
    int size = data.size();
    double min_element = pow(10,300); 
    for (int ii=0; ii<size; ii++)
    {
        if (data[ii] <= min_element) {min_element = data[ii];}
    }
    return min_element; 
}

double maxElement2D(vector<vector<double> > data)
{
    int size_1 = data.size();
    double max_element = 0.;
    for (int ii=0; ii<size_1; ii++)
    {
        int size_2 = data[ii].size();
        for (int jj=0; jj<size_2; jj++)
        {
            if (data[ii][jj] > max_element) {max_element = data[ii][jj];}
        }
    }
    return max_element;
}

double minElement2D(vector<vector<double> > data)
{
    int size_1 = data.size();
    double min_element = pow(10,300);
    for (int ii=0; ii<size_1; ii++)
    {
        int size_2 = data[ii].size();
        for (int jj=0; jj<size_2; jj++)
        {
            if (data[ii][jj] < min_element) {min_element = data[ii][jj];}
        }
    }
    return min_element;
}