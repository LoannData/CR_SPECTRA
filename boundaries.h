#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;




void boundaries(int xi, int ei, int NX, int NE, double Pcr, double Ip, int xi_t, int ei_t, double Pcrt[5][5], double Ipt[5][5])
{
    if (xi == 0 && (ei != 0 && ei != NE-1))
    {
        Pcr = Pcrt[xi_t+1][ei_t]; 
        Ip  = Ipt[xi_t+1][ei_t];
    }

    if (xi == NX-1 && (ei !=0 && ei != NE-1))
    {
        Pcr = Pcrt[xi_t-1][ei_t];
        Ip  = Ipt[xi_t-1][ei_t];

    }

    if (ei == 0 && (xi != 0 && xi != NX-1))
    {
        Pcr = Pcrt[xi_t][ei_t];
        Ip  = Ipt[xi_t][ei_t];

    }

    if (ei == NE-1 && (xi != 0 && xi != NX-1))
    {
        Pcr = Pcrt[xi_t][ei_t];
        Ip  = Ipt[xi_t][ei_t];

    }

    if (xi == 0 && ei == 0)
    {
        Pcr = Pcrt[xi_t+1][ei_t]; 
        Ip  = Ipt[xi_t+1][ei_t];
    }

    if (xi == NX-1 && ei == 0)
    {
        Pcr = Pcrt[xi_t-1][ei_t]; 
        Ip  = Ipt[xi_t-1][ei_t];
    }

    if (xi == 0 && ei == NE-1)
    {
        Pcr = Pcrt[xi_t+1][ei_t]; 
        Ip  = Ipt[xi_t+1][ei_t];
    }

    if (xi == NX-1 && ei == NE-1)
    {
        Pcr = Pcrt[xi_t-1][ei_t]; 
        Ip  = Ipt[xi_t-1][ei_t];
    }
}  






/*void iterate(vector<vector<vector<double>>>)
{

}*/