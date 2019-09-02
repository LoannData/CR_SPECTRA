#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;



double pressureSolver(double dt, double dX, double dE, double Pcr[5][5], double Ip[5][5], double VA[5][5], double Gd[5][5], double Db[5][5], double E,
                      double Q)
                      {
                          double V1, V3, V4, V2A, V2B;
                          int i = 2; int e = 2; 

                          V1 = - VA[i][e]*(dt/dX)*0.5*(Pcr[i+1][e] - Pcr[i-1][e]);     
                          V2A = 0.5*((Db[i][e]/Ip[i][e]) + (Db[i+1][e]/Ip[i+1][e]))*(dt/pow(dX,2))*(Pcr[i+1][e] - Pcr[i][e]);
                          V2B = -0.5*((Db[i][e]/Ip[i][e]) + (Db[i-1][e]/Ip[i-1][e]))*(dt/pow(dX,2))*(Pcr[i][e] - Pcr[i-1][e]);
                          V3 = E/3.*(dt/(dE*dX))*(0.25*(VA[i+1][e] - VA[i-1][e])*(Pcr[i][e+1] - Pcr[i][e-1]) - 0.25*(VA[i][e+1] - VA[i][e-1])*(Pcr[i+1][e] - Pcr[i-1][e]));
                          V4 = 0.5*4./3*(dt/dX)*Pcr[i][e]*(VA[i+1][e] - VA[i-1][e]);

                          // Full term
                          double Pcr_new = Pcr[e][i] + V1 + V2A + V2B + V3 + V4 + Q; 
                          // No energy dependance 
                          //double Pcr_new = Pcr[e][i] + V1 + V2A + V2B + V4 + Q; 

                          return Pcr_new;
                      }

double wavesSolver(double dt, double dX, double dE, double Pcr[5][5], double Ip[5][5], double VA[5][5], double Gd[5][5], double Db[5][5], double B,
                   double Q)
                      {
                          double V1, V2, V3, V4, V5; 
                          int i = 2; int e = 2;

                          V1 = - dt/dX*VA[i][e]*0.5*(Ip[i+1][e] - Ip[i-1][e]);
                          V2 = - dt/dX*Ip[i][e]*0.5*(VA[i+1][e] - VA[i-1][e]);
                          double ff = - 0.5*VA[i][e]/dX*0.5*(Pcr[i+1][e] - Pcr[i-1][e])/(pow(B,2)/(8*pi)); 
                          V3 = ff*dt;
                          V4 = - Gd[i][e]*Ip[i][e]*dt;
                          V5 = Q*dt;
                          
                          //Full term
                          double Ip_new = Ip[e][i]  + V1 + V2 + V3 + V4 + V5;

                          return Ip_new;
                      }





