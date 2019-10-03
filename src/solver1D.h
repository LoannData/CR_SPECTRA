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
                          // Pressure advection term
                          V1   = - VA[i][e]*(dt/dX)*0.5*(Pcr[i+1][e] - Pcr[i-1][e]);  
                          // Pressure diffusion terms    
                          V2A  = 0.5*((Db[i][e]/Ip[i][e]) + (Db[i+1][e]/Ip[i+1][e]))*(dt/pow(dX,2))*(Pcr[i+1][e] - Pcr[i][e]);
                          V2B  = -0.5*((Db[i][e]/Ip[i][e]) + (Db[i-1][e]/Ip[i-1][e]))*(dt/pow(dX,2))*(Pcr[i][e] - Pcr[i-1][e]);
                          // Adiabatic losses
                          V3   = E/3.*(dt/(dE*dX))*(0.25*(VA[i+1][e] - VA[i-1][e])*(Pcr[i][e+1] - Pcr[i][e-1]) - 0.25*(VA[i][e+1] - VA[i][e-1])*(Pcr[i+1][e] - Pcr[i-1][e]));
                          // Alfven velocity fluctuation in space 
                          V4   = 0.5*4./3*(dt/dX)*Pcr[i][e]*(VA[i+1][e] - VA[i-1][e]);
                          // Full term
                          double Pcr_new = Pcr[e][i] + solver_PcrDiffusion*(V2A + V2B) + solver_energy*V3 + V4 + Q; 

                          //Pcr_new = Pcr_new + solver_PcrAdvection*V1;
                          return Pcr_new;
                      }

double advectPressure(double dt, double dX, double dE, double Pcr[5][5], double Ip[5][5], double VA[5][5], double Gd[5][5], double Db[5][5], double E,
                      double Q)
                      {
                          /*double V1;
                          int i = 2; int e = 2; 
                          // Pressure advection term
                          V1   = - VA[i][e]*(dt/dX)*0.5*(Pcr[i+1][e] - Pcr[i-1][e]);  
                          double Pcr_new = Pcr[e][i] + solver_PcrAdvection*V1;*/
                          
                          int i = 2; int e = 2;
                          double Pcr_new = pow(1 + VA[i][e]*dt/dX, -1)*(Pcr[e][i] + VA[i][e]*dt/dX*Pcr[e][i-1])*solver_PcrAdvection;

                          return Pcr_new;
                      }

double wavesSolver(double dt, double dX, double dE, double Pcr[5][5], double Ip[5][5], double VA[5][5], double Gd[5][5], double Db[5][5], double B,
                   double Q)
                      {
                          double V1, V2, V3, V4, V5; 
                          int i = 2; int e = 2;
                          // Waves advection term
                          V1 = - dt/dX*VA[i][e]*0.5*(Ip[i+1][e] - Ip[i-1][e]);
                          // Alfven velocity fluctuation in space 
                          V2 = - dt/dX*Ip[i][e]*0.5*(VA[i+1][e] - VA[i-1][e]);
                          // Waves Growth term
                          double ff = - 0.5*VA[i][e]/dX*0.5*(Pcr[i+1][e] - Pcr[i-1][e])/(pow(B,2)/(8*pi)); 
                          //V3 = ff*dt; // Terme de taux de croissance linéaire 
                          V3 = ff*dt*(log10(Ip[e][i] + 1)/Ip[e][i]); // Terme de taux de croissance avec une fin logarithmique
                          // Waves Damping term
                          V4 = - Gd[i][e]*Ip[i][e]*dt;
                          // Waves Source term
                          V5 = Q*dt;
                          //Full term
                          double Ip_new = Ip[e][i]  + V2 + solver_Ipgrowth*V3 + solver_damping*V4 + solver_waves_source*V5;

                          //Ip_new  = Ip_new + solver_IpAdvection*V1;
                          return Ip_new;
                      }

double advectWaves(double dt, double dX, double dE, double Pcr[5][5], double Ip[5][5], double VA[5][5], double Gd[5][5], double Db[5][5], double B,
                   double Q)
                      {
                          /*double V1; 
                          int i = 2; int e = 2;
                          // Waves advection term
                          V1 = - dt/dX*VA[i][e]*0.5*(Ip[i+1][e] - Ip[i-1][e]);

                          double Ip_new  = Ip[i][e] + solver_IpAdvection*V1;*/

                          int i = 2; int e = 2;
                          double Ip_new = pow(1 + VA[i][e]*dt/dX, -1)*(Ip[e][i] + VA[i][e]*dt/dX*Ip[e][i-1])*solver_IpAdvection;

                          return Ip_new;
                      }





