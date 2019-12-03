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





// New functions !!! 
void thetaDiffusionSolver(vector<vector<double>> &u, vector<vector<double>> &Pcr_new,   double dt, vector<double> X, int NE, vector<vector<double>> I, vector<vector<double>> Db, vector<vector<double>> &Pcr_background)
{
    //cout<<"Yo!"<<endl;
    double theta = 1.;
    //int nx = X.size();
    
    //vector<vector<double>> Pcr_new(nx);
    //for (int xi = 0; xi < nx; xi++)
    //{
    //    Pcr_new[xi].resize(NE);
    //}
    //cout << u.size()<<endl;

    for (int ei = 0; ei < NE; ei++)
    {
        int NX = X.size()-1; 
        vector<double> R(NX+1, 0.);
        vector<double> aa(NX+1, 0.);
        vector<double> bb(NX+1, 0.);
        vector<double> cc(NX+1, 0.);
        vector<double> loc_Pcr(NX+1,0.);



        double dx, F1, alpha_m, alpha_p, an, bn, cn, rn;
        for (int xi = 1; xi < NX; xi++)
        {

            dx = 0.5*(X[xi+1] - X[xi-1]);
            F1 = theta*dt/pow(dx, 2.);
            alpha_m = 0.5*(Db[xi][ei]/I[xi][ei] + Db[xi-1][ei]/I[xi-1][ei]); 
            alpha_p = 0.5*(Db[xi][ei]/I[xi][ei] + Db[xi+1][ei]/I[xi+1][ei]);

            //cout<<"ei = "<<ei<<", xi = "<<xi<<"dx = "<<dx<<", F1 = "<<F1<<", alpha_m = "<<alpha_m<<", alpha_p = "<<alpha_p<<"I = "<<I[ei][xi]<<"Db = "<<Db[ei][xi]<<endl;

            an = -F1*alpha_m;
            bn = F1*alpha_p + F1*alpha_m + 1;
            cn = -F1*alpha_p;
            rn = u[xi][ei] + (1-theta)*dt/pow(dx,2.)*(alpha_p*(u[xi+1][ei] - u[xi][ei]) - alpha_m*(u[xi][ei] - u[xi-1][ei]));

            aa[xi] = an;
            bb[xi] = bn;
            cc[xi] = cn;

            R[xi] = rn;

        }
        aa.erase(aa.begin()+0);
        cc.erase(cc.begin()+NX);

        // Cas xi = NX
        dx = X[1] - X[0]; 
        F1 = theta*dt/pow(dx,2);
        alpha_m = 0.5*(Db[NX][ei]/I[NX][ei] + Db[NX-1][ei]/I[NX-1][ei]);
        alpha_p = alpha_m;
        aa[NX] = -F1*alpha_m;
        bb[NX] = F1*alpha_p + F1*alpha_m +1; 

        // Cas xi = 0
        dx = X[1] - X[0];
        F1 = theta*dt/pow(dx,2);
        alpha_p = 0.5*(Db[0][ei]/I[0][ei] + Db[1][ei]/I[1][ei]);
        alpha_m = alpha_p;
        cc[0] = -F1*alpha_p;
        bb[0] = F1*alpha_p + F1*alpha_m +1; 

        //cout<<aa[0]<<"   "<<cc[NX-1]<<endl;

        loc_Pcr = TDMA(aa, bb, cc, R);

        loc_Pcr[0] = Pcr_background[0][ei];//loc_Pcr[1];
        loc_Pcr[NX-1] = Pcr_background[NX-1][ei];//loc_Pcr[NX-1];
        loc_Pcr[NX] = Pcr_background[NX][ei];

        

        for (int xi = 0; xi < NX+1; xi++){Pcr_new[xi][ei] = loc_Pcr[xi];}
    }

    //return Pcr_new;

}



void advectionSolverX(vector<vector<double>> &u_old, vector<vector<double>> &u_new, double dt, vector<double> X, int NE, vector<vector<double>> V)
{
    // Note : Problem with 2nd order scheme ! 
    int NX = X.size()-1;
    double ddx, ux_p, ux_m, a_p, a_m;
    int order = 1;
    for (int ei = 0; ei < NE; ei++)
    {   
        if (order == 1)
        {
            for (int xi = 1; xi < NX; xi++)
            {
                ddx = (X[xi+1] - X[xi-1])/2.;

                a_p = max(V[xi][ei], 0.);
                a_m = min(V[xi][ei], 0.);

                ux_p = (u_old[xi+1][ei] - u_old[xi][ei])/ddx;
                ux_m = (u_old[xi][ei] - u_old[xi-1][ei])/ddx;

                u_new[xi][ei] = u_old[xi][ei] - dt*(a_p*ux_m + a_m*ux_p);
            }
            // Cas xi = 0
            u_new[0][ei]  = u_old[0][ei];
            // Cas xi = NX
            u_new[NX][ei] = u_old[NX][ei];
        }

        if (order == 2)
        {
            for (int xi = 2; xi < NX-1; xi++)
            {
                ddx = (X[xi+1] - X[xi-1])/2.;

                a_p = max(V[xi][ei], 0.);
                a_m = min(V[xi][ei], 0.);

                ux_m = 0.5*(3*u_old[xi][ei] - 4*u_old[xi-1][ei] + u_old[xi-2][ei])/ddx;
                ux_p = 0.5*(-u_old[xi+2][ei] + 4*u_old[xi+1][ei] - 3*u_old[xi][ei])/ddx;

                u_new[xi][ei] = u_old[xi][ei] - dt*(a_p*ux_m + a_m*ux_p);
            }
            u_new[0][ei] = u_old[0][ei];
            u_new[1][ei] = u_old[1][ei];

            u_new[NX][ei] = u_old[NX][ei];
            u_new[NX-1][ei] = u_old[NX-1][ei];

        }

    }

}

void advectionSolverE(vector<vector<double>> &u_old, vector<vector<double>> &u_new, double dt, vector<double> E, int NX, vector<vector<double>> V, vector<vector<double>> u_background)
{
    // Note : No 2nd order scheme ! 
    // Note : For the boundaries, we use absorbing ones. If energy leaves from the borders, it disapear of ever ...
    // We need to find a good method for the boundaries ...  
    int NE = E.size()-1;
    double dde, ux_p, ux_m, a_p, a_m, ad0;
    int order = 1;
    for (int xi = 0; xi < NX; xi++)
    {
        if (order == 1)
        {
            for (int ei = 1; ei < NE; ei++)
            {
                dde = (E[ei+1] - E[ei-1])/2.;



                a_p = max(V[xi][ei], 0.);
                a_m = min(V[xi][ei], 0.);

                ux_p = (u_old[xi][ei+1] - u_old[xi][ei])/dde;
                ux_m = (u_old[xi][ei] - u_old[xi][ei-1])/dde;

                //cout<<a_p<<" "<<a_m<<" "<<ux_p<<" "<<ux_m<<endl;

                u_new[xi][ei] = u_old[xi][ei] - dt*(a_p*ux_m + a_m*ux_p);
                if (ei == 1){ad0 = - dt*(a_p*ux_m + a_m*ux_p);}
                //cout<<dt<<"  "<<a_p<<" "<<a_m<<" "<<ux_p<<" "<<ux_m<<" "<<- dt*(a_p*ux_m + a_m*ux_p)<<endl;
                if (u_new[xi][ei] < u_background[xi][ei]){u_new[xi][ei] = u_background[xi][ei];}
            }

            // Cas ei = 0
            u_new[xi][0]  = u_old[xi][0] + ad0;
            // Cas ei = NE
            u_new[xi][NE] = u_old[xi][NE] - dt*(a_p*ux_m + a_m*ux_p);

        }
    }
}

void sourceSolver(vector<vector<double>> &u_old, vector<vector<double>> &u_new, double dt, vector<vector<double>> source, double factor)
    {
        int NX = u_old.size();
        int NE;
        double sum;
        for (int xi = 0; xi < NX; xi++)
        {
            NE = u_old[xi].size();
            for (int ei = 0; ei < NE; ei++)
            {
                sum = factor*source[xi][ei]*u_old[xi][ei];
                u_new[xi][ei] = u_old[xi][ei] + sum*dt;
            }
        }
    }


void sourceSolverDamp(vector<vector<double>> &u_old, vector<vector<double>> &u_new, double dt, vector<vector<double>> source, vector<vector<double>> background, double factor)
    {
        int NX = u_old.size();
        int NE;
        double sum;
        for (int xi = 0; xi < NX; xi++)
        {
            NE = u_old[xi].size();
            for (int ei = 0; ei < NE; ei++)
            {
                sum = factor*source[xi][ei]*(u_old[xi][ei] - background[xi][ei]);
                u_new[xi][ei] = u_old[xi][ei] + sum*dt;
            }
        }
    }


void sourceGrowthRateSolver(vector<vector<double>> &u_old, vector<vector<double>> &u_new, vector<vector<double>> &v_old, vector<double> X, double dt, vector<vector<double>> V, vector<double> B, double factor)
    {
        int NX = u_old.size();
        int NE;
        double sum;
        double dudx; 
        for (int xi = 0; xi < NX; xi++)
        {
            NE = u_old[xi].size();
            for (int ei = 0; ei < NE; ei++)
            {
                dudx = 0.;
                if (xi > 0 and xi < NX-1){dudx = (v_old[xi+1][ei]-v_old[xi-1][ei])/(X[xi+1] - X[xi-1]);}
                sum = -0.5*V[xi][ei]*dudx/(pow(B[xi],2)/(8*pi))*(log10(u_old[xi][ei]+1)/u_old[xi][ei]);
                //sum = factor*source[xi][ei]*u_old[xi][ei];
                u_new[xi][ei] = u_old[xi][ei] + sum*dt;
            }
        }
    }
                         /* // Waves Growth term
                          double ff = - 0.5*VA[i][e]/dX*0.5*(Pcr[i+1][e] - Pcr[i-1][e])/(pow(B,2)/(8*pi)); 
                          //V3 = ff*dt; // Terme de taux de croissance linéaire 
                          V3 = ff*dt*(log10(Ip[e][i] + 1)/Ip[e][i]); // Terme de taux de croissance avec une fin logarithmique*/


void CRsInjectionSourceSolver(vector<vector<double>> &u_old, vector<vector<double>> &u_new, double dt, vector<double> Pcr_ini, vector<double> Finj_temp, vector<double> vec_theta)
    {
        int NX = u_old.size();
        int NE;
        for (int xi = 0; xi < NX; xi++)
        {   
            NE = u_old[xi].size();
            for (int ei = 0; ei < NE; ei++)
            {
                u_new[xi][ei] = u_old[xi][ei] + dt*Pcr_ini[ei]*Finj_temp[ei]*vec_theta[xi];
            }
        }
    }





void NotMove(vector<vector<double>> u, vector<vector<double>> u_new)
{
    for (int ei = 0; ei < u.size(); ei++)
    {
        for (int xi = 0; xi < u[ei].size(); xi++)
        {
            u_new[xi][ei] = u[xi][ei];
        }
    }
}


