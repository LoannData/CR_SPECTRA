#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
#include <ctime>
#include <chrono>
#include <thread>
#include <sys/resource.h>
#include <omp.h>
using namespace std;


#include "./mathematics.h"
#include "./constants.h"
#include "./read2D.h"
#include "./freader.h"
#include "./fwritter.h"
#include "./out.h"
#include "./logmaker.h"
#include "./tools.h"
#include "./solver1D.h"
#include "./cr_source.h"

// Main function 
int main()
{

    std::clock_t start;
    double duration;
/*
    // Script which defines the limit memory (here 64 MB)
    const rlim_t kStackSize = 64 * 1024 * 1024; 
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }
    //--------------------------------------------------------//
    */

    /*std::string parameters = "./parameters.dat";
    std::string snx = search(parameters, "NX"); int nx = stoi(snx);
    std::string sne = search(parameters, "NE"); int ne = stoi(sne);
    std::string sni = search(parameters,"ni");        double ni = stod(sni);
    std::string sX  = search(parameters,"X");         double Xi = stod(sX); 
    double nt = ni/Xi;
    std::string smn = search(parameters,"mn");        double m_neutral = stod(smn); 
    std::string sT   = search(parameters,"T");        double T  = stod(sT);*/

    


    //cout<<"DIMENSION = "<<DIMENSION<<endl;
    cout<<"NX = "<<nx<<endl;
    cout<<"NE = "<<ne<<endl;
    //cout<<"GRID = "<<GRID<<endl;


    // We get the initial conditions in the ISM
    vector<vector<double>> rVA  = parse2DCsvFile("./data_ini/Alfven.dat", 1);
    vector<vector<double>> rDb  = parse2DCsvFile("./data_ini/DBohm.dat", 1);
    vector<vector<double>> rIp  = parse2DCsvFile("./data_ini/Ip.dat", 1);
    vector<vector<double>> rIm  = parse2DCsvFile("./data_ini/Im.dat", 1);
    vector<vector<double>> rPcr = parse2DCsvFile("./data_ini/Pcr.dat",1);
    vector<vector<double>> rPe  = parse2DCsvFile("./data_ini/Pe.dat",1);
    vector<vector<double>> rGd  = parse2DCsvFile("./data_ini/damping.dat",1);





    //We prepare some informations in order to read our input data 
    int NX = pow(2, nx);
    int NE = pow(2, ne); 
    int index_emin = 2; 
    int index_emax = index_emin + NE;

    vector<double> X = readAxis("X", NX);
    vector<double> E = readAxis("E", NE);

    vector<double> log10E = E;
    for (int ei = 0; ei < E.size(); ei++)
    {
        log10E[ei] = log10(E[ei]);
    }

    vector<double> B = readAxis("B", NX);



    vector<vector<double>> VA(NX), Db(NX), Ip(NX), Im(NX), Pcr(NX), Pe(NX), Gd(NX); 
    for (int xi = 0; xi < NX; xi++)
    {
        VA[xi].resize(NE); 
        Db[xi].resize(NE);
        Ip[xi].resize(NE);
        Im[xi].resize(NE);
        Pcr[xi].resize(NE);
        Pe[xi].resize(NE);
        Gd[xi].resize(NE);
    }

    for (int xi = 0; xi < NX; xi++)
    {
        for (int ei = 0; ei < NE; ei++)
        {
            VA[xi][ei]  = rVA[xi][index_emin + ei]; 
            Db[xi][ei]  = rDb[xi][index_emin + ei];
            Ip[xi][ei]  = rIp[xi][index_emin + ei]; 
            Im[xi][ei]  = rIm[xi][index_emin + ei]; 
            Pcr[xi][ei] = rPcr[xi][index_emin + ei];
            Pe[xi][ei] = rPe[xi][index_emin + ei];
            Gd[xi][ei]  = rGd[xi][index_emin + ei];
            //cout<<"Db["<<xi<<", "<<ei<<"] = "<<Db[xi][ei]<<endl;
        }
    }


    // Initialisation des vecteurs de variance d'espace 
    vector<double> dX(NX), dE(NE); 
    for (int xi=1; xi<NX-1; xi++){dX[xi] = 0.5*abs(X[xi+1] - X[xi-1]);}
    dX[0] = dX[1]; dX[NX-1] = dX[NX-2]; 
    for (int e=1; e<NE-1; e++){dE[e] = 0.5*abs(log10E[e+1] - log10E[e-1]);}
    dE[0] = dE[1]; dE[NE-1] = dE[NE-2];

    // CFL Conditions 
    double C1 = 0.5;
    double C2 = 0.5;
    double C3 = 0.5;
    double C4 = 0.5;
    double C5 = 0.5;

    vector<double> Vdminvec(NE), Dminvec(NE), Gammadminvec(NE);

    vector<vector<double>> abs_VA = VA, dVAdlog10E = VA, cdVAdX = VA, c2dVAdX = VA;
    for (int xi = 0; xi < VA.size(); xi++)
    {
        for (int ei = 0; ei < VA[xi].size(); ei++)
        {
            abs_VA[xi][ei] = abs(VA[xi][ei]);
            if (ei > 0 and ei < VA[xi].size()-1)
            {
                dVAdlog10E[xi][ei] = -1./3*1./log(10)*(VA[xi][ei+1] - VA[xi][ei-1])/(log10E[ei+1] - log10E[ei-1]);
            }
            if (xi > 0 and xi < VA.size()-1)
            {
                cdVAdX[xi][ei] = (VA[xi+1][ei] - VA[xi-1][ei])/(X[xi+1] - X[xi-1])/(3.*log(10.)); //1e-11 -> Test value
                c2dVAdX[xi][ei] = -(VA[xi+1][ei] - VA[xi-1][ei])/(X[xi+1] - X[xi-1]);
            }
            cdVAdX[0][ei] = cdVAdX[1][ei];
            cdVAdX[VA.size()-1][ei] = cdVAdX[VA.size()-2][ei];
            c2dVAdX[0][ei] = c2dVAdX[1][ei];
            c2dVAdX[VA.size()-1][ei] = c2dVAdX[VA.size()-2][ei];
        }
        dVAdlog10E[xi][0] = 0.;
        dVAdlog10E[xi][VA[xi].size()-1] = 0.;

    }


    double maxE  = maxElement1D(log10E);
    double minW0 = minElement1D(B); minW0 = pow(minW0,2)/(8*pi);
    double maxVd = maxElement2D(abs_VA);
    double maxDb = maxElement2D(Db);
    double minIp = minElement2D(Ip);
    double maxIp = maxElement2D(Ip);
    double maxD  = maxDb/minIp;// cout<<"Max D = "<<maxD<<endl;
    double maxGd = maxElement2D(Gd);
    double mindX = minElement1D(dX); double mindE = minElement1D(dE); 

    vector<double> cfl;
    cfl.push_back(C1*mindX/maxVd); 
    //cfl.push_back(C2*pow(mindX,2)/maxD);
    cfl.push_back(C3*mindE*mindX/(maxE*maxVd));
    cfl.push_back(C4*2*mindX/maxVd);
    cfl.push_back(C5/(maxGd));

    cout<<"CFL Values : "<<endl;
    cout<<"Advection : "<<C1*mindX/maxVd<<" s"<<endl;
    cout<<"Diffusion : "<<C2*pow(mindX,2)/maxD<<" s"<<endl;
    cout<<"Energy    : "<<C3*mindE*mindX/(maxE*maxVd)<<" s"<<endl;
    cout<<"Growth    : "<<C4*2*mindX/maxVd<<" s"<<endl;
    cout<<"Damping   : "<<C5/maxGd<<" s"<<endl;
    double dt = minElement1D(cfl);
    cout<<"Time-step = "<<dt<<" s"<<endl;


    // SIMULATION ===========================================================================================================//
    // Simulation output and time variables 
    double Tmax = setTmax();
    double time = 0.;
    int time_index = 0;
    vector<double> dat_output = getOutput();
    int out_Pcr_index = 0;
    int out_Pe_index = 0;
    int out_Ip_index = 0;
    int out_Im_index = 0;

    // Main variables which we solve 
    vector<vector<double>> Pcr_new, Pcr_old, Pcr_background;
    vector<vector<double>> Pe_new, Pe_old, Pe_background;
    vector<vector<double>> Ip_new, Ip_old, Ip_background;
    vector<vector<double>> Im_new, Im_old, Im_background;


    // Initialization of the variables 
    Pcr_background = Pcr;
    Pe_background = Pe;
    Ip_background = Ip;
    Im_background = Im; 
    Pcr_new = Pcr; 
    Pe_new = Pe; 
    Ip_new  = Ip;  
    Im_new  = Im;


    // Variables for CRs injection in the simulation 
    vector<double> Finj_temp(NE);
    vector<double> Finj_temp_e(NE); // Finj of e-  
    vector<double> Pcr_ini_temp(NE); 
    vector<double> Pe_ini_temp(NE);
    vector<double> ttesc(NE);
    vector<double> ttesc_e(NE); // Escape time array of e- 
    vector<double> vec_theta(NX);
    vector<double> vec_theta_e(NX); // vec theta of e-
    double temp_theta, r_snr;


    for (int j=0; j<NE; j++)
    {
        Pcr_ini_temp[j] = Pcr_ini(E[j]);
        ttesc[j] = tesc(E[j]);
        Pe_ini_temp[j] = Pcr_ini(E[j]);                    // A modifier, concerne les e- 
        ttesc_e[j] = tesc_e(E[j]);                           // A modifier, concerne les e- -> C'est fait   
    }



    int xi, ei; 
    //int nb = nproc;
    int g;

    while (time < Tmax)
    {   
        start = std::clock();

        Pcr_old = Pcr_new;
        Pe_old  = Pe_new;
        Ip_old  = Ip_new;
        Im_old  = Im_new;

        for (g=0; g<NE; g++)
        {
            Finj_temp[g] = Finj(time, dt, E[g], ttesc[g]);
            Finj_temp_e[g] = Finj(time, dt, E[g], ttesc_e[g]);
        }
        for (g=0; g<NX; g++)
        {
            vec_theta[g] = theta(X[g], time, r_snr);
            vec_theta_e[g] = theta(X[g], time, r_snr);
        }

        r_snr = RSNR(time);

        
        //----------------------------------------------------------------------//
        // This solver do nothing, just lose your time                          //
        //----------------------------------------------------------------------//
        //NotMove(Ip_old, Ip_new);


        //----------------------------------------------------------------------//
        // Spatially variable diffusion solver, implicit scheme for Pcr and Pe  //
        //----------------------------------------------------------------------//
        if (solver_PcrDiffusion == 1)
        {thetaDiffusionSolver(Pcr_old, Pcr_new, dt, X, NE, Ip_new, Im_new, Db, Pcr_background);        Pcr_old = Pcr_new;}
        if (solver_PeDiffusion  == 1)
        {thetaDiffusionSolver(Pe_old, Pe_new, dt, X, NE, Ip_new, Im_new, Db, Pe_background);        Pe_old = Pe_new;}


        //----------------------------------------------------------------------//
        // Explicit Advection solver for Pcr and Pe by Alfvén velocity          //
        //----------------------------------------------------------------------//
        if (solver_PcrAdvection == 1)
        {advectionSolverX(Pcr_old, Pcr_new, dt, X, NE, VA, 0);                                  Pcr_old = Pcr_new;}
        if (solver_PeAdvection == 1)
        {advectionSolverX(Pe_old, Pe_new, dt, X, NE, VA, 0);                                  Pe_old = Pe_new;}

        //----------------------------------------------------------------------//
        // Explicit Advection solver for I by Alfvén velocity                   //
        // -> This term seems ok                                                //
        //----------------------------------------------------------------------//
        if (solver_IpAdvection == 1)
        {advectionSolverX(Ip_old, Ip_new, dt, X, NE, VA,  1);                                   Ip_old = Ip_new;}

        if (solver_ImAdvection == 1)
        {advectionSolverX(Im_old, Im_new, dt, X, NE, VA, -1);                                   Im_old = Im_new;}


        //----------------------------------------------------------------------//
        // Explicit Advection solver for Pcr and Pe by the energy derivative of //
        // Alfvén velocity. Which seems to have no effects on the diffusion     //
        // Note : This term needs more studies ...                              //
        //----------------------------------------------------------------------//
        if (solver_PcrAdvection2 == 1)
        {advectionSolverX(Pcr_old, Pcr_new, dt, X, NE, dVAdlog10E, 0);                          Pcr_old = Pcr_new;}
        if (solver_PeAdvection2 == 1)
        {advectionSolverX(Pe_old, Pe_new, dt, X, NE, dVAdlog10E, 0);                          Pe_old = Pe_new;}

        //----------------------------------------------------------------------//
        // Explicit Advection solver for Pcr and Pe in energy cdVAdX
        // This value needs to be studied more in details 
        //----------------------------------------------------------------------//
        if (solver_PcrAdvectionE == 1)
        {advectionSolverE(Pcr_old, Pcr_new, dt, log10E, NX, cdVAdX, Pcr_background);         Pcr_old = Pcr_new;}
        if (solver_PeAdvectionE == 1)
        {advectionSolverE(Pe_old, Pe_new, dt, log10E, NX, cdVAdX, Pe_background);         Pe_old = Pe_new;}

        //----------------------------------------------------------------------//
        // Source term from synchrotron radiation of e-. It contains a pure 
        // source term and an energy advective term
        //----------------------------------------------------------------------//
        if (solver_PeAdvectionE2 == 1)
        {advectionSolverE2(Pe_old, Pe_new, dt, log10E, NX, B, E, Pe_background);         Pe_old = Pe_new;}


        //----------------------------------------------------------------------//
        // Source term effect due to the dependance of the Alfvén velocity to   //
        // the space (for Pcr and Pe)                                           //
        //----------------------------------------------------------------------//
        if (solver_PcrSource1 == 1)
        {sourceSolver(Pcr_old, Pcr_new, dt, c2dVAdX, 4./3);                                  Pcr_old = Pcr_new;}
        if (solver_PeSource1 == 1)
        {sourceSolver(Pe_old, Pe_new, dt, c2dVAdX, 4./3);                                  Pe_old = Pe_new;}

        //----------------------------------------------------------------------//
        // Source term effect due to the dependance of the Alfvén velocity to   //
        // the space                                                            //
        //----------------------------------------------------------------------//
        if (solver_IpSource1 == 1)
        {sourceSolver(Ip_old, Ip_new, dt, c2dVAdX, 1.);                                      Ip_old = Ip_new;}
        if (solver_ImSource1 == 1)
        {sourceSolver(Im_old, Im_new, dt, c2dVAdX, 1.);                                      Im_old = Im_new;}

        //----------------------------------------------------------------------//
        // Source term effect due to production of self-turbulence - damping    //
        // -> This term seems ok                                                //
        //----------------------------------------------------------------------//
        if (solver_IpDampGrowth == 1)
        {sourceGrowthDampRateSolver(Ip_old, Ip_new, Pcr_old, Gd, Ip_background, X, dt, VA, B,  1); Ip_old = Ip_new;}
        if (solver_ImDampGrowth == 1)
        {sourceGrowthDampRateSolver(Im_old, Im_new, Pcr_old, Gd, Im_background, X, dt, VA, B, -1); Im_old = Im_new;}

        //----------------------------------------------------------------------//
        // CRs and e- injection term from SNRs                                  // 
        //----------------------------------------------------------------------//
        if (solver_PcrSource2 == 1)
        {CRsInjectionSourceSolver(Pcr_old, Pcr_new, dt, Pcr_ini_temp, Finj_temp, vec_theta); Pcr_old = Pcr_new; }
        if (solver_PeSource2 == 1)
        {CRsInjectionSourceSolver(Pe_old, Pe_new, dt, Pe_ini_temp, Finj_temp_e, vec_theta_e); Pe_old = Pe_new; }




        if (time > dat_output[out_Pcr_index] && out_Pcr_index < dat_output.size()) 
        {
            out_Ip_index  = writeXE("Ip", out_Ip_index, Ip_new, NX, NE);
            out_Im_index  = writeXE("Im", out_Im_index, Im_new, NX, NE);
            out_Pcr_index = writeXE("Pcr", out_Pcr_index, Pcr_new, NX, NE);
            out_Pe_index = writeXE("Pe", out_Pe_index, Pe_new, NX, NE);
        }

        
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        duration = duration / nproc;
        showLog_0(time, Tmax, getLogOutput(), time_index, duration);
        time += dt;
        time_index += 1;

    }
}
