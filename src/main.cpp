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

    std::string parameters = "./parameters.dat";
    std::string snx = search(parameters, "NX"); int nx = stoi(snx);
    std::string sne = search(parameters, "NE"); int ne = stoi(sne);
    std::string sni = search(parameters,"ni");        double ni = stod(sni);
    std::string sX  = search(parameters,"X");         double Xi = stod(sX); 
    std::string smn = search(parameters,"mn");        double m_neutral = stod(smn); 
    std::string sT   = search(parameters,"T");        double T  = stod(sT);




    //cout<<"DIMENSION = "<<DIMENSION<<endl;
    cout<<"NX = "<<nx<<endl;
    cout<<"NE = "<<ne<<endl;
    //cout<<"GRID = "<<GRID<<endl;


    // We get the initial conditions in the ISM
    vector<vector<double>> rVA  = parse2DCsvFile("./data_ini/Alfven.dat", 1);
    vector<vector<double>> rDb  = parse2DCsvFile("./data_ini/DBohm.dat", 1);
    vector<vector<double>> rIp  = parse2DCsvFile("./data_ini/Ip.dat", 1);
    vector<vector<double>> rPcr = parse2DCsvFile("./data_ini/Pcr.dat",1);
    vector<vector<double>> rGd  = parse2DCsvFile("./data_ini/damping.dat",1);





    //We prepare some informations in order to read our input data 
    int NX = pow(2, nx);
    int NE = pow(2, ne); 
    int index_emin = 2; 
    int index_emax = index_emin + NE;

    vector<double> X = readAxis("X", NX);
    vector<double> E = readAxis("E", NE);

    vector<double> B = readAxis("B", NX);



    vector<vector<double>> VA(NX), Db(NX), Ip(NX), Pcr(NX), Gd(NX); 
    for (int xi = 0; xi < NX; xi++)
    {
        VA[xi].resize(NE); 
        Db[xi].resize(NE);
        Ip[xi].resize(NE);
        Pcr[xi].resize(NE);
        Gd[xi].resize(NE);
    }

    for (int xi = 0; xi < NX; xi++)
    {
        for (int ei = 0; ei < NE; ei++)
        {
            VA[xi][ei]  = rVA[xi][index_emin + ei]; 
            Db[xi][ei]  = rDb[xi][index_emin + ei];
            Ip[xi][ei]  = rIp[xi][index_emin + ei]; 
            Pcr[xi][ei] = rPcr[xi][index_emin + ei];
            Gd[xi][ei]  = rGd[xi][index_emin + ei];
            //cout<<"Db["<<xi<<", "<<ei<<"] = "<<Db[xi][ei]<<endl;
        }
    }


    // Initialisation des vecteurs de variance d'espace 
    vector<double> dX(NX), dE(NE); 
    for (int xi=1; xi<NX-1; xi++){dX[xi] = 0.5*abs(X[xi+1] - X[xi-1]);}
    dX[0] = dX[1]; dX[NX-1] = dX[NX-2]; 
    for (int e=1; e<NE-1; e++){dE[e] = 0.5*abs(E[e+1] - E[e-1]);}
    dE[0] = dE[1]; dE[NE-1] = dE[NE-2];

    // CFL Conditions 
    double C1 = 0.5;
    double C2 = 0.5;
    double C3 = 0.5;
    double C4 = 0.5;
    double C5 = 0.5;

    vector<double> Vdminvec(NE), Dminvec(NE), Gammadminvec(NE);

    double maxE  = maxElement1D(E);
    double minW0 = minElement1D(B); minW0 = pow(minW0,2)/(8*pi);
    double maxVd = maxElement2D(VA);
    double maxDb = maxElement2D(Db);
    double minIp = minElement2D(Ip);
    double maxD  = maxDb/minIp;
    double maxGd = maxElement2D(Gd);
    double mindX = minElement1D(dX); double mindE = minElement1D(dE); 

    vector<double> cfl;
    cfl.push_back(C1*mindX/maxVd); 
    cfl.push_back(C2*pow(mindX,2)/maxD);
    cfl.push_back(C3*mindE*mindX/(maxE*maxVd));
    cfl.push_back(C4*2*mindX/maxVd);
    cfl.push_back(C5/maxGd);

    cout<<"CFL Values : "<<endl;
    cout<<"Advection : "<<C1*mindX/maxVd<<" s"<<endl;
    cout<<"Diffusion : "<<C2*pow(mindX,2)/maxD<<" s"<<endl;
    cout<<"Energy    : "<<C3*mindE*mindX/(maxE*maxVd)<<" s"<<endl;
    cout<<"Growth    : "<<C4*2*mindX/maxVd<<" s"<<endl;
    cout<<"Damping   : "<<C5/maxGd<<" s"<<endl;
    double dt = minElement1D(cfl);
    cout<<"Time-step = "<<dt<<" s"<<endl;


    // SIMULATION
    double Tmax = setTmax();
    double time = 0.;
    int time_index = 0;
    vector<double> dat_output = getOutput();

    int out_Pcr_index = 0;
    int out_Ip_index = 0;

    vector<vector<double>> Pcr_new, Pcr_old, Ip_new, Ip_old, Pcr_background;

    Pcr_new = Pcr; Ip_new  = Ip; Pcr_background = Pcr;


    /*vector<vector<double>> Qcrs(NX); 
    for (int xi = 0; xi < NX; xi++)
    {
        Qcrs[xi].resize(NE); 
    }
    Qcrs   = CRsource1D(X, E, "classic");*/ 

    double box_Pcr[5][5], box_Ip[5][5], box_VA[5][5], box_Gd[5][5], box_Db[5][5], box_X[5], box_E[5]; 
    int loc_xi, loc_ei;
    int lxi, lei;
    double P1, P2, Ip1, Ip2;
    double Pcr_new_temp, Ip_new_temp;
    double Pcr_background_temp, B_temp, Qcrs_temp, dX_temp, dE_temp, E_temp;
    int ei_temp, xi_temp, t_NX = NX, t_NE = NE;
    int xi_box, ei_box;
    double Qwaves;

    vector<double> Finj_temp(NE); 
    vector<double> Pcr_ini_temp(NE); 
    double temp_theta; 


    for (int j=0; j<NE; j++)
    {
        Pcr_ini_temp[j] = Pcr_ini(E[j], ni, Xi, m_neutral, T);
    }

    int xi, ei; 
    int nb = nproc;
    int g;

    while (time < Tmax)
    {   
        start = std::clock();

        
        Pcr_old = Pcr_new;
        Ip_old  = Ip_new;

        
        #pragma omp parallel num_threads(nb)
        #pragma omp for schedule(static, int(double(NE/nb))) private (g)
        for (g=0; g<NE; g++)
        {
            Finj_temp[g] = Finj(time, dt, E[g], ni, Xi);
        }

        //cout<<Ip_new[NX-1][0]<<endl;

        #pragma omp parallel num_threads(nb)
        #pragma omp for schedule(static, int(double(NX/nb))) private(xi, ei, xi_box, ei_box, box_X, box_E, box_Pcr, box_Ip, box_Gd, box_Db, box_VA, loc_xi, loc_ei, lxi, lei, P1, P2, Ip1, Ip2, Pcr_background_temp, B_temp, Qcrs_temp, dX_temp, dE_temp, E_temp, ei_temp, xi_temp, Qwaves, Pcr_new_temp, Ip_new_temp, t_NX, t_NE)
        for (xi = 0; xi<NX; xi++)
        {
            temp_theta = theta(X[xi], time, ni, Xi, m_neutral, T);
            for (ei = 0; ei<NE; ei++)
            {
                //We load all the values we need 
                dX_temp = dX[xi];
                dE_temp = dE[ei];
                E_temp = E[ei];
                B_temp = B[xi];
                Pcr_background_temp = Pcr_background[xi][ei];
                for (lxi = 0; lxi < 5; lxi++)
                {
                    for (lei = 0; lei < 5; lei++)
                    {   

                        if (xi-2+lxi >= 0 and xi+2+lxi < NX){xi_box = xi-2+lxi;}
                        if (xi-2+lxi < 0)                   {xi_box = 0;}          //Absorbing left layer 
                        if (xi+2+lxi >= NX)                 {xi_box = NX-1;}       //Absorbing right layer 
                        if (ei-2+lei >= 0 and ei+2+lei < NE){ei_box = ei-2+lei;}
                        if (ei-2+lei < 0)                   {ei_box = 0;}          //Absorbing low energy layer 
                        if (ei+2+lei >= NE)                 {ei_box = NE-1;}       //Absorbing high energy layer 
                                box_Pcr[lxi][lei] = Pcr_old[xi_box][ei_box];
                                box_Ip[lxi][lei]  = Ip_old[xi_box][ei_box]; 
                                box_VA[lxi][lei]  = VA[xi_box][ei_box]; 
                                box_Gd[lxi][lei]  = Gd[xi_box][ei_box]; 
                                box_Db[lxi][lei]  = Db[xi_box][ei_box];
                        box_E[lei] = E[ei_box];
                    }
                    box_X[lxi] = X[xi_box];
                }
                    // EQ. Solving 
                    Qwaves = box_Gd[2][2]*box_Ip[2][2];
                    Qcrs_temp = temp_theta*Finj_temp[ei]*Pcr_ini_temp[ei];

                    Pcr_new_temp = pressureSolver(dt, dX_temp, dE_temp, box_Pcr, box_Ip, box_VA, box_Gd, box_Db, E_temp, Qcrs_temp);
                    Ip_new_temp  = wavesSolver(dt, dX_temp, dE_temp, box_Pcr, box_Ip, box_VA, box_Gd, box_Db, B_temp, Qwaves);
                    if (Pcr_new_temp < Pcr_background_temp) {Pcr_new_temp = Pcr_background_temp;}

                    if (xi == 0){Pcr_new_temp = Pcr_old[xi+1][ei]; Ip_new_temp = Ip_old[xi+1][ei];}
                    if (xi == NX-1){Pcr_new_temp = Pcr_new[xi-1][ei]; Ip_new_temp = Ip_new[xi-1][ei];}

                Pcr_new[xi][ei] = Pcr_new_temp;
                Ip_new[xi][ei]  = Ip_new_temp;
            }
            
        }


        if (time > dat_output[out_Pcr_index] && out_Pcr_index < dat_output.size()) 
        {
            out_Ip_index  = writeXE("Ip", out_Ip_index, Ip_new, NX, NE);
            out_Pcr_index = writeXE("Pcr", out_Pcr_index, Pcr_new, NX, NE);
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        showLog_0(time, Tmax, getLogOutput(), time_index, duration);
        time += dt;
        time_index += 1;

    }
}
