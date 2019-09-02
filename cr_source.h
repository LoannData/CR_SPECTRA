#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;

// Temporaire
//#include "./constants.h"





double RSNR(double time, double ni, double X, double m_neutral, double T)
{
    //cout<<Esn<<endl;
    //Esn = Esn/1e51;
    double n = ni/X; //Total density of the medium [cm^-3]
    double vej8 = 10*pow(Esn/Mej, 0.5); 
    double C06 = 9.79e5*sqrt(5./3/(m_neutral/mp))*sqrt(T*kbolz*6.24e11)/1e6;
    double tfree = 0.3*pow(Esn, -0.5)*Mej*pow(n, -1./3); // [kyr]
    //cout<<tfree<<endl;
    double tsf = 3.61e4*pow(Esn,3./14)/pow(n, 4./7)*1e-3; // [kyr]
    double tPDS = tsf/exp(1.); // [kyr]
    double tMCS = min(61*pow(vej8,3)/(pow(xi_n,9./14)*pow(n,3./7)*pow(Esn, 3./14)), 476/pow(xi_n*phi_c, 9./14)*tPDS); // [kyr]
    double tmerge = 153*pow(pow(Esn, 1./14)*pow(n, 1./7)*pow(xi_n, 3./14)/(beta*C06), 10./7)*tPDS; // [kyr]
    double RSNR, R0, R1; 




    //cout<<"time [kyr] = "<<time/kyr<<endl;
    //cout<<"tfree [kyr] = "<<tfree<<endl;
    if (time/kyr < tfree)
    {
        R0 = 5.0*pow(Esn/n, 1./5)*pow(1 - (0.05*pow(Mej, 5./6)/(pow(Esn, 0.5)*pow(n, 1./3)*tfree)), 2./5)*pow(tfree, 2./5);
        RSNR =  R0*pow(time/kyr/tfree, 1.);
    }
    if (time/kyr > tfree and time/kyr < tPDS)
    {
        RSNR = 5.0*pow(Esn/n, 1./5)*pow(1 - (0.05*pow(Mej, 5./6)/(pow(Esn, 0.5)*pow(n, 1./3)*time/kyr)), 2./5)*pow(time/kyr, 2./5);
    }
    if (time/kyr > tPDS and time/kyr < tMCS)
    {
        R0 = 5.0*pow(Esn/n, 1./5)*pow(1 - (0.05*pow(Mej, 5./6)/(pow(Esn, 0.5)*pow(n, 1./3)*tPDS)), 2./5)*pow(tPDS, 2./5);
        RSNR = R0*pow(time/kyr/tPDS, 3./10);
    }
    if (time/kyr > tMCS and time/kyr < tmerge)
    {
        R0 = 5.0*pow(Esn/n, 1./5)*pow(1 - (0.05*pow(Mej, 5./6)/(pow(Esn, 0.5)*pow(n, 1./3)*tPDS)), 2./5)*pow(tPDS, 2./5);
        R1 = R0*pow(time/kyr/tPDS, 3./10);
        RSNR = R1*pow(time/kyr/tMCS, 1./4);
    }

    return 4.*RSNR*pc; 
}

double tesc(double E, double ni, double X)
{   
    double n = ni/X; //Total density of the medium [cm^-3]
    double b = gam - 2; 
    double A; 
    if (b != 0) {A = 4*sqrt(pi*n)*e/(125*c)*xhi_cr*pow(xhi_0*Esn*1e51/n, 3./5)*b/(1+b)/E*(pow(Emin,b)/(pow(E,b) - pow(Emin,b)));}
    if (b == 0) {A = 4*sqrt(pi*n)*e/(125*c)*xhi_cr*pow(xhi_0*Esn*1e51/n, 3./5)*1/(E*log(E/Emin));}
    return pow(A, 5./4);
}

double Resc(double E, double ni, double X, double m_neutral, double T)
{
    double A = tesc(E, ni, X); 
    return RSNR(A, ni, X, m_neutral, T);
}

double sigm(double E, double ni, double X)
{
    return tesc(E, ni, X)/100.;               // Corresponds approximately to the width of the escape time divided by 100
}


// CRs time injection function. Very important ! 
double Finj(double t, double dt, double E, double ni, double X)
{
    double A = exp(-0.5*pow((t - tesc(E, ni, X))/sigm(E, ni, X),2));
    //double A = exp(-0.5*pow((t - 0.01*kyr)/(0.001*kyr),2));
    double B;

    double time_min = 1e-6*kyr; 
    double time_max = 2*tesc(E, ni, X) - time_min;

    int Na = 100;                              // Constant in order to easily and rapidly normalize the injection function  
    int Nc = int((time_max - time_min)/dt);

    double ratio = double(Nc)/double(Na);

    double C_a = 0;
    for (int ti = 0; ti < Na; ti++)
    {
        C_a += exp(-0.5*pow((t - tesc(E, ni, X))/sigm(E, ni, X),2));
        //C_a += exp(-0.5*pow((t - 0.01*kyr)/(0.001*kyr),2));

    } 
    
    //cout<<A<<endl;

    if (C_a == 0){return 0.;}
    if (C_a > 0) {B = A/(C_a*ratio); return B;}
}




/*double Finj(double t, double E, double ni, double X)
{
    double A = exp(-0.5*pow((t - tesc(E, ni, X))/sigm(E, ni, X),2));

    vector<double> time(1000); 
    double time_0   = 1e-6*kyr; 
    double time_max = 1e4*kyr;
    for (int i = 0; i<time.size(); i++)
    {
        time[i] = time_0 + (time_max - time_0)*i/time.size(); 
    }
    double intA = 0;
    for (int i = 1; i<time.size(); i++)
    {
        intA = intA + exp(-0.5*pow((time[i] - tesc(E, ni, X))/sigm(E, ni, X),2))*(time[i] - time[i-1]); 
    }
    //cout<<intA<<endl;
    //cout<<A/intA<<endl;
    if (intA > 0){return A/intA;}
    if (intA == 0){return 0;}
    //return A/intA;
}*/

double theta(double z, double t, double ni, double X, double m_neutral, double T)
{
    double Rsnr = RSNR(t, ni, X, m_neutral, T);
    //cout<<Rsnr/pc<<endl;
    if (z < Rsnr){return 1;}
    if (z > Rsnr){return 0;}
}

double dNdE(double E)
{
    double spec;
    if (gam != 2)
    {
        spec = (2 - gam)*Esn*1e51/(pow(Emax,2-gam) - pow(Emin,2-gam))*pow(E,-gam);
    }
    if (gam == 2)
    {
        spec = Esn/(log(Emax) - log(Emin))*pow(E,-gam);
    }
    return spec;
}

double ff(double E, double ni, double X, double m_neutral, double T)
{
    double R = Resc(E, ni, X, m_neutral, T);
    double A = 3*pow(c, 3.)/(16*pow(pi, 2.)*pow(R, 3.))/E*dNdE(E);
    return A;
}

double Pcr_ini(double E, double ni, double X, double m_neutral, double T)
{
    double A = 4*pi/(3*pow(c, 3))*pow(E, 4.)*ff(E, ni, X, m_neutral, T); 
    return A;
}




/*double Qcrs(double E, double t, double z, double ni, double X, double m_neutral, double T)
{

    double A = Finj(t, E, ni, X); 
    double B = theta(z, t, ni, X, m_neutral, T);
    double C = Pcr_ini(E, ni, X, m_neutral, T); 
    return A*B*C;
}*/






/*double dNdE(double WCR, double g, double Emin, double Emax, double E)
{
    if (g > 2.)
    {
        return (2-g)*WCR*pow(E,-g)/(pow(Emax, 2-g) - pow(Emin, 2-g));
    }
}

double n(double X, double Xmin, double Xmax)
{
    if (X < Xmin){return 0.;}
    if (X > Xmax){return 0.;}
    if (X >= Xmin && X <= Xmax)
    {
        return 1/(2./3*pi*pow(Xmax - Xmin, 3));
    }
}

double ff(double X, double E, double Xmin, double Xmax, double Emin, double Emax, double WCR, double g)
{
    return n(X, Xmin, Xmax)*pow(c, 3)/(4*pi*pow(E,2))*dNdE(WCR, g, Emin, Emax, E);
}


vector<vector<double>> CRsource1D(vector<double> X, vector<double> E, std::string type)
{
    int NX = X.size();
    int NE = E.size();

    if (type == "classic")
    {

        // We define our CR spectra distribution 
        vector<vector<double>> Pcr(NX); 
        for (int xi = 0; xi < NX; xi++){Pcr[xi].resize(NE);}

        double g = 2.2;
        double WCR = 1.e51*0.1;

        double Xcr_min = 0.*pc;
        double Xcr_max = 24.72*pc;
        double Ecr_min = 1e-4*GeV;
        double Ecr_max = 2e4*GeV;

        for (int xi = 0; xi < NX; xi++)
        {
            for (int ei = 0; ei < NE; ei++)
            {
                double c2 = ff(X[xi], E[ei], Xcr_min, Xcr_max, Ecr_min, Ecr_max, WCR, g);
                Pcr[xi][ei] = 4*pi/3.*c*(pow(E[ei],4)/pow(c,4))*c2; 
            }
        }



    return Pcr;
    }
}


double CRinjection1Dfunction(double time, std::string type)
{
    if (type == "dirac_ini")
    {
        if (time == 0.){return 1.;}
        else {return 0.;}
    }
}*/






/*int main()
{

    double Xmin = 0.*pc;
    double Xmax = 500.*pc;
    int NX = 500;

    vector<double> X(NX);
    for (int xi = 0; xi < NX; xi++)
    {
        X[xi] = Xmin + (xi+1)*(Xmax - Xmin)/NX; 
    }




    double Emin = 1.*GeV;
    double Emax = 10.*TeV;
    int NE = 10;

    vector<double> E(NE); 
    double lEmin = log10(Emin); 
    double lEmax = log10(Emax);

    for (int ei = 0; ei < NE; ei++)
    {
        E[ei] = pow(10.,(lEmin + (ei+1)*(lEmax - lEmin)/NE)); 
        //E[ei] = Emin + ei*(Emax - Emin)/NE;
    }




    vector<vector<double>> data = CRsource1D(X, E, "classic");


    return 0;
}*/