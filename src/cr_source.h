#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;


std::string parameters = "./parameters.dat";
std::string snx = search(parameters, "NX"); int nx = stoi(snx);
std::string sne = search(parameters, "NE"); int ne = stoi(sne);
std::string sni = search(parameters,"ni");  double ni = stod(sni);
std::string sX  = search(parameters,"X");   double Xi = stod(sX); 
double nt = ni/Xi;
std::string smn = search(parameters,"mn");  double m_neutral = stod(smn); 
std::string sT   = search(parameters,"T");  double T  = stod(sT);
std::string scenter = search(parameters, "center"); double x_center = stod(scenter);
std::string sBcenter = search(parameters, "B"); double Bcenter = stod(sBcenter);





double GetTSed()
    {
        double tfree  = 0.3*pow(Esn,-0.5)*Mej*pow(nt,-1./3)*kyr; // [s]
        return tfree;
    }




// R_sh function 
double RSNR(double time)
    {
        //###############################################################################
        //# FUNCTIONS IN ORDER TO MAKE OUR SNR EXPAND IN THE ISM                        #
        //###############################################################################
        double vej8  = 10.*pow(Esn/Mej, 0.5);

        //# We define the characteristic times of our problem
        double tini   = 1e-4*kyr; // [s]
        double tfree  = 0.3*pow(Esn,-0.5)*Mej*pow(nt,-1./3)*kyr; // [s]
        double tPDS   = exp(-1.)*3.61e4*pow(Esn,3./14)/(pow(xi_n,5./14)*pow(nt,4./7))*yr; // [s]
        double tMCS   = min(61*pow(vej8,3)/(pow(xi_n,9./14)*pow(nt,3./7)*pow(Esn,3./14)), 476./pow(xi_n*phi_c,9./14))*tPDS; // [s]
        double tmerge = 153.*pow(pow(Esn,1./14)*pow(nt,1./7)*pow(xi_n,3./14)/(bbeta*C06), 10./7)*tPDS; // [s]
        double tmax = min(tMCS, tmerge); // [s]

        // We define the characteristic radii of our problem
        vector<double> t;
        vector<double> r;
        double R_free  = 5.0*pow(Esn/nt,1./5)*pow(1 - (0.05*pow(Mej,5./6))/(pow(Esn,0.5)*pow(nt,1./3)*(tfree/kyr)), 2./5)*pow(tfree/kyr,2./5)*pc; // [cm]
        double R_ini   = R_free*(tini/tfree); // [cm]
        double R_PDS   = 5.0*pow(Esn/nt,1./5)*pow(1 - (0.05*pow(Mej,5./6))/(pow(Esn,0.5)*pow(nt,1./3)*(tPDS/kyr)), 2./5)*pow(tPDS/kyr,2./5)*pc; // [cm]
        double R_MCS   = R_PDS*pow(tMCS/tPDS, 3./10); // [cm]
        double R_merge = R_MCS*pow(tmerge/tMCS, 1./4); // [cm]
        if (tMCS < tmerge)
        {
            t = {tini, tfree, tPDS, tMCS, tmerge};
            r = {R_ini, R_free, R_PDS, R_MCS, R_merge};
        }
        if (tMCS >= tmerge)
        {
            t = {tini, tfree, tPDS, tmerge};
            r = {R_ini, R_free, R_PDS, R_merge};
        }
    
        vector<double> logt; logt.resize(t.size());
        vector<double> logr; logr.resize(r.size());
        for (int i = 0; i < logt.size(); i++)
        {
            //cout<<"log10(t[i]) = "<<log10(t[i])<<endl;
            //cout<<"log10(r[i]) = "<<log10(r[i])<<endl;
            logt[i] = log10(t[i]);
            logr[i] = log10(r[i]);
        }

        double logr_new, r_new; 
        double logtime = log10(time);
        logr_new = InterpolatingSpline(logt, logr, logtime);
        //cout<<"log_rnew = "<<logr_new<<endl; //-----------------------------------------------------------------------------------
        r_new = pow(10, logr_new);
        //cout<<"r_new = "<<logr_new<<endl;
        return r_new;
    }


// Shock velocity in time 
double u_sh(double time)
    {
        double dt = time/1e4;
        double R_sup = RSNR(time+dt);
        double R_inf = RSNR(time-dt);
        return (R_sup - R_inf)/(2*dt);
    }


// Get the EMAX value for the tesc calculation 
double GetEM()
    {
        int precision = 1000;
        vector<double> Emax; Emax.resize(precision);
        double eps = 1e-4;
        double x0 = 10.*GeV;
        double aa,bb,cc;
        double bbeta = gam-2;
        vector<double> cst; cst.resize(3);

        vector<double> loc_time; loc_time.resize(precision);
        vector<double> log_time; log_time.resize(precision);
        double tmin = 0.01*yr;
        double tmax = 2000*kyr;

        for (int i = 0; i < precision; i++)
            {
                log_time[i] = log10(tmin) + i*(log10(tmax) - log10(tmin))/precision;
                loc_time[i] = pow(10, log_time[i]);
            }





        if (gam-2 != 0)
            {
                for (int i = 0; i < precision; i++)
                    {
                        aa = Emin;
                        bb = bbeta;
                        cc = (bbeta/(1+bbeta))*e*sqrt(4*pi*nt*mp)/(10*c)*xhi_cr*pow(u_sh(loc_time[i]),2)*RSNR(loc_time[i]);
                        //cout<<"cc = "<<RSNR(loc_time[i])<<endl;
                        cst[0] = aa; cst[1] = bb; cst[2] = cc;
                        Emax[i] = NewtonRaphson(f2, df2dx, x0, eps, cst);
                        //cout<<"Emax = "<<Emax[i]/GeV<<" GeV"<<endl; //----------------------------------------------------------------------
                    }

            }
        if (gam-2 == 0)
            {
                for (int i = 0; i < precision; i++)
                    {
                        aa = Emin;
                        bb = e*sqrt(4*pi*nt*mp)/(10*c)*xhi_cr*pow(u_sh(loc_time[i]),2.)*RSNR(loc_time[i]);
                        cst[0] = aa; cst[1] = bb;
                        Emax[i] = NewtonRaphson(f1, df1dx, x0, eps, cst);
                    }
            }
        
        double EM = GetMax(Emax);
        return EM;
    }


// CRs Escape time function !!! 
double tesc(double E)
    {
        double tSed = GetTSed();
        double EM = GetEM();
        double ddelt = -1./(2.*delta); 

        //cout<<"EM = "<<EM/GeV<<" GeV"<<endl;

        //tSed*((E**2/cst.c**2 - cst.mp**2*cst.c**2)/(EMAX**2/cst.c**2 - cst.mp**2*cst.c**2))**(-1./(2.*delta))

        
        //cout<<"(pow(E/c,2) - pow(mp*c,2)) = "<<(pow(E/c,2) - pow(mp*c,2))<<endl;
        //cout<<"(pow(EM/c,2) - pow(mp*c,2)) = "<<(pow(EM/c,2) - pow(mp*c,2))<<endl;
        //cout<<"(pow(E/c,2) - pow(mp*c,2))/(pow(EM/c,2) - pow(mp*c,2) = "<<( pow(E/c,2) - pow(mp*c,2) )/( pow(EM/c,2) - pow(mp*c,2) )<<endl;
        //cout<<"pow( (pow(E/c,2) - pow(mp*c,2)) / (pow(EM/c,2) - pow(mp*c,2)) , ddelt) = "<<pow( ( pow(E/c,2) - pow(mp*c,2) ) / (pow(EM/c,2) - pow(mp*c,2) ) , ddelt )<<endl;
        //cout<<"-1./(2.*delta) = "<<-1./(2.*delta)<<endl;
        //cout<<"E = "<<E/GeV<<" GeV, tesc = "<<tSed*pow( (pow(E/c,2) - pow(mp*c,2)) / (pow(EM/c,2) - pow(mp*c,2)) , ddelt )/kyr<<" kyr"<<endl;
        return tSed*pow( (pow(E/c,2) - pow(mp*c,2)) / (pow(EM/c,2) - pow(mp*c,2)) , ddelt );
    }


double tesc_e(double E)
    {
        double tSed = GetTSed();
        double EM = GetEM();
        double ddelt = -1./(2.*delta); 
        double tesc_p;
        double tloss_e; 

        tesc_p = tSed*pow( (pow(E/c,2) - pow(mp*c,2)) / (pow(EM/c,2) - pow(mp*c,2)) , ddelt );
        tloss_e = 4*E/(c*sig_T*pow(Bcenter,2)*pow(1 + E/(me*pow(c,2)),2));


        //cout<<"Tesc = "<<tesc_p<<" , tloss = "<<tloss_e<<endl;


        return min(tesc_p, tloss_e);
    }


//####################################################################################################//

double Resc(double E)
{
    double A = tesc(E); 
    return RSNR(A);
}

double sigm(double E, double ttesc)
{
    return ttesc/injection_function_width;           // Corresponds approximately to the width of the escape time divided
}


// CRs time injection function. Very important ! 
double Finj(double t, double dt, double E, double ttesc)
{
    double A = exp(-0.5*pow((t - ttesc)/sigm(E, ttesc),2));
    //double A = exp(-0.5*pow((t - 0.01*kyr)/(0.001*kyr),2));
    double B;

    double time_min = t_start_injection; 
    double time_max = t_end_injection*ttesc - time_min;

    int Na = injection_function_norm;                              // Constant in order to easily and rapidly normalize the injection function  
    int Nc = int((time_max - time_min)/dt);

    double ratio = double(Nc)/double(Na);

    double C_a = 0;
    double loc_dt = (time_max - time_min)/double(Na); 
    for (int ti = 0; ti < Na; ti++)
    {
        //C_a += exp(-0.5*pow((t - ttesc)/sigm(E, ttesc),2));
        C_a += exp(-0.5*pow((time_min + ti*loc_dt - ttesc)/sigm(E, ttesc),2))*loc_dt;

    } 
    
    //cout<<ratio<<endl;

    if (C_a == 0){return 0.;}
    if (C_a > 0) 
    {
        //B = A/(C_a*ratio); 
        B = A/C_a;
        if (A > C_a*ratio) {return 0.;}
        if (A < C_a*ratio) {return B;}
    }
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

/*double theta(double z, double t, double Rsnr)
{
    //double Rsnr = RSNR(t);
    

    //cout<<Rsnr/pc<<endl;
    if (z < Rsnr){return 1;}
    if (z > Rsnr){return 0;}
}*/

double theta(double z, double t, double Rsnr)
{
    //double Rsnr = RSNR(t);
    double r_width = Rsnr/r_snr_thickness;
    double center  = x_center;
    double loc_z = z - center;

    //return (1 - tanh((z - Rsnr)/r_width))/2.;
    return 0.5*(erf((Rsnr - loc_z)/(r_width)) + erf((Rsnr + loc_z)/(r_width)));
}

double dNdE(double E)
{
    double spec;
    double Emax = GetEM();
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

double ff(double E)
{
    double R = Resc(E);
    double A = 3*pow(c, 3.)/(16*pow(pi, 2.)*pow(R, 3.))/E*dNdE(E);
    return A;
}

double Pcr_ini(double E)
{
    double A = 4*pi/(3*pow(c, 3))*pow(E, 4.)*ff(E); 
    return A;
}








