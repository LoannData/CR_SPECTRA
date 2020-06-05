#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;





//===============================================================//
// PROTONS SOURCE
//===============================================================//

/// Get the EMAX value for the tesc calculation 
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



/// Protons escape time function according the model of Celli et al. (2019)
/// but also alternatives models of CR escape during radiative stages
double tesc(double E)
    {
        double tSed = GetTSed();
        double EM = GetEM();
        double ddelt = -1./(delta); 

        double tPDS   = exp(-1.)*3.61e4*pow(Esn,3./14)/(pow(xi_n,5./14)*pow(nt,4./7))*yr; // [s]
        
        double loc_tesc = tSed*pow( (pow(E/c,2)) / (pow(EM/c,2)) , ddelt);
        double tfree  = 0.3*pow(Esn,-0.5)*Mej*pow(nt,-1./3)*kyr; // [s]

        double v_sh = u_sh(loc_tesc);

        double t_sh_lib   = tSed;
        double v_sh_110 = u_sh(t_sh_lib); 
        double delta_t    = 100.*yr;
        // We calculate the time at which u_sh = 110 km/s 
        while(v_sh_110 > 110*km)
        {
            t_sh_lib += delta_t;
            v_sh_110 = u_sh(t_sh_lib); 
        }

        if (E > EM || E < Emin)
        {
            if (injection_law_protons == 0){return pow(10,20)*kyr;}
            if (injection_law_protons == 1){return loc_tesc;}
        }



        if (tesc_model == 0) // no radiative escape model 
        {
            return loc_tesc;
        }
        if (tesc_model == 1) // if t >= tPDS
        {
            return min(tPDS, loc_tesc); 
        }
        if (tesc_model == 2) // if vsh <= 110 km.s^-1
        {
            if (loc_tesc > tPDS)
            {
                return t_sh_lib;
            }
            else 
            {
                return loc_tesc;
            }
        }
    }


/// This function return the value of the escape shock radius at a given energy of protons
double Resc(double E)
{
    double A = tesc(E); 
    return RSNR(A);
}

/// This function defines the variance of the time injection function Finj of CRs 
double sigm(double E, double ttesc)
{
    return ttesc/injection_function_width;           // Corresponds approximately to the width of the escape time divided
}


/// CRs time injection function. 
double Finj(double t, double dt, double E, double ttesc)
{
    if (injection_shape_time == 0) // Dirac type injection 
    {
        if (ttesc > t - dt/2. && ttesc < t + dt/2.){
            if (verbose == 1) {cout<<"E = "<<E/GeV<<" GeV, t_min = "<<(t - dt/2.)/kyr<<" < "<<ttesc/kyr<<" < t_max = "<<(t + dt/2.)/kyr<<" kyr"<<endl;} 
            return 1./dt; }
        return 0.;
    }

    if (injection_shape_time == 1) // Gaussian type injection 
    {
        double A = exp(-0.5*pow((t - ttesc)/sigm(E, ttesc),2));
        //double A = exp(-0.5*pow((t - 0.01*kyr)/(0.001*kyr),2));
        double B;

        double time_min = t_start_injection; 
        double time_max = t_end_injection*ttesc;// - time_min;

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
}


/// This function define the spatial shape of the CRs distribution at a given time 
double theta(double z, double t, double Rsnr)
{
    //double Rsnr = RSNR(t);
    double r_width = Rsnr/r_snr_thickness;
    double center  = x_center;
    double loc_z = z - center;

    //return (1 - tanh((z - Rsnr)/r_width))/2.;
    return 0.5*(erf((Rsnr - loc_z)/(r_width)) + erf((Rsnr + loc_z)/(r_width)));
}

/// This function defines the injection spectrum of both electrons of protons escaping from the SNR
double dNdE(double E)
{
    double spec;
    double Emax = GetEM();
    if (gam != 2)
    {
        if (injection_law_protons == 0){spec = (2 - gam)*Esn*xhi_cr*1e51/(pow(Emax,2-gam) - pow(Emin,2-gam))*pow(E,-gam);}
        if (injection_law_protons == 1){spec = (2 - gam)*Esn*xhi_cr*1e51/(pow(Emax,2-gam) - pow(Emin,2-gam))*pow(E,-gam)*exp(-inj_exp_alpha*E/Emax);}
        
    }
    if (gam == 2)
    {
        if (injection_law_protons == 0){spec = Esn*xhi_cr/(log(Emax) - log(Emin))*pow(E,-gam);}
        if (injection_law_protons == 1){spec = Esn*xhi_cr/(log(Emax) - log(Emin))*pow(E,-gam)*exp(-inj_exp_alpha*E/Emax);}
    }
    return spec;
}

/// This function defines the distribution function of both electrons and protrons escaping from the SNR
double ff(double E)
{
    double R = Resc(E);
    double A = 3*pow(c, 3.)/(16*pow(pi, 2.)*pow(R, 3.))/pow(E,2)*dNdE(E);
    return A;
}

/// This function defines the initial pressure distribution of both electrons and protons escaping from the SNR
double Pcr_ini(double E)
{
    double A = 4*pi/(3*pow(c, 3))*pow(E, 4.)*ff(E); 
    return A;
}








//===============================================================//
// ELECTRONS SOURCE
// Note : Some of electrons functions are in common with those 
// used in the proton source
//===============================================================//

/// This function allows to calculate the maximum energy provided to the electrons escaping from the SNR according to 
/// the values in the constans.h file.
double GetEM_e()
    {
        double t_sed  = 0.3*pow(Esn,-0.5)*Mej*pow(nt,-1./3)*kyr; // [s]
        double Bsat = B_sat(t_sed);
        double EM = 9*pow(me, 4)*pow(c, 7)/(4*16*pow(e, 4)*pow(Bsat,2)*t_sed);
        //cout<<"EM electrons = "<<EM/GeV<<endl;
        double EM_p = GetEM();

        return min(EM, EM_p);
    }

/// This function calulate the escape time of the electrons as a function of the energy according 
/// both model of escape from Celli et al. (2019) and alternative models of escape in radiative phases
double tesc_e(double E)
    {
        double tSed = GetTSed();
        double EM = GetEM_e();
        double EM_p = GetEM();
        double ddelt = -1./(delta); 

        double tPDS   = exp(-1.)*3.61e4*pow(Esn,3./14)/(pow(xi_n,5./14)*pow(nt,4./7))*yr; // [s]
        
        //double loc_tesc = tSed*pow( (pow(E/c,2) - pow(mp*c,2)) / (pow(EM/c,2) - pow(mp*c,2)) , ddelt);
        double loc_tesc = tSed*pow( (pow(E/c,2)) / (pow(EM_p/c,2)) , ddelt);
        double tfree  = 0.3*pow(Esn,-0.5)*Mej*pow(nt,-1./3)*kyr; // [s]

        double v_sh = u_sh(loc_tesc);

        double t_sh_lib   = tSed;
        double v_sh_110 = u_sh(t_sh_lib); 
        double delta_t    = 100.*yr;
        // We calculate the time at which u_sh = 110 km/s 
        while(v_sh_110 > 110*km)
        {
            t_sh_lib += delta_t;
            v_sh_110 = u_sh(t_sh_lib); 
        }

        if (E > EM || E < Emin)
        {
            if (injection_law_electrons == 0){return pow(10,20)*kyr;}
            if (injection_law_electrons == 1){return loc_tesc;}
        }


        if (tesc_model == 0) // no radiative escape model 
        {
            return loc_tesc; 
        }
        if (tesc_model == 1) // if t >= tPDS
        {
            return min(tPDS, loc_tesc); 
        }
        if (tesc_model == 2) // if vsh <= 110 km.s^-1
        {
            if (loc_tesc > tPDS)
            {
                return t_sh_lib;
            }
            else 
            {
                return loc_tesc;
            }
        }
    }

/// This function return the value of the escape shock radius at a given energy of protons
double Resc_e(double E)
{
    double A = tesc_e(E); 
    return RSNR(A);
}



/// This function defines the injection spectrum of both electrons of protons escaping from the SNR
double dNdE_e(double E)
{
    double spec;
    double Emax = GetEM_e();
    if (gam != 2)
    {
        if (injection_law_electrons == 0){spec = (2 - gam_e)*Esn*xhi_cr*1e51/(pow(Emax,2-gam_e) - pow(Emin,2-gam_e))*pow(E,-gam_e);}
        if (injection_law_electrons == 1){spec = (2 - gam_e)*Esn*xhi_cr*1e51/(pow(Emax,2-gam_e) - pow(Emin,2-gam_e))*pow(E,-gam_e)*exp(-inj_exp_alpha_e*E/Emax);}
        
    }
    if (gam == 2)
    {
        if (injection_law_electrons == 0){spec = Esn*xhi_cr/(log(Emax) - log(Emin))*pow(E,-gam_e);}
        if (injection_law_electrons == 1){spec = Esn*xhi_cr/(log(Emax) - log(Emin))*pow(E,-gam_e)*exp(-inj_exp_alpha_e*E/Emax);}
    }
    return spec;
}

/// This function defines the distribution function of both electrons and protrons escaping from the SNR
double ff_e(double E)
{
    double R = Resc_e(E);
    double A = 3*pow(c, 3.)/(16*pow(pi, 2.)*pow(R, 3.))/pow(E,2)*dNdE_e(E);
    return A;
}

/// This function defines the initial pressure distribution of both electrons and protons escaping from the SNR
double Pe_ini(double E)
{
    double A = 4*pi/(3*pow(c, 3))*pow(E, 4.)*ff_e(E); 
    return A;
}











