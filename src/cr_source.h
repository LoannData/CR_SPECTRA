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
std::string scenter_index = search(parameters, "center_index"); int x_center_index = stod(scenter_index);
std::string sBcenter = search(parameters, "B"); double Bcenter = stod(sBcenter);




/// This function return the value of t_sedov according to the SNR 
/// properties specified in the constants.h file 
double GetTSed()
    {
        double tfree  = 0.3*pow(Esn,-0.5)*Mej*pow(nt,-1./3)*kyr; // [s]
        return tfree;
    }



double RSNR(double time)
    {
        /*! 
            This function provides the value of the shock radius as a function of the time according to the SNR expansion model from 
            Cioffi et al. 1998 and Truelove & McKee 19... The SNR shock has different propagation stages which are separated by 
            characteristic times defined as below : 
            \f{eqnarray*}{ 

                    t_{\mathrm{ini},\mathrm{kyr}}   & = & 10^{-4} \\ 
                    t_{\mathrm{free}, \mathrm{kyr}} & = & 0.3 E_\mathrm{SN}^{-1/2} M_\mathrm{ej} n_t^{-1/3} \\ 
                    t_{\mathrm{PDS}, \mathrm{yr}}   & = & e^{-1} \times 3.61 \times 10^4 E_\mathrm{SN}^{3/14}/(\xi_n^{5/14}n_t^{4/7}) \\ 
                    t_{\mathrm{MCS}, \mathrm{yr}}   & = & \mathrm{min}\left[
                        61 V_{\mathrm{ej},8}^3/(\xi_n^{9/14}*n_t^{3/7}*E_\mathrm{SN}^{3/14}),
                        476 t_{\mathrm{PDS}, \mathrm{yr}} / (\xi_n \phi_c)^{9/14}
                    \right] \\ 
                    t_{\mathrm{merge}, \mathrm{yr}} & = & 153 \left(E_\mathrm{SN}^{1/14} n_t^{1/7} \xi_n^{3/14}/ (\beta C_{06}) \right)^{10/7} t_{\mathrm{PDS}, \mathrm{yr}} \\ 
                    t_\mathrm{max} & = & \mathrm{max}(t_\mathrm{MCS}, t_\mathrm{merge}) 

            \f}    
            where \f$ t_\mathrm{ini}\f$ corresponds to the initial time of the numerical computation of the spline,  
            \f$ t_\mathrm{free} \f$ correponds to the transition time between the free expansion phase and the Sedov-Taylor phase, 
            \f$ t_\mathrm{PDS} \f$ corresponds to the transition time between the Sedov-Taylor phase and the Pressure Driven Snowplow phase, 
            \f$ t_\mathrm{MCS} \f$ correponds to the transition time between the Pressure Driven Snowplow phase and the Momentum Conserving Snowplow phase and 
            finally \f$ t_\mathrm{merge} \f$ means the time at which the shock pressure becomes of the order of the ISM pressure. 
            In some ISM phases, the shock can merge before entering in the Momentum Conserving Snowplow phase explaining \f$ t_\mathrm{max} \f$.
            The SNR shock evolves with time according the following radii 
            \f{eqnarray*}{ 
                R_\mathrm{ini} & = &  R_\mathrm{free} (t_\mathrm{ini}/t_\mathrm{free}) \\ 

            \f}
          
            \param time Time in [s].
            \sa RSNR(), InterpolatingSpline()
        */

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
            t.push_back(1e-6*kyr); 
            t.push_back(tini); 
            t.push_back(tfree); 
            t.push_back(tPDS); 
            t.push_back(tMCS); 
            t.push_back(tmerge); 
            //t = {1e-6*kyr, tini, tfree, tPDS, tMCS, tmerge};
            r.push_back(1e-6*pc);
            r.push_back(R_ini);
            r.push_back(R_free);
            r.push_back(R_PDS);
            r.push_back(R_MCS);
            r.push_back(R_merge);
            //r = {1.e-6*pc ,R_ini, R_free, R_PDS, R_MCS, R_merge};
        }
        if (tMCS >= tmerge)
        {
            t.push_back(1e-6*kyr); 
            t.push_back(tini); 
            t.push_back(tfree); 
            t.push_back(tPDS); 
            t.push_back(tmerge);
            //t = {1e-6*kyr ,tini, tfree, tPDS, tmerge};
            r.push_back(1e-6*pc);
            r.push_back(R_ini);
            r.push_back(R_free);
            r.push_back(R_PDS);
            r.push_back(R_merge);
            //r = {1.e-6*pc, R_ini, R_free, R_PDS, R_merge};
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

        //{cout<<"r_new = "<<r_new; cout<<"logr_new = "<<logr_new<<endl;}

        if (isnan(r_new) ){return 1e-6*pc;}
        if (!(isnan(r_new))){return r_new;}
    }


/// Shock velocity in time 
double u_sh(double time)
    {
        double dt = time/1e4;
        double R_sup = RSNR(time+dt);
        double R_inf = RSNR(time-dt);
        return (R_sup - R_inf)/(2*dt);
    }


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
            if (injection_cutoff == 0){return pow(10,20)*kyr;}
            if (injection_cutoff == 1){return loc_tesc;}
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


double B_sat(double time)
    {
        double B_ISM = Bcenter;

        // See Ohira et al. (2012), eq. 20
        double t_sed  = 0.3*pow(Esn,-0.5)*Mej*pow(nt,-1./3)*kyr; // [s]
        double R_sed  = 5.0*pow(Esn/nt,1./5)*pow(1 - (0.05*pow(Mej,5./6))/(pow(Esn,0.5)*pow(nt,1./3)*(t_sed/kyr)), 2./5)*pow(t_sed/kyr,2./5)*pc; // [cm]
        
        double Bfree = eta_gfree*eta_acc*c*t_sed*Eknee/(3*e*pow(R_sed,2));

        //cout<<"Bfree = "<<Bfree*1e6<<endl; 

        double alpha_B, t_B;
        if (oh_model == 1){alpha_B = alpha - 1./5;}
        if (oh_model == 2){alpha_B = 9./10;}
        if (oh_model == 3){alpha_B = 3./5;}

        t_B = t_sed*pow(Bfree/B_ISM, 1./alpha_B);

        //cout<<"tesc_p = "<<time/kyr<<" kyr"<<endl;

        //cout<<Bfree<<", "<<B_ISM<<", "<<Bfree*pow(time/t_sed, -alpha_B)<<endl;

        if (time <= t_sed)              {return Bfree;}
        if (time > t_sed & time <= t_B) {return Bfree*pow(time/t_sed, -alpha_B);}
        if (time > t_B)                 {return B_ISM;}
    }


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
            if (injection_cutoff == 0){return pow(10,20)*kyr;}
            if (injection_cutoff == 1){return pow(10,20)*kyr;}
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



//####################################################################################################//


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
        if (injection_cutoff == 0){spec = (2 - gam)*Esn*xhi_cr*1e51/(pow(Emax,2-gam) - pow(Emin,2-gam))*pow(E,-gam);}
        if (injection_cutoff == 1){spec = (2 - gam)*Esn*xhi_cr*1e51/(pow(Emax,2-gam) - pow(Emin,2-gam))*pow(E,-gam)*exp(-inj_exp_alpha*E/Emax);}
        
    }
    if (gam == 2)
    {
        if (injection_cutoff == 0){spec = Esn*xhi_cr/(log(Emax) - log(Emin))*pow(E,-gam);}
        if (injection_cutoff == 1){spec = Esn*xhi_cr/(log(Emax) - log(Emin))*pow(E,-gam)*exp(-inj_exp_alpha*E/Emax);}
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








