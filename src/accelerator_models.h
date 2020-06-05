#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;




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

/// Magentic field saturation model 
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