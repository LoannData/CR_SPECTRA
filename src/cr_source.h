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
std::string sni = search(parameters,"ni");        double ni = stod(sni);
std::string sX  = search(parameters,"X");         double Xi = stod(sX); 
double nt = ni/Xi;
std::string smn = search(parameters,"mn");        double m_neutral = stod(smn); 
std::string sT   = search(parameters,"T");        double T  = stod(sT);



// TriDiagonal matrix inversion function
void InverseTrigonalMatrix(vector<vector<double>> &T)
    {
    int n = T.size();

    

    double a[n];
    double b[n];
    double c[n];
    for (int i = 0; i < n; i++)
        {
            b[i] = T[i][i];
            if (i < n-1){c[i] = T[i][i+1];}
            if (i > 0){a[i] = T[i][i-1];}
        }
    

    double theta_mm = 0.;
    double theta_m = 1.;
    double theta[n];
    theta[0] = b[0]*theta_m;
    theta[1] = b[1]*theta[0] - a[1]*c[0]*theta_m;
    for (int i = 2; i < n; i++)
        {
            theta[i] = b[i]*theta[i-1] - a[i]*c[i-1]*theta[i-2];
        }
    
    double phi_pp = 0;
    double phi_p  = 1;
    double phi[n+2];
    phi[n] = phi_p;
    phi[n+1] = phi_pp;
    phi[n-1] = b[n-1]*phi_p;
    phi[n-2] = b[n-2]*phi[n-1] - c[n-2]*a[n-1]*phi_p;
    for (int i = n-3; i > -1; i--)
        {
            phi[i] = b[i]*phi[i+1] - c[i]*a[i+1]*phi[i+2];
        }
    
    vector<vector<double>> Tinv;
    double Tij;
    Tinv.resize(n);
    for (int i = 0; i < n; i++)
        {
            
            Tinv[i].resize(n);
            for (int j = 0; j < n; j++)
                {
                    //cout<<i<<", "<<j<<endl;
                    if (i < j)
                        {
                            double p = 1.;
                            double loc_theta;
                            for (int ei = i; ei < j; ei++)
                                {
                                    p = p*c[ei];
                                }
                            if (i-1 == -1){loc_theta = theta_m;}
                            if (i-1 > 1){loc_theta = theta[i-1];}
                            Tij = pow(-1,i+j)*p*loc_theta*phi[j+1]/theta[n-1];
                        }
                    if (i == j)
                        {
                            double loc_theta;
                            if (i-1 == -1){loc_theta = theta_m;}
                            if (i-1 > -1){loc_theta = theta[i-1];}
                            Tij = loc_theta*phi[j+1]/theta[n-1];
                        }
                    if (i > j)
                        {
                            double p = 1.;
                            double loc_theta;
                            for (int ei = j+1; ei < i+1; ei++)
                                {
                                    p = p*a[ei];
                                }
                            
                            if (j-1 == -1){loc_theta = theta_m;}
                            if (j-1 > -1){loc_theta = theta[j-1];}
                            
                            Tij = pow(-1, i+j)*p*loc_theta*phi[i+1]/theta[n-1];
                            //cout<<"Hello"<<endl;
                        }
                    Tinv[i][j] = Tij;
                }
        }
        T = Tinv;
    }

void ProductMatrix(vector<vector<double>> A, vector<vector<double>> B, vector<vector<double>> &C)
    {
        int A_l = A.size();
        int A_c = A[0].size();

        int B_l = B.size();
        int B_c = B[0].size();

        int C_l = A_l;
        int C_c = B_c;

        
        C.resize(C_l);
        for (int i = 0; i < C_l; i++)
            {
                C[i].resize(C_c);
                for (int j = 0; j < C_c; j++)
                {
                    double s = 0.;
                    for (int k = 0; k < A_c; k++)
                        {
                            s = s + A[i][k]*B[k][j];
                        }
                    C[i][j] = s;
                }
            }
    }

double InterpolatingSpline(vector<double> X, vector<double> Y, double x)
    {
        int N = X.size();

        
        double h[N-1];
        for (int i = 0; i < N-1; i++)
            {
                h[i] = X[i+1] - X[i];
            }
        
        
        double lowF = 0.;
        double highF = 0.;
        vector<vector<double>> F; 
        F.resize(N);
        for (int i = 0; i < N; i++)
            {
                F[i].resize(1);
            }
        F[0][0] = lowF;
        for (int i = 1; i < N-1; i++)
            {
                F[i][0] = (Y[i+1] - Y[i])/h[i] - (Y[i] - Y[i-1])/h[i-1];
            }
        F[N-1][0] = highF;

        
        vector<vector<double>> R;
        R.resize(N);
        for (int i = 0; i < N; i++)
            {
                R[i].resize(N);
            }
        R[0][0] = 1.;
        R[N-1][N-1] = 1.;
        for (int i = 1; i < N-1; i++)
            {
                R[i][i] = (h[i-1] + h[i])/3.;
                R[i][i+1] = h[i]/6.;
                R[i][i-1] = h[i-1]/6.;
            }
        
        
        vector<vector<double>> Rinv;
        InverseTrigonalMatrix(R);
        Rinv = R;

        vector<vector<double>> M; 
        
        ProductMatrix(Rinv, F, M);
        
        double Ci[N-1]; 
        double Cip[N-1];
        for (int i = 0; i < N-1; i++)
            {
                Ci[i]  = (Y[i+1] - Y[i])/h[i] - h[i]*(M[i+1][0] - M[i][0])/6.;
                Cip[i] = Y[i] - M[i][0]*pow(h[i],2)/6.;
            }
        
        double loc_a, loc_b, loc_c;
        for (int k = 0; k < N-1; k++)
            {
                if (x >= X[k] and x <= X[k+1])
                    {
                        loc_a = M[k][0]*pow(X[k+1] - x,3)/(6*h[k]);
                        loc_b = M[k+1][0]*pow(x - X[k], 3)/(6*h[k]);
                        loc_c = Ci[k]*(x-X[k]) + Cip[k];
                    }
            }
        
        return loc_a + loc_b + loc_c;
    }


// Functions in order to calculate Emax(t) 
// With Newton-Raphson method 
double f1(double x, vector<double> cst)
    {
        double a = cst[0];
        double c = cst[1];
        return x*log(x/a) - c;
    }

double df1dx(double x, vector<double> cst)
    {
        double a = cst[0];
        return log(x/a);        
    }

double f2(double x, vector<double> cst)
    {
        double a = cst[0];
        double b = cst[1];
        double c = cst[2];
        return x*(pow(x/a, b) -1) - c;
    }

double df2dx(double x, vector<double> cst)
    {
        double a = cst[0];
        double b = cst[1];
        return b*pow(x/a,b) + pow(x/a,b) - 1;
    }


double NewtonRaphson(double f(double, vector<double>), double df(double, vector<double>), double x0, double eps, vector<double> cst)
    {
        double exp = 1e30;
        int niter  = 0;
        double x = x0; 
        double xold;
        while (exp > eps || niter < 50)
        {
            xold = x;
            x = x - f(x, cst)/df(x, cst);
            exp = abs(xold-x)/x;
            niter ++;
        }
        return x;
    }

double GetMax(vector<double> V)
    {
        double max = 0.;
        for (int i = 0; i < V.size(); i++)
            {
                if (V[i] >= max)
                    {
                        max = V[i];
                    }
            }
        return max;
    }






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
        double tmerge = 153.*pow(pow(Esn,1./14)*pow(nt,1./7)*pow(xi_n,3./14)/(beta*C06), 10./7)*tPDS; // [s]
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
            logt[i] = log10(t[i]);
            logr[i] = log10(r[i]);
        }

        double logr_new, r_new; 
        logr_new = InterpolatingSpline(logt, logr, time);
        r_new = pow(10, logr_new);
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
                        cst[0] = aa; cst[1] = bb; cst[2] = cc;
                        Emax[i] = NewtonRaphson(f2, df2dx, x0, eps, cst);
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
        return tSed*pow((pow(E,2)/pow(c,2) - pow(mp*c,2))/(pow(EM/c,2) - pow(mp*c,2)), -1./(2.*delta));
    }


//####################################################################################################//

double Resc(double E)
{
    double A = tesc(E); 
    return RSNR(A);
}

double sigm(double E)
{
    return tesc(E)/injection_function_width;           // Corresponds approximately to the width of the escape time divided
}


// CRs time injection function. Very important ! 
double Finj(double t, double dt, double E, double ttesc)
{
    double A = exp(-0.5*pow((t - ttesc)/sigm(E),2));
    //double A = exp(-0.5*pow((t - 0.01*kyr)/(0.001*kyr),2));
    double B;

    double time_min = t_start_injection; 
    double time_max = t_end_injection*ttesc - time_min;

    int Na = injection_function_norm;                              // Constant in order to easily and rapidly normalize the injection function  
    int Nc = int((time_max - time_min)/dt);

    double ratio = double(Nc)/double(Na);

    double C_a = 0;
    for (int ti = 0; ti < Na; ti++)
    {
        C_a += exp(-0.5*pow((t - ttesc)/sigm(E),2));
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

double theta(double z, double t)
{
    double Rsnr = RSNR(t);
    //cout<<Rsnr/pc<<endl;
    if (z < Rsnr){return 1;}
    if (z > Rsnr){return 0;}
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








