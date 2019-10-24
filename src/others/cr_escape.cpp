#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std; 


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

double GetTesc(double E, double delta, vector<double> cst)
    {
        double tSed = cst[2];
        double EM = cst[3]; 
        double c = cst[0];
        double mp = cst[1];

        return tSed*pow((pow(E,2)/pow(c,2) - pow(mp*c,2))/(pow(EM/c,2) - pow(mp*c,2)), -1./(2.*delta));
    }





int main()
    {

// Main constants 
double pi = 3.14720;
double c  = 2.997e10; 
double yr = 362.25*86400;
double kyr = yr*1e3;
double pc   = 3.086e18;   // 1 pc in [cm]
double GeV  = 0.00160218;
double e = 4.80326e-10; 
double mp = 1.6726219e-24;

// We define the constants of our problem
double nt    = 0.35;  // [atom/cm^-3] Density, here WNM
double xhi_m = 1; // Metallicity (1 for solar abundances)
double xhi_cr= 0.1;
double E51   = 1; // [10^51 erg] Total energy released by the SNR
double Mej   = 1; // [Msum] Total mass released by the SNR
double C06   = 1; // 
double beta  = 1; //
double phi_c = 1; // Ratio of the thermal conductivity to the Spitzer (1962) value
double vej8  = 10.*pow(E51/Mej, 0.5);

//###############################################################################
//# FUNCTIONS IN ORDER TO MAKE OUR SNR EXPAND IN THE ISM                        #
//###############################################################################
//# We define the characteristic times of our problem
double tini   = 1e-4*kyr; // [s]
double tfree  = 0.3*pow(E51,-0.5)*Mej*pow(nt,-1./3)*kyr; // [s]
double tPDS   = exp(-1.)*3.61e4*pow(E51,3./14)/(pow(xhi_m,5./14)*pow(nt,4./7))*yr; // [s]
double tMCS   = min(61*pow(vej8,3)/(pow(xhi_m,9./14)*pow(nt,3./7)*pow(E51,3./14)), 476./pow(xhi_m*phi_c,9./14))*tPDS; // [s]
double tmerge = 153.*pow(pow(E51,1./14)*pow(nt,1./7)*pow(xhi_m,3./14)/(beta*C06), 10./7)*tPDS; // [s]
double tmax = min(tMCS, tmerge); // [s]


// We define the characteristic radii of our problem
vector<double> t;
vector<double> r;
double R_free  = 5.0*pow(E51/nt,1./5)*pow(1 - (0.05*pow(Mej,5./6))/(pow(E51,0.5)*pow(nt,1./3)*(tfree/kyr)), 2./5)*pow(tfree/kyr,2./5)*pc; // [cm]
double R_ini   = R_free*(tini/tfree); // [cm]
double R_PDS   = 5.0*pow(E51/nt,1./5)*pow(1 - (0.05*pow(Mej,5./6))/(pow(E51,0.5)*pow(nt,1./3)*(tPDS/kyr)), 2./5)*pow(tPDS/kyr,2./5)*pc; // [cm]
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

int loc_size = 100;
vector<double> logt_new; logt_new.resize(loc_size);
vector<double> logr_new; logr_new.resize(loc_size);
for (int i = 0; i < loc_size; i++)
{
    logt_new[i] = logt[0] + i*(logt[logt.size()-1] - logt[0])/loc_size;
    logr_new[i] = InterpolatingSpline(logt, logr, logt_new[i]);
}

vector<double> t_new; t_new.resize(logt_new.size());
vector<double> r_new; r_new.resize(logr_new.size());
for (int i = 0; i < loc_size; i++)
{
    t_new[i] = pow(10, logt_new[i]);
    r_new[i] = pow(10, logr_new[i]);
}



// We calculate the shock velocity 
vector<double> u_sh; u_sh.resize(t_new.size());
for (int i = 1; i < u_sh.size(); i++)
    {
        u_sh[i] = (r_new[i] - r_new[i-1])/(t_new[i] - t_new[i-1]);
    }
u_sh[0] = u_sh[1] - (u_sh[2] - u_sh[1]);



double gamma = 2.2;
double bbeta = gamma - 2.;
double Emin = 0.1*GeV;

vector<double> Emax; Emax.resize(t_new.size());
double eps = 1e-4;
double x0 = 10.*GeV;
double a,b,cc;
vector<double> cst; cst.resize(3);

if (beta != 0)
    {
        for (int i = 0; i < t_new.size(); i++)
            {
                a = Emin;
                b = bbeta;
                cc = (bbeta/(1+bbeta))*e*sqrt(4*pi*nt*mp)/(10*c)*xhi_cr*pow(u_sh[i],2)*r_new[i];
                cst[0] = a; cst[1] = b; cst[2] = cc;
                Emax[i] = NewtonRaphson(f2, df2dx, x0, eps, cst);
            }
    }
if (beta == 0)
    {
        for (int i = 0; i < t_new.size(); i++)
            {
                a = Emin;
                b = e*sqrt(4*pi*nt*mp)/(10*c)*xhi_cr*pow(u_sh[i],2.)*r_new[i];
                cst[0] = a; cst[1] = b;
                Emax[i] = NewtonRaphson(f1, df1dx, x0, eps, cst);
            }
    }


//for (int i = 0; i < Emax.size(); i++)
//    {
//        cout<<"t_new = "<<t_new[i]/kyr<<" kyr, Emax = "<<Emax[i]/GeV<<" GeV"<<endl;
//    }

double EM = GetMax(Emax);
double delta = 3.;

double tSed = tfree;


vector<double>  cst2 = {c, mp, tSed, EM};
double Etest = 100.*GeV;
double tesc_test = GetTesc(Etest, delta, cst2);

cout<<"Etest = "<<Etest/GeV<<" GeV, tesc = "<<tesc_test/kyr<<" kyr"<<endl;





    return 0;    
    }