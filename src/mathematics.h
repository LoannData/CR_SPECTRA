#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;


// Thomas Algorithm (Trigonal matrix inversion O(N))
vector<double> TDMA(vector<double> a, vector<double> b, vector<double> c, vector<double> d)
    {
        int n = d.size();
        vector<double> w(n-1, 0);
        vector<double> g(n, 0);
        vector<double> p(n, 0);

        w[0] = c[0]/b[0];
        g[0] = d[0]/b[0];

        for (int i = 1; i < n-1; i++)
            {
                w[i] = c[i]/(b[i] - a[i-1]*w[i-1]);
            }
        for (int i = 1; i < n; i++)
            {
                g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1]);
            }
        p[n-1] = g[n-1]; 
        for (int i = n-1; i > 0; i--)
            {
                p[i-1] = g[i-1] - w[i-1]*p[i]; 
            }
        return p;

    }



// TriDiagonal matrix inversion function O(N^3)
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

        //cout<<"x = "<<x<<endl;
        
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
        
        //cout<<"x = "<<x<<endl;
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
            //cout<<"x = "<<x<<endl;
            //cout<<
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

void transpose(vector<vector<double>> &b, vector<vector<double>> &c)
{
    if (b.size() == 0)
        return;

    vector<vector<double> > trans_vec(b[0].size(), vector<double>());

    for (int i = 0; i < b.size(); i++)
    {
        for (int j = 0; j < b[i].size(); j++)
        {
            trans_vec[j].push_back(b[i][j]);
        }
    }

    c = trans_vec;    // <--- reassign here
}