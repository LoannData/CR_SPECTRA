#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;


void thetaDiffusionSolver(vector<vector<double> > &u, vector<vector<double> > &Pcr_new,   double dt, vector<double> X, int NE, vector<vector<double> > Ip, vector<vector<double> > Im, vector<vector<double> > Db, vector<vector<double> > &Pcr_background)
{
    double theta = 1.;
    int NX = X.size() - 1;

    vector<double> R(NX+1, 0.);
    vector<double> aa(NX+1, 0.);
    vector<double> bb(NX+1, 0.);
    vector<double> cc(NX+1, 0.);
    vector<double> loc_Pcr(NX+1, 0.);
    vector<double> loc_empty(NX+1, 0.);

    double u_c, u_l, u_r; 
    double dx, F1, alpha_m, alpha_p, an, bn, cn, rn;

    int ei, xi;

    //cout<<"Hello !"<<endl;

    //#pragma omp parallel num_threads(nproc)
    #pragma omp for schedule(static, int(double(NE/nproc))) private(R, aa, bb, cc, loc_Pcr, u_c, u_l, u_r, dx, F1, alpha_m, alpha_p, an, bn, cn, rn, ei, xi)
    for (ei = 0; ei < NE; ei++)
    {

        R       = loc_empty;
        aa      = loc_empty;
        bb      = loc_empty;
        cc      = loc_empty;
        loc_Pcr = loc_empty; 


        if (set_background == 1){
        for (xi = 0; xi < NX+1; xi++)
        {
            loc_Pcr[xi] = Pcr_background[xi][ei];
            R[xi] = Pcr_background[xi][ei];
        }}

        
        for (xi = 1; xi < NX; xi++)
        {

            dx = 0.5*(X[xi+1] - X[xi-1]);
            F1 = theta*dt/pow(dx, 2.);
            alpha_m = 0.5*(Db[xi][ei]/(Ip[xi][ei]+Im[xi][ei]) + Db[xi-1][ei]/(Ip[xi-1][ei]+Im[xi-1][ei])); 
            alpha_p = 0.5*(Db[xi][ei]/(Ip[xi][ei]+Im[xi][ei]) + Db[xi+1][ei]/(Ip[xi+1][ei]+Im[xi+1][ei]));

            an = -F1*alpha_m;
            bn = F1*alpha_p + F1*alpha_m + 1;
            cn = -F1*alpha_p; 

            rn = u[xi][ei] + (1-theta)*dt/pow(dx,2.)*(alpha_p*(u[xi+1][ei] - u[xi][ei]) - alpha_m*(u[xi][ei] - u[xi-1][ei]));

            aa[xi] = an;
            bb[xi] = bn;
            cc[xi] = cn;

            R[xi] = rn;

        }
        aa.erase(aa.begin()+0);
        cc.erase(cc.begin()+NX);

        // Cas xi = NX
        dx = X[NX] - X[NX-1]; 
        F1 = theta*dt/pow(dx,2);
        alpha_m = 0.5*(Db[NX][ei]/(Ip[NX][ei]+Im[NX][ei]) + Db[NX-1][ei]/(Ip[NX-1][ei]+Im[NX-1][ei]));
        alpha_p = 0.5*(Db[NX][ei]/(Ip[NX][ei]+Im[NX][ei]) + Db[0][ei]/(Ip[0][ei]+Im[0][ei]));//alpha_m;
        aa[NX] = -F1*alpha_m;
        bb[NX] = F1*alpha_p + F1*alpha_m +1; 
        R[NX]  = u[NX][ei] + (1-theta)*dt/pow(dx,2.)*(alpha_p*(u[0][ei] - u[NX][ei]) - alpha_m*(u[NX][ei] - u[NX-1][ei]));
        //R[NX] = u[NX-1][ei] - Pcr_background[NX-1][ei];

        // Cas xi = 0
        dx = X[1] - X[0];
        F1 = theta*dt/pow(dx,2);
        alpha_p = 0.5*(Db[0][ei]/(Ip[0][ei]+Im[0][ei]) + Db[1][ei]/(Ip[1][ei]+Im[1][ei]));
        alpha_m = 0.5*(Db[NX][ei]/(Ip[NX][ei]+Im[NX][ei]) + Db[NX-1][ei]/(Ip[NX-1][ei]+Im[NX-1][ei]));//alpha_p;
        cc[0] = -F1*alpha_p;
        bb[0] = F1*alpha_p + F1*alpha_m +1; 
        R[0]  = u[0][ei] + (1-theta)*dt/pow(dx,2.)*(alpha_p*(u[1][ei] - u[0][ei]) - alpha_m*(u[0][ei] - u[NX][ei]));

        loc_Pcr = TDMA(aa, bb, cc, R);

        if (set_background == 1){
        for (xi = 0; xi < NX+1; xi++)
        {
            Pcr_new[xi][ei] = max(loc_Pcr[xi], Pcr_background[0][ei]);
        }}
        if (set_background == 0){
        for (xi = 0; xi < NX+1; xi++)
        {
            Pcr_new[xi][ei] = loc_Pcr[xi]; 
        }}
    }
}




void advectionSolverX(vector<vector<double> > &u_old, vector<vector<double> > &u_new, double dt, vector<double> X, int NE, vector<vector<double> > V, int sign, int border)
{
    // Upwind advection scheme
    // CF : https://en.wikipedia.org/wiki/Upwind_scheme
    // Note : Problem with 2nd order scheme ! 
    int NX = X.size()-1;
    int order = 1;
    double ddx, ux_p, ux_m, a_p, a_m, V_loc;
    
    //#pragma omp parallel num_threads(nproc)
    #pragma omp for schedule(static, int(double(NE/nproc))) private(V_loc, ddx, ux_p, ux_m, a_p, a_m)
    for (int ei = 0; ei < NE; ei++)
    {   
        if (order == 1)
        {
            for (int xi = 1; xi < NX; xi++)
            {
                ddx = (X[xi+1] - X[xi-1])*0.5;
                V_loc = V[xi][ei];
                if (V_loc > 0){if (sign == -1){V_loc = - V_loc;}}
                if (V_loc < 0){if (sign ==  1){V_loc = - V_loc;}}
                a_p = max(V_loc, 0.);
                a_m = min(V_loc, 0.);
                ux_p = (u_old[xi+1][ei] - u_old[xi][ei])/ddx;
                ux_m = (u_old[xi][ei] - u_old[xi-1][ei])/ddx;
                u_new[xi][ei] = u_old[xi][ei] - dt*(a_p*ux_m + a_m*ux_p);
            }

            if (border == 0){
            // Absorbing boundaries
            u_new[0][ei]  = u_old[0][ei]; // Cas xi = 0
            u_new[NX][ei] = u_old[NX][ei]; // Cas xi = NX
            }

            if (border == 1){
            // Looped boundaries 
            // Case xi = NX
            ddx = (X[NX] - X[NX-1]);
            V_loc = V[NX][ei];
            if (V_loc > 0){if (sign == -1){V_loc = - V_loc;}}
            if (V_loc < 0){if (sign ==  1){V_loc = - V_loc;}}
            a_p = max(V_loc, 0.);
            a_m = min(V_loc, 0.);
            ux_p = (u_old[0][ei] - u_old[NX][ei])/ddx;
            ux_m = (u_old[NX][ei] - u_old[NX-1][ei])/ddx;
            u_new[NX][ei] = u_old[NX][ei] - dt*(a_p*ux_m + a_m*ux_p);
            
            // Case xi = 0
            ddx = (X[1] - X[0]);
            V_loc = V[0][ei];
            if (V_loc > 0){if (sign == -1){V_loc = - V_loc;}}
            if (V_loc < 0){if (sign ==  1){V_loc = - V_loc;}}
            a_p = max(V_loc, 0.);
            a_m = min(V_loc, 0.);
            ux_p = (u_old[1][ei] - u_old[0][ei])/ddx;
            ux_m = (u_old[0][ei] - u_old[NX][ei])/ddx;
            u_new[0][ei] = u_old[0][ei] - dt*(a_p*ux_m + a_m*ux_p);
            }
            
        }
    }
}

void advectionSolverE(vector<vector<double> > &u_old, vector<vector<double> > &u_new, double dt, vector<double> E, int NX, vector<vector<double> > V, vector<vector<double> > u_background)
{
    // Note : No 2nd order scheme ! 
    // Note : For the boundaries, we use absorbing ones. If energy leaves from the borders, it disapear of ever ...
    // We need to find a good method for the boundaries ...  
    int NE = E.size()-1;
    double dde, ux_p, ux_m, a_p, a_m, ad0, delta;
    int xi, ei;

    //#pragma omp parallel num_threads(nproc)
    #pragma omp for schedule(static, int(double(NX/nproc))) private(xi, ei, dde, ux_p, ux_m, a_p, a_m, ad0, delta)
    for (xi = 0; xi < NX; xi++)
    {
            for (ei = 1; ei < NE; ei++)
            {
                
                dde = (E[ei+1] - E[ei-1])/2.;
                a_p = max(V[xi][ei], 0.);
                a_m = min(V[xi][ei], 0.);
                ux_p = (u_old[xi][ei+1] - u_old[xi][ei])/dde;
                ux_m = (u_old[xi][ei]   - u_old[xi][ei-1])/dde;
                u_new[xi][ei] = u_old[xi][ei] - dt*(a_p*ux_m + a_m*ux_p);
                if (ei == 1){ad0 = - dt*(a_p*ux_m + a_m*ux_p);}
                //cout<<dt<<"  "<<a_p<<" "<<a_m<<" "<<ux_p<<" "<<ux_m<<" "<<- dt*(a_p*ux_m + a_m*ux_p)<<endl;
                if (set_background == 1){
                    if (u_new[xi][ei] < u_background[xi][ei]){u_new[xi][ei] = u_background[xi][ei];}}
            }
            //u_new[xi][0]  = u_old[xi][0] + ad0; // Cas ei = 0
            //u_new[xi][NE] = u_old[xi][NE] - dt*(a_p*ux_m + a_m*ux_p); // Cas ei = NE

            // Border conditions, extrapolation 
            // Case where ei = 0
            u_new[xi][0] = pow(u_new[xi][1],2.)/u_new[xi][2];

            // Case where ei = NE
            u_new[xi][NE] = pow(u_new[xi][NE-1],2.)/u_new[xi][NE-2];
    }
}


void advectionSolverE1(vector<vector<double> > &u_old, vector<vector<double> > &u_new, double dt, vector<double> E, vector<double> BB, vector <double> EE, int NX, vector<vector<double> > u_background)
{
    // Note : No 2nd order scheme ! 
    // Note : For the boundaries, we use absorbing ones. If energy leaves from the borders, it disapear of ever ...
    // We need to find a good method for the boundaries ...  
    int NE = E.size()-1;
    double dde, ux_p, ux_m, a_p, a_m, ad0, V_loc;
    int xi, ei;

    //#pragma omp parallel num_threads(nproc)
    #pragma omp for schedule(static, int(double(NX/nproc))) private(xi, ei, dde, ux_p, ux_m, a_p, a_m, ad0, V_loc)
    for (xi = 0; xi < NX; xi++)
    {
            for (ei = 1; ei < NE; ei++)
            {
                V_loc = - 2*c*sig_T*pow(2*BB[xi], 2)*EE[ei]/(4*pi*me*me*pow(c, 4))/log(10.); 
                
                dde = (E[ei+1] - E[ei-1])/2.;
                a_p = max(V_loc, 0.);
                a_m = min(V_loc, 0.);
                ux_p = (u_old[xi][ei+1] - u_old[xi][ei])/dde;
                ux_m = (u_old[xi][ei]   - u_old[xi][ei-1])/dde;
                u_new[xi][ei] = u_old[xi][ei] - dt*(a_p*ux_m + a_m*ux_p);

                //if (dt*(a_p*ux_m + a_m*ux_p) < 0.){
                //cout<<"dde = "<<dde<<", dt = "<<dt<<endl;
                //cout<<"E_ei+1 = "<<u_old[xi][ei+1]<<", E_ei = "<<u_old[xi][ei]<<", E_ei-1 = "<<u_old[xi][ei-1]<<endl;
                //cout<<"V_loc = "<<V_loc<<", a_p = "<<a_p<<", a_m = "<<a_m<<endl;
                //cout<<"ux_p = "<<ux_p<<", ux_m = "<<ux_m<<endl;
                //cout<<"a_p*ux_m + a_m*ux_p = "<<a_p*ux_m + a_m*ux_p<<", dt*(a_p*ux_m + a_m*ux_p) = "<<dt*(a_p*ux_m + a_m*ux_p)<<endl;
                //cin.get();}
                if (ei == 1){ad0 = - dt*(a_p*ux_m + a_m*ux_p);}
                //cout<<dt<<"  "<<a_p<<" "<<a_m<<" "<<ux_p<<" "<<ux_m<<" "<<- dt*(a_p*ux_m + a_m*ux_p)<<endl;
                if (set_background == 1){
                    if (u_new[xi][ei] < u_background[xi][ei]){u_new[xi][ei] = u_background[xi][ei];}}
                //else {if (u_new[xi][ei] < 1e-60){u_new[xi][ei] = 1e-60;}}
            }
            //u_new[xi][0]  = u_old[xi][0] + ad0; // Cas ei = 0
            //u_new[xi][NE] = u_old[xi][NE] - dt*(a_p*ux_m + a_m*ux_p); // Cas ei = NE

            // Border conditions, extrapolation 
            // Case where ei = 0
            u_new[xi][0] = pow(u_new[xi][1],2.)/u_new[xi][2];

            // Case where ei = NE
            u_new[xi][NE] = pow(u_new[xi][NE-1],2.)/u_new[xi][NE-2];
    }
}


void advectionSolverE2(vector<vector<double> > &u_old, vector<vector<double> > &u_new, double dt, vector<double> E, int NX, vector<double> BB, vector <double> EE, vector<vector<double> > u_background)
{
    // Note : No 2nd order scheme ! 
    // Note : For the boundaries, we use absorbing ones. If energy leaves from the borders, it disapear of ever ...
    // We need to find a good method for the boundaries ...  
    int NE = E.size()-1;
    double dde, ux_p, ux_m, a_p, a_m, ad0;
    //double C = c*sig_T/(4*log(10.));
    //double C_0 = c*sig_T/(4.*me*pow(c, 2.));
    double C2, C3;
    int xi, ei;

    //#pragma omp parallel num_threads(nproc)
    #pragma omp for schedule(static, int(double(NX/nproc))) private(xi ,ei, C3)
    for (xi = 0; xi < NX; xi++)
    {
            for (ei = 0; ei < NE+1; ei++)
            {
                C3 = 2*c*sig_T*pow(2*BB[xi], 2)*EE[ei]/(4*pi*me*me*pow(c, 4));
                if (source_terms_exact == 0)
                {
                    if (set_background == 1){u_new[xi][ei] = max(u_old[xi][ei]*(1 - C3*dt), u_background[xi][ei]);}
                    else {u_new[xi][ei] = u_old[xi][ei]*(1 - C3*dt);} 
                }
                if (source_terms_exact == 1)
                {
                    if (set_background == 1){u_new[xi][ei] = max(u_old[xi][ei]*exp(-C3*dt), u_background[xi][ei]);}
                    else {u_new[xi][ei] = u_old[xi][ei]*exp(-C3*dt);}
                }
            }
    }
}



void sourceSolver(vector<vector<double> > &u_old, vector<vector<double> > &u_new, double dt, vector<vector<double> > source, double factor)
    {
        int NX = u_old.size();
        int NE = u_old[0].size();;
        double sum;
        int ei, xi;

        //#pragma omp parallel num_threads(nproc)
        #pragma omp for schedule(static, int(double(NX/nproc))) private(sum, ei, xi)
        for (xi = 0; xi < NX; xi++)
        {
            for (ei = 0; ei < NE; ei++)
            {
                if (source_terms_exact == 0)
                {
                    sum = factor*source[xi][ei]*u_old[xi][ei];
                    u_new[xi][ei] = u_old[xi][ei] + sum*dt;
                }
                if (source_terms_exact == 1)
                {
                    u_new[xi][ei] = u_old[xi][ei]*exp(factor*source[xi][ei]*dt);
                }
            }
        }
    }

/*void sourceGrowthDampRateSolver(vector<vector<double> > &u_old, vector<vector<double> > &u_new, vector<vector<double> > v_old, vector<vector<double> > source, vector<vector<double> > background,  vector<double> X, double dt, vector<vector<double> > V, vector<double> B, int factor)
    {
        int NX = u_old.size();
        int NE = u_old[0].size();
        int xi, ei;
        double sum;
        double dudx; 
        double w0;
        double u_max = 1e4;
        //double tau_sat = - log(0.1)/0.01;

        //#pragma omp parallel num_threads(nproc)
        #pragma omp for schedule(static, int(double(NX/nproc))) private(sum, dudx, w0, xi)
        for (xi = 0; xi < NX; xi++)
        {
            w0 = B[xi]*B[xi]/(8*pi);
            for (ei = 0; ei < NE; ei++)
            {
                    dudx = 0.;
                    if (xi > 0 and xi < NX-1){dudx = (v_old[xi+1][ei]-v_old[xi-1][ei])/(X[xi+1] - X[xi-1]);}
                    if (factor ==  1){if (dudx > 0.){dudx = 0.;}} // Condition Foward waves
                    if (factor == -1){if (dudx < 0.){dudx = 0.;}} // Condition Backward waves
                    //sum = abs(0.5*V[xi][ei]*dudx/w0)*(log10(u_old[xi][ei]+1)/u_old[xi][ei]) - abs(source[xi][ei]*(u_old[xi][ei] - background[xi][ei]));///;
                    sum = abs(0.5*V[xi][ei]*dudx/w0)*(exp(-u_old[xi][ei]*ttau_sat)/u_old[xi][ei]) - abs(source[xi][ei]*(u_old[xi][ei] - background[xi][ei]));///;
                    u_new[xi][ei] = u_old[xi][ei] + sum*dt;
                    if (u_new[xi][ei] < background[xi][ei]) {u_new[xi][ei] = background[xi][ei];}
                    if (u_new[xi][ei] > u_max) {u_new[xi][ei] = u_max;}
            }
        }
    }*/

void sourceGrowthDampRateSolver(vector<vector<double> > &u_old, vector<vector<double> > &u_new, vector<vector<double> > v_old, vector<vector<double> > source, vector<vector<double> > background,  vector<double> X, double dt, vector<vector<double> > V, vector<double> B, int factor)
    {
        int NX = u_old.size();
        int NE = u_old[0].size();
        int xi, ei;
        double sum;
        double dudx; 
        double w0;
        double u_max = 1e4;
        //double tau_sat = - log(0.1)/0.01;

        //#pragma omp parallel num_threads(nproc)
        #pragma omp for schedule(static, int(double(NX/nproc))) private(sum, dudx, w0, xi)
        for (xi = 0; xi < NX; xi++)
        {
            w0 = B[xi]*B[xi]/(8*pi);
            for (ei = 0; ei < NE; ei++)
            {
                    dudx = 0.;
                    if (xi > 0 and xi < NX-1){dudx = (v_old[xi+1][ei]-v_old[xi-1][ei])/(X[xi+1] - X[xi-1]);}
                    if (factor ==  1){if (dudx > 0.){dudx = 0.;}} // Condition Foward waves
                    if (factor == -1){if (dudx < 0.){dudx = 0.;}} // Condition Backward waves


                    sum = abs(0.5*V[xi][ei]*dudx/w0)*exp(-u_old[xi][ei]*ttau_sat) - abs(source[xi][ei])*(u_old[xi][ei] - background[xi][ei]);  


                    //sum = abs(0.5*V[xi][ei]*dudx/w0)*(exp(-u_old[xi][ei]*ttau_sat)/u_old[xi][ei]) - abs(source[xi][ei]*(u_old[xi][ei] - background[xi][ei]));///;

                    u_new[xi][ei] = u_old[xi][ei] + sum*dt;
                    if (u_new[xi][ei] < background[xi][ei]) {u_new[xi][ei] = background[xi][ei];}
                    if (u_new[xi][ei] > u_max) {u_new[xi][ei] = u_max;}
            }
        }
    }




void CRsInjectionSourceSolver(vector<vector<double> > &u_old, vector<vector<double> > &u_new, double dt, vector<double> Pcr_ini, vector<double> Finj_temp, vector<double> vec_theta)
    {
        int NX = u_old.size();
        int NE = u_old[0].size();

        //#pragma omp parallel num_threads(nproc)
        #pragma omp for schedule(static, int(double(NX/nproc))) //private(NE)
        for (int xi = 0; xi < NX; xi++)
        {   
            //NE = u_old[xi].size();
            for (int ei = 0; ei < NE; ei++)
            {
                u_new[xi][ei] = u_old[xi][ei] + dt*Pcr_ini[ei]*Finj_temp[ei]*vec_theta[xi]; 
            }
        }
    }




/*
void dilute_solver(vector<vector<double> > &u_old, vector<vector<double> > &u_new, vector<vector<double> > u_background, double r_new, double r_old)
{
    int NX = u_old.size();
    int NE = u_old[0].size();
        #pragma omp for schedule(static, int(double(NX/nproc))) 
        for (int xi = 0; xi < NX; xi++)
        {   
            //NE = u_old[xi].size();
            for (int ei = 0; ei < NE; ei++)
            {
                u_new[xi][ei] = max(u_old[xi][ei]*pow(r_old/r_new,2.), u_background[xi][ei]);
            }
        }
}
*/


void dilute_solver(vector<vector<double> > &u_old, vector<vector<double> > &u_new, vector<vector<double> > u_background, double r_new, double r_old, vector<double> nn, vector<double> ni, vector<double> mn, vector<double> mi, vector<double> B, vector<double> T)
{
    int NX = u_old.size();
    int NE = u_old[0].size();
    double rho_g, Cma, Pc, Pg, Pm;
    double n_tot; 
    double gamma_m = 2; 
    double gamma_g = 5./3; 
    double gamma_c = 4./3; 
        #pragma omp for schedule(static, int(double(NX/nproc))) private(rho_g, Cma, Pc, Pg, Pm, n_tot)
        for (int xi = 0; xi < NX; xi++)
        {   
            rho_g = mn[xi]*nn[xi] + mi[xi]*ni[xi]; 
            n_tot = nn[xi] + ni[xi]; 
            Pg = n_tot*kbolz*T[xi]; 
            Pm = pow(B[xi],2)/(8*pi); 

            

            //cout<<"Cma = "<<Cma<<", Pg = "<<Pg<<", Pm = "<<Pm<<", rho_g = "<<rho_g<<" T = "<<T[xi]<<endl;

            //NE = u_old[xi].size();
            for (int ei = 0; ei < NE; ei++)
            {
                Cma = pow((gamma_g*Pg + gamma_c*u_background[xi][ei] + gamma_m*Pm)/rho_g, 0.5); 
                //u_new[xi][ei] = max(u_old[xi][ei]*pow(r_old/r_new,2.), u_background[xi][ei]);

                u_new[xi][ei] = u_old[xi][ei]/(1 + (u_old[xi][ei] - u_background[xi][ei])/(rho_g*pow(Cma,2.))); 

                //cout<<"r = "<<u_new[xi][ei]/u_old[xi][ei]<<endl;
                //cout<<"r = "<<(u_old[xi][ei] - u_background[xi][0])/(rho_g*pow(Cma,2.))<<endl;
            }
        }
}






void NotMove(vector<vector<double> > u, vector<vector<double> > u_new)
{
    for (int ei = 0; ei < u.size(); ei++)
    {
        for (int xi = 0; xi < u[ei].size(); xi++)
        {
            u_new[xi][ei] = u[xi][ei];
        }
    }
}

//=============================================//
// OTHER SOLVERS OPTIONNAL                     //
//=============================================//
void electron_source(vector<vector<double> > &u_old, vector<vector<double> > &u_new, vector<double> E, double dt, int NE, int NX, vector<vector<double> > u_background)
{ 
    double Q;
    int xi, ei;

    //#pragma omp parallel num_threads(nproc)
    #pragma omp for schedule(static, int(double(NX/nproc))) private(xi ,ei, Q)
    for (xi = 0; xi < NX; xi++)
    {
            for (ei = 0; ei < NE; ei++)
            {
                Q = 1e-14*pow(E[ei]/(10.*GeV), -3.1);


                if (set_background == 1){
                u_new[xi][ei] = max(u_old[xi][ei] + Q*dt, u_background[xi][ei]);}
                else {u_new[xi][ei] = u_old[xi][ei] + Q*dt;} 
                //cout<<dt<<" "<<BB[xi]<<" "<<E[ei]<<" "<<sig_T<<endl;
            }
    }
}
