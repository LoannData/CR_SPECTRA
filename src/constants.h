#include <math.h>   //Basic mathematic functions 

const double pi   = 3.1472;
const double yr   = 3.154e+7;   // 1 yr in [s]
const double kyr  = 1e3*yr;     // 1 kyr in [s]
const double pc   = 3.086e18;   // 1 pc in [cm]
const double kpc  = 1.e3*pc;    // 1 kpc in [cm]
const double GeV  = 0.00160218; // 1 GeV in [erg] 
const double TeV  = 1e3*GeV;    // 1 TeV in [erg]
const double eV   = GeV*1e-9;   // 1 eV in [erg]
const double MeV  = 1e-3*GeV;
const double mp   = 1.6726219e-24; // Proton mass [g]
const double mn   = 1.6749286e-24; // Neutron mass [g]
const double mHI  = mp; // HI mass [g]
const double mHII  = mp; // HII mass [g]
const double mHeI = 2*mp + 2*mn; // HeI mass [g]
const double mCII = 6*mp + 6*mn; // CII mass [g]
const double mH2  = 2*mp;        // H2 mass [g]
const double e = 4.80326e-10; // e charge in [statC]
const double c = 29979245800.; // light celerity in [cm/s]
const double kbolz = 1.3807e-16; // Boltzmann constant in CGS
const double kms = 1e5; // 1 km/s in cm/s 

// Solver Options (1 : On, 0 : Off)
const int solver_PcrAdvection = 1; // Advective term of the CR Pressure
const int solver_PcrDiffusion = 1; // Diffusive term of the CR Pressure
const int solver_IpAdvection  = 1; // Advective term of the waves energy density
const int solver_Ipgrowth     = 1; // Self-generated turbulence 
const int solver_energy       = 1; // Terms of energy dependance, adiabatic losses
const int solver_damping      = 1; // Term of waves damping 
const int solver_waves_source = 1; // Term of waves source (background turbulence)

// Run & Output parameters
// (Note ! For more options, you can directly edit the ./src/out.h file)
const int nproc = 4;                    // Number of processors for the run 
const double t_data_out_min = 0.*kyr;   // Instant of the first output data 
const double t_data_out_max = 10.*kyr; // Instant of the last output data
const int number_out_data   = 500;     // Total number of output data
const int time_distrib_of_data = 0;     // Time distribution of output data (0 : linspace, 1 : log10-space)
const double log_first_data = 1.001;    // 
const int delta_log_output = 1000;      // Number of time-step between two LogOutput
const double Tmax = 10.1*kyr;           // Define the limit time of your simulation 



// SNR Properties
const double Esn      = 1;      // 1e51 erg : total energy released by SNR
const double Mej      = 1;      // Msun : total mass released by SNR in sun mass units 
const double xi_n     = 1;      // For solar abundances 
const double phi_c    = 1;      // Actual thermal cond. / the Sptitzer (1962) value 
const double beta     = 2; 
const double C06      = 1.;
const double xhi_cr   = 0.1;   // Efficiency of CRs acceleration 
const double xhi_0    = 2.026; 
const double gam      = 2.2;      // CRs injection energy power law index 
const double Emin     = 0.1*GeV;  // Minimum accelered CRs during the Sedov phase 
const double delta    = 3; 
//const double Emax     = 2e5*GeV;  // Maximum CRs energy
const double t_start_injection        = 1e-6*kyr; // Time start CRs injection function 
const double t_end_injection          = 2;   // [in tesc[E] units] Time end CRs injection function (number of tesc)
const double injection_function_width = 100. ; // Corresponds approximately to the width of the escape time divided 
const int injection_function_norm     = 100; // Constant in order to easily and rapidly normalize the injection function 