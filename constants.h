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



// SNR Properties
const double Esn   = 1; // 1e51 erg : total energy released by SNR
const double Mej   = 1;    // Msun : total mass released by SNR in sun mass units 
const double xi_n  = 1;    // For solar abundances 
const double phi_c = 1;    // Actual thermal cond. / the Sptitzer (1962) value 
const double beta  = 2; 
const double xhi_cr = 0.1; // Efficiency of CRs acceleration 
const double xhi_0 = 2.026; 
const double gam = 2.2; // CRs injection energy power law index 
const double Emin = 0.1*GeV; // Minimum accelered CRs during the Sedov phase 
const double Emax = 2e5*GeV; // Maximum CRs energy