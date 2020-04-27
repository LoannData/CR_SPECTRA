#include <math.h>   //Basic mathematic functions 

const double pi   = 3.1472;
const double yr   = 3.154e+7;   // 1 yr in [s]
const double kyr  = 1e3*yr;     // 1 kyr in [s]
const double km   = 1e5;        // 1 km in [cm]
const double pc   = 3.086e18;   // 1 pc in [cm]
const double kpc  = 1.e3*pc;    // 1 kpc in [cm]
const double GeV  = 0.00160218; // 1 GeV in [erg] 
const double TeV  = 1e3*GeV;    // 1 TeV in [erg]
const double eV   = GeV*1e-9;   // 1 eV in [erg]
const double MeV  = 1e-3*GeV;
const double mp   = 1.6726219e-24; // Proton mass [g]
const double mn   = 1.6749286e-24; // Neutron mass [g]
const double me   = 9.1095e-28;    // Electron mass [g]
const double mHI  = mp; // HI mass [g]
const double mHII  = mp; // HII mass [g]
const double mHeI = 2*mp + 2*mn; // HeI mass [g]
const double mCII = 6*mp + 6*mn; // CII mass [g]
const double mH2  = 2*mp;        // H2 mass [g]
const double e = 4.80326e-10; // e charge in [statC]
const double c = 29979245800.; // light celerity in [cm/s]
const double kbolz = 1.3807e-16; // Boltzmann constant in CGS
const double kms = 1e5; // 1 km/s in cm/s
const double sig_T = 6.65e-25; // [cm^2] Thomson cross section (see Schlickeiser (2002) p.80) 

// Solver Options (1 : On, 0 : Off)
const int solver_PcrAdvection  = 1; // Advective term of the CR Pressure (the classical one -> V_A*grad ...)
const int solver_PcrDiffusion  = 1; // Diffusive term of the CR Pressure
const int solver_PcrAdvection2 = 1; // Explicit Advection solver for Pcr by the energy derivative of Alfvén velocity.
const int solver_PcrAdvectionE = 1; // Explicit Advection solver for Pcr in energy cdVAdX
const int solver_PcrSource1    = 1; // Source term effect due to the dependance of the Alfvén velocity to the space
const int solver_PcrSource2    = 1; // Source term effect due to the CR injection from the source in the system  

const int solver_PeAdvection   = 1; // Advective term of the e- Pressure (the classical one -> V_A*grad ...)
const int solver_PeDiffusion   = 1; // Diffusive term of the e- Pressure
const int solver_PeAdvection2  = 1; // Explicit Advection solver for e- by the energy derivative of Alfvén velocity.
const int solver_PeAdvectionE  = 1; // Explicit Advection solver for e- in energy cdVAdX

const int solver_PeAdvectionE1 = 1; // Sychrotron radiations of e- (looses of energy) - Advection term 
const int solver_PeAdvectionE2 = 1; // Synchrotron radiations of e- (looses of energy) - source term

const int solver_PeSource1     = 1; // Source term effect due to the dependance of the Alfvén velocity to the space (for e-)
const int solver_PeSource2     = 1; // Source term effect due to the e- injection from the source in the system 

const int solver_IpAdvection   = 1; // Advective term of the foward waves 
const int solver_ImAdvection   = 1; // Advective term of the backward waves
const int solver_IpSource1     = 1; // Source term effect applied on foward waves due to the dependance of the Alfvén velocity to the space
const int solver_ImSource1     = 1; // Source term effect applied on backward waves due to the dependance of the Alfvén velocity to the space
const int solver_IpDampGrowth  = 1; // Source term effect due to production of self-turbulence - damping applied on foward waves 
const int solver_ImDampGrowth  = 1; // Source term effect due to production of self-turbulence - damping applied on backward waves 

const int solver_Dilution      = 0; // Time dilution term according to the SNR shock evolution in the flux tube approx. ie. R_sh^2(t0) P(t0) = R_sh^2(t1) P(t1) if conserved energy 
                                    // Still experimental !!! 

// Other parameters 
const int set_background = 1; // If need to perform some test without background conditions (1 : On, 0 : off)
const double step_implicit = 500.*yr; // Time step of the simulation if only implicit solvers are used and maximum time step value.
const int source_terms_exact = 1; // The way to solve the source terms : 1 -> Exact solutions, 0 -> 1st order numerical solution  

// Run & Output parameters
// (Note ! For more options, you can directly edit the ./src/out.h file)
const int nproc = 1;                    // Number of processors for the run 

const int output_freq = 0;               // Model of output frequency (0 : n output between t_start and t_end, 1 : 1 output each n timestep) n = number_out_data
const double t_data_out_min = 0.*kyr;   // Instant of the first output data 
const double t_data_out_max = 200.*kyr; //200.*kyr; // Instant of the last output data
const int number_out_data   = 200;     // Total number of output data
const int time_distrib_of_data = 0;     // Time distribution of output data (0 : linspace, 1 : log10-space)
const double log_first_data = 1.001;    // 
const int delta_log_output = 100;      // Number of time-step between two LogOutput
const double Tmax = 200.1*kyr;//200.1*kyr;           // Define the limit time of your simulation 




// Growth waves saturation rate
const double ttau_sat = - log(0.1)/0.1; // Has the form - log(a)/b where b : characteristic max value, a : suppression factor after b, ttau_sat = 0 -> Linear growth 

// SNR Properties (Model from Cioffi et al. 2012 & ...)
//const double snr_position_x = 500*pc; // Position of the center of the source on the grid
const double Esn      = 1;      // 1e51 erg : total energy released by SNR
const double Mej      = 1;      // Msun : total mass released by SNR in sun mass units 
const double xi_n     = 1;      // For solar abundances 
const double phi_c    = 1;      // Actual thermal cond. / the Sptitzer (1962) value 
const double bbeta    = 2; 
const double C06      = 1.;
const double xhi_cr   = 0.1;    // Efficiency of CRs acceleration 
const double xhi_0    = 2.026; 
const double gam      = 2.2;      // CRs injection energy power law index 
const double Emin     = 0.1*GeV;  // Minimum accelered CRs during the Sedov phase 
const double delta    = 2;       // From Celli et al. (2019) - see Brahimi et al. (2020)
//const double Emax     = 2e5*GeV;  // Maximum CRs energy
const double t_start_injection        = 1e-6*kyr; // Time start CRs injection function 
const double t_end_injection          = 2;   // [in tesc[E] units] Time end CRs injection function (number of tesc)
const double injection_function_width = 5. ; // Corresponds approximately to the width of the escape time divided 
                                             // (take care : too high value may affect the calculation and create noise)
                                             // Max value -> value such as the width of Finj >> dt, max recommended value = 10 for NX,NE = 10, 7 
const int injection_function_norm     = 100; // Constant in order to easily and rapidly normalize the injection function 
const double r_snr_thickness          = 100; // = R_SNR(t)/by the value you chose, it allows to smooth the injection shape of CRs  
const double electron_injection_rate  = 1e-2; // Corresponds to the energy injected in the electron spectrum compared to the energy injected in the proton spectrum

const int tesc_model  = 1; // CR escape time model (1 : All CRs escape at the begining of the radiative phase, 2 : If v_sh < 110 km/s, all CRs escape)


// Electrons escape model (from Ohira et al. 2018)

// Model 1 : eta_g = eta_g,free -> 1
// Model 2 : B^2 = u_sh^3       -> 2
// Model 3 : B^2 = u_sh^2       -> 3 (Best model)
const int oh_model = 3;

const double eta_gfree = 1;
const double eta_acc = 10;
const double Eknee = 1e15*eV;
const double alpha = 2.6; 
