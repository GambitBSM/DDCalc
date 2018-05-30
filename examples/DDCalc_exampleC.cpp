/**********************************************************************
 * DDCALC EXAMPLE PROGRAM (C++)
 * This program shows how to use the DDCalc module from C++, making
 * use of the interface defined in the DDCalc.hpp header file.
 * 
 * For various types of DM-nucleon interactions, as well as for various 
 * direct detection experiments, it computes the expected signal rates,
 * any by comparing to the observed number of events, also the 
 * log(likelihood) value for each of the examples.
 * In order to convert those likelihood values into an upper bound on
 * e.g. the scattering cross section, please have a look at
 * DDCalc_exclusionC.cpp.
 * 
 * 
 * Run:
 *   ./DDCalc_exampleC 
 * 
 * 
 *       A. Scaffidi     U of Adelaide    2015    
 *       C. Savage       Nordita          2015
 *       P. Scott        Imperial Collge  2016
 *       F. Kahlhoefer   RWTH Aachen      2018
 *       S. Wild     	 DESY		  2017-2018
 *       ddcalc@projects.hepforge.org
 * 
 **********************************************************************/

// All the DDCalc routines used below are declared in (or included from)
// the DDCalc.hpp header file.

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include "DDCalc.hpp"


// MAIN PROGRAM --------------------------------------------------------

int main()
{  
  
  double mDM, sigmap_SI, sigman_SI, sigmap_SD, sigman_SD, DM_spin, fp, fn, ap, an;
  
  // These three sets of indices refer to instance of the three types that
  // are the bedrock of DDCalc.  (Almost) every calculation needs to be
  // provided with an instance of each of these to do its job.  Passing the
  // index of one of them to DDCalc tells it which one in its internal cache
  // to use for the calculation requested. You can create as many
  // different instances as you want, corresponding to e.g. different
  // detectors/analyses, WIMP models and DM halo models; the factory
  // funcions create the instances in DDCalc and return you the index of
  // the resulting object.
  int WIMP;
  int Halo;
  int Detector;  


  /* Initialise a DM Halo object to default values.*/
  Halo = DDCalc::InitHalo();
    
  /* Initialise a WIMP objects to default values.*/
  WIMP = DDCalc::InitWIMP();	

  /* Optionally set the Standard Halo Model parameters to values different from the default choices:
       rho     Local dark matter density [GeV/cm^3]
       vrot    Local disk rotation speed [km/s]
       v0      Maxwell-Boltzmann most probable speed [km/s]
       vesc    Galactic escape speed [km/s] */
  DDCalc::SetSHM(Halo, 0.3, 235.0, 235.0, 550.0);



  /* **************************************************************************************************************** */
  /* Example 1: XENON1T (2017) analysis, with standard SI/SD interactions specified by WIMP-nucleon cross sections.   */
  Detector = DDCalc::XENON1T_2017_Init();	// Initalize the XENON1T_2017 detector.
  mDM = 100.0;                           	// DM Mass in GeV.
  sigmap_SI = 4.0e-9;				// SI WIMP-proton cross section in pb.
  sigman_SI = -0.3e-9;				// SI WIMP-neutron cross section in pb.
						//   The negative value indicates that the corresponding WIMP-nucleon coupling is negative.
  sigmap_SD = 2.0e-5;				// SD WIMP-proton cross section in pb.
  sigman_SD = 8.0e-5;				// SD WIMP-neutron cross section in pb.
  DDCalc::SetWIMP_msigma(WIMP, mDM, sigmap_SI, sigman_SI, sigmap_SD, sigman_SD); 
  DDCalc::CalcRates(Detector,WIMP,Halo); 	// This performs the actual calculation of the rates.


  printf("********************************************************************************\n");
  printf("Example 1: Mixed SI and SD interactions at XENON1T (2017),\n");
  printf("           specified by the WIMP-nucleon cross sections.\n\n");
  printf("mDM       = %.5e GeV\n", mDM);
  printf("sigmap_SI = %.5e pb\n", sigmap_SI);
  printf("sigman_SI = %.5e pb\n", sigman_SI);
  printf("sigmap_SD = %.5e pb\n", sigmap_SD);
  printf("sigman_SD = %.5e pb\n\n", sigman_SD);
  printf("Expected number of signal events:     %.5e\n",DDCalc::Signal(Detector));
  printf("Observed number of events:            %d\n",DDCalc::Events(Detector));
  printf("Expected number of background events: %.5e\n",DDCalc::Background(Detector));
  printf("Log(likelihood):                      %.5e\n",DDCalc::LogLikelihood(Detector));
  printf("********************************************************************************\n\n");
  /* **************************************************************************************************************** */




  /* **************************************************************************************************************** */
  /* Example 2: PICO-60 (2017) analysis, with standard SI and SD interactions
                specified by the couplings fp, fn, ap, an.          
                fp and fn are the coefficients of the DM DM N N term in the Lagrangian, assuming a Majorana DM particle.
                The definition of ap and an follow the convention in Jungman&Kamionkowski, 'SUPERSYMMETRIC DARK MATTER'.
  */ 
  Detector = DDCalc::PICO_60_2017_Init();	// Initalize the PICO_60_2017 detector.
  mDM = 30.0;                           	// DM Mass in GeV.
  fp = 1.0e-8;					// SI WIMP-proton coupling fp, in units [GeV^(-2)].
  fn = 3.0e-8;					// SI WIMP-neutron coupling fn, in units [GeV^(-2)].
  ap = 1.0e-2;					// SD WIMP-proton coupling ap, unitless.
  an = 0.0;					// SD WIMP-neutron coupling an, unitless.
  DDCalc::SetWIMP_mfa(WIMP, mDM, fp, fn, ap, an); 
  DDCalc::CalcRates(Detector,WIMP,Halo); 	// This performs the actual calculation of the rates.

  printf("********************************************************************************\n");
  printf("Example 2: Mixed SI and SD interactions at PICO-60 (2017),\n");
  printf("           specified by the WIMP-nucleon couplings.\n\n");
  printf("mDM  = %.5e GeV\n", mDM);
  printf("fp   = %.5e GeV^(-2)\n", fp);
  printf("fn   = %.5e GeV^(-2)\n", fn);
  printf("ap   = %.5e\n", ap);
  printf("an   = %.5e\n\n", an);
  printf("Expected number of signal events:     %.5e\n",DDCalc::Signal(Detector));
  printf("Observed number of events:            %d\n",DDCalc::Events(Detector));
  printf("Expected number of background events: %.5e\n",DDCalc::Background(Detector));
  printf("Log(likelihood):                      %.5e\n",DDCalc::LogLikelihood(Detector));
  printf("********************************************************************************\n\n");
  /* **************************************************************************************************************** */




  /* **************************************************************************************************************** */
  /* Example 3: CRESST (2017) analysis, with a couple of non-relativistic scattering operators                        */
  /*            set to a non-zero value.  									      */
  /*            In particular, the analysis of this experiment involves several bins in energy.                       */

  Detector = DDCalc::CRESST_II_Init();				// Initalize the CRESST_II detector.
  mDM = 3.0;                           				// DM Mass in GeV.
  DM_spin = 0.5;						// DM Spin.
  DDCalc::SetWIMP_NREffectiveTheory(WIMP, mDM, DM_spin);
	// This defines a WIMP within the non-relativistic effective theory of DM-nucleon interactions, with a given mass in [GeV] and spin.
	// Initially, all the couplings corresponding to the various operators are set to zero.

  // The operator coefficients are set by DDCalc::SetNRCoefficient(WIMP, OpIndex, tau, value).
  //    OpIndex is an integer specifing the operator, ranging from 1 to 18 (excluding 2 and 16). Additional possible values are -1 and -4,
  //      corresponding to the momentum suppressed operators O_1 * (q^2/mp^2) and O_4 * (q^2/mp^2), respectively.
  //    tau is an integer, with 0 corresponding to the isoscalar part of the operator, and 1 to the isovector part.
  //    value specifies the operator coefficient corresponding to OpIndex and tau, in units [GeV^(-2)].
  DDCalc::SetNRCoefficient(WIMP, 11, 0, 5.0e-4);	// This sets the isoscalar part of O_11 to 5.0e-4 GeV^(-2).
  DDCalc::SetNRCoefficient(WIMP, 11, 1, -2.0e-4);       // This sets the isovector part of O_11 to -2.0e-4 GeV^(-2).
  DDCalc::CalcRates(Detector,WIMP,Halo); 	// This performs the actual calculation of the rates.


  printf("********************************************************************************\n");
  printf("Example 3: Non-standard scattering operators at CRESST-II.\n");
  printf("           This is also an example for an experiment with several bins.\n\n");
  printf("mDM       = %.5e GeV\n", mDM);
  printf("DM spin   = %.5e\n", DM_spin);
  printf("\t\tExpected signal\t\tObserved\tExpected background\n");
  printf("All bins\t%.5e\t\t%d\t\t%.5e\n", DDCalc::Signal(Detector), DDCalc::Events(Detector), DDCalc::Background(Detector));
  for( int i_bin = 1; i_bin <= DDCalc::Bins(Detector); i_bin = i_bin + 1 ) {
     printf("Bin %d\t\t%.5e\t\t%d\t\t%.5e\n", i_bin, DDCalc::BinSignal(Detector, i_bin), DDCalc::BinEvents(Detector, i_bin), DDCalc::BinBackground(Detector, i_bin));
     }
  printf("Log(likelihood):                      %.5e\n",DDCalc::LogLikelihood(Detector));
  printf("********************************************************************************\n\n");
  /* **************************************************************************************************************** */



  // Clean up all the objects
  DDCalc::FreeAll();

} 


// END FILE ------------------------------------------------------------

  
