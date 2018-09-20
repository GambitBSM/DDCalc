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
  
  int WIMP1, WIMP2;
  int Halo;
  int Detector1, Detector2;  
  int opi;
  double coeff;


  Halo = DDCalc::InitHalo();
  WIMP1 = DDCalc::InitWIMP();	
  WIMP2 = DDCalc::InitWIMP();	

  Detector1 = DDCalc::XENON1T_2018_Init();			
  Detector2 = DDCalc::XENON1T_2018_Init();			

  mDM = 50.0;                          				
  DM_spin = 0.5;						
  DDCalc::SetWIMP_NREffectiveTheory(WIMP1, mDM, DM_spin);
  DDCalc::SetWIMP_NREFT_CPT(WIMP2, mDM, DM_spin);

  opi = 13;
  coeff = 5.0e-4;

  DDCalc::SetNRCoefficient(WIMP1, opi, 0, coeff);	
  DDCalc::CalcRates(Detector1,WIMP1,Halo); 	

  DDCalc::SetNRCoefficient(WIMP2, opi, 0, coeff/2);	
  DDCalc::SetNRCoefficient(WIMP2, opi, 1, coeff/2);	
  DDCalc::CalcRates(Detector2,WIMP2,Halo); 	

  printf("Signal: \t%.5e\t%.5e\n", DDCalc::Signal(Detector1), DDCalc::Signal(Detector2));

  // Clean up all the objects
  DDCalc::FreeAll();

} 


// END FILE ------------------------------------------------------------

  
