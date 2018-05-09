/**********************************************************************
 * DDCALC EXAMPLE PROGRAM (C++)
 * This program shows how to use the DDCalc module from C++, making
 * use of the interface defined in the DDCalc.hpp header file.
 * 
 * Run:
 *   ./DDCalc_exampleC 
 * 
 * 
 *       A. Scaffidi     U of Adelaide    2015    
 *       C. Savage       Nordita          2015
 *       P. Scott        Imperial Collge  2016
 *       F. Kahlhoefer   DESY		  2017
 *       S. Wild     	 DESY		  2017
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
  
  double mDM;
  
  // These three sets of indices refer to instance of the three types that
  // are the bedrock of DDCalc.  (Almost) every calculation needs to be
  // provided with an instance of each of these to do its job.  Passing the
  // index of one of them to DDCalc tells it which one in its internal cache
  // to use for the calculation requested. You can make have as many
  // different instances as you want, corresponding to e.g. different
  // detectors/analyses, WIMP models and DM halo models; the factory
  // funcions create the instances in DDCalc and return you the index of
  // the resulting object.
  int WIMP;
  int Halo;
  int XENON, LUX;  


  /* Initialise a DM Halo object to default values.*/
  Halo = DDCalc::InitHalo();
    
  /* Initialise a WIMP object to default values.*/
  WIMP = DDCalc::InitWIMP();


  /* Initialise detector objects */
  LUX      = DDCalc::LUX_2016_Init();


  /* Optionally set the Standard Halo Model parameters:
       rho     Local dark matter density [GeV/cm^3]
       vrot    Local disk rotation speed [km/s]
       v0      Maxwell-Boltzmann most probable speed [km/s]
       vesc    Galactic escape speed [km/s] */
  DDCalc::SetSHM(Halo, 0.3, 235.0, 235.0, 550.0);
 
  mDM = 100.0;
  DDCalc::SetWIMP_mfa(WIMP,mDM,1.0,-2.0,0.0,0.0);

  /* Do rate calculations using the specified WIMP and halo parameters.*/
  DDCalc::CalcRates(LUX,WIMP,Halo);
    
  /* Event quantities. */
  printf("%-20s  %e     \n", "Expected signal     ",DDCalc::Signal(LUX));
    
  // Clean up all the objects
  DDCalc::FreeAll();

} 


// END FILE ------------------------------------------------------------

 
