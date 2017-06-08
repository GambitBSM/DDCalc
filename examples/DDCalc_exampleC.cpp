/**********************************************************************
 * DDCALC EXAMPLE PROGRAM (C++)
 * This program shows how to use the DDCalc module from C++, making
 * use of the interface defined in the DDCalc.hpp header file.
 * 
 * Run:
 *   ./DDCalc_exampleC [--mG|--mfa|--msigma]
 * where the optional flag specifies the form in which the WIMP-nucleon
 * couplings will be provided (default: --msigma).
 * 
 * 
 *       A. Scaffidi     U of Adelaide    2015    
 *       C. Savage       Nordita          2015
 *       P. Scott        Imperial Collge  2016
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


// CONSTANTS -----------------------------------------------------------

// These constants will be used to specify the type of input parameters.
//const int TYPE_MG     = 1;  // Four-fermion effective couplings G
const int TYPE_MFA    = 2;  // Effective couplings f (SI), a (SD)
//const int TYPE_MSIGMA = 3;  // WIMP-nucleon cross-sections


// UTILITY FUNCTIONS DECLARATIONS --------------------------------------

/* Provide prototypes for utility functions used by this example
   program.  The function definitions are given later in this file,
   after the main() routine. */

// Print description of input
void WriteDescription(const int type);

// Get WIMP mass and couplings
bool GetWIMPParams(const int type, double& M, double& xpSI, double& xnSI,
                   double& xpSD, double& xnSD);


// MAIN PROGRAM --------------------------------------------------------

int main(int argc, char* argv[])
{  
  
  int type;
  double M,xpSI,xnSI,xpSD,xnSD,GpSI,GnSI,GpSD,GnSD,fp,fn,ap,an,
         sigmapSI,sigmanSI,sigmapSD,sigmanSD;
  
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
  int MyDetector, XENON, LUX, SCDMS, SIMPLE;  

  // Parse command line options
  // Notably, determining how WIMP parameters will be specified.
  // Default command line option (no argument) will give type = TYPE_MFA.
  type = TYPE_MFA;
  for (int i=1; i<argc; i++)
  {
    if (std::string(argv[i]) == "--mfa")
      type = TYPE_MFA;
    else if (std::string(argv[i]) == "--help")
    {
      std::cout << "Usage:" << std::endl;
      std::cout << "  ./DDCalc_exampleC [--mfa]" << std::endl;
      std::cout << "where the optional flag specifies the form in which the WIMP-" << std::endl;
      std::cout << "nucleon couplings will be provided (default: --mfa)." << std::endl;
      exit(0);
    } else
    {
      std::cout << "WARNING:  Ignoring unknown argument '" << argv[i] << "'." << std::endl;
    }
  }
  
  /* Note that we never have to initialise DDCalc as a whole, we just
     have to create Detectors, WIMPs and Halos, then manipulate their
     parameters and hand them back to DDCalc to do calculations on. */

  /* Initialise a DM Halo object to default values.  See below for how to
     modify these values. */
  Halo = DDCalc::InitHalo();
    
  /* Initialise a WIMP object to default values.  Actually, this isn't
     necessary here, as we set the WIMP properties from the commandline
     later -- but here's how you would make a default version if needed: */
  WIMP = DDCalc::InitWIMP();

  /* As we are responsible adults, we also choose our own detectors below.
    Here is what you'd do if you wanted to just rely on the default
    (currently LUX 2013): */
  MyDetector = DDCalc::InitDetector(true);

  /* Explicitly create detector objects for all the experiments to be
     used (set up isotopes, efficiencies, array sizing, etc.)  The   
     argument indicates if extra sub-interval calculations should
     be performed.  Those calculations are required for maximum gap
     analyses, but are unnecessary for calculating total rates and
     likelihoods.  If .FALSE. is given, a no-background-subtraction
     p-value can still be calculated, but a Poisson is used instead
     of the maximum gap.  We show some maximum gap results below, so
     we must use true here (the flag is ignored for experiments
     that do not have the event energies necessary for a maximum gap
     analysis). */
  XENON    = DDCalc::XENON100_2012_Init(true);
  LUX      = DDCalc::LUX_2013_Init(true);
  SCDMS    = DDCalc::SuperCDMS_2014_Init(true);
  SIMPLE   = DDCalc::SIMPLE_2014_Init(true);

  /* Can optionally specify a minimum recoil energy to be included in
     the rate calculations [keV].  Note the efficiency curves already
     account for detector and analysis thresholds regardless of this
     setting, so setting this to 0 keV (the default behavior when
     initialization is performed) does not imply that very low energy
     recoils actually contribute to the signal.
     EXAMPLE: Uncomment to set a minimum recoil energy of 3 keV for LUX: */
  //DDCalc::SetEmin(LUX,3.0)

  /* Optionally set the Standard Halo Model parameters:
       rho     Local dark matter density [GeV/cm^3]
       vrot    Local disk rotation speed [km/s]
       v0      Maxwell-Boltzmann most probable speed [km/s]
       vesc    Galactic escape speed [km/s]
     This example uses the default values (and is thus optional). */
  //DDCalc::SetSHM(0.4, 235.0, 235.0, 550.0)
  
  // INPUT LOOP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Loop over input to this example program.
  // GetWIMPParams is defined below.
  while(GetWIMPParams(type,M,xpSI,xnSI,xpSD,xnSD))
  {

    std::cout << std::endl;
    
    /* Set the WIMP parameters.
       Specify the WIMP-nucleon couplings, with
       the WIMP mass [GeV] always the first argument:
         * DDCalc::SetWIMP_mfa(m,fp,fn,ap,an)
           The standard couplings fp,fn [GeV^-2] & ap,an [unitless]
       In the above, 'p' is for proton, 'n' is for neutron. */
    switch (type)
    {
      case TYPE_MFA:
        DDCalc::SetWIMP_mfa(WIMP,M,xpSI,xnSI,xpSD,xnSD);
        break;
    }
    
    /* Get the WIMP parameters with the same signatures and units as
       above.  The only difference is that WIMP-nucleon cross-sections
       are always positive. */
    DDCalc::GetWIMP_mfa(WIMP,M,fp,fn,ap,an);
    
    /* Print out the above WIMP mass, couplings, and cross sections. */
    printf("%s %- #12.5g\n","WIMP mass [GeV]     ",M);
    std::cout << std::endl;
    printf("%-28s %11s %11s %11s %11s\n","WIMP-nucleon couplings",
           " proton-SI "," neutron-SI"," proton-SD "," neutron-SD");
    printf("%-28s %- #11.5g %- #11.5g %- #11.5g %- #11.5g\n",
           "  f & a [GeV^-2,unitless]",fp,fn,ap,an);
    std::cout << std::endl;
    
    /* Do rate calculations using the specified WIMP and halo parameters.
       Does the calculations necessary for predicted signals, likelihoods
       and/or maximum gap statistics. */
    DDCalc::CalcRates(XENON,WIMP,Halo);
    DDCalc::CalcRates(LUX,WIMP,Halo);
    DDCalc::CalcRates(SCDMS,WIMP,Halo);
    DDCalc::CalcRates(SIMPLE,WIMP,Halo);
    DDCalc::CalcRates(MyDetector,WIMP,Halo);
    
    /* Header */
    printf("%-20s  %11s  %11s  %11s  %11s  %11s\n","",
           " XENON 2012"," LUX 2013  ","SuCDMS 2014","SIMPLE 2014",
           " (special) ");
    //printf("%-20s  %11s  %11s  %11s  %11s  %11s\n","",
    //       "-----------","-----------","-----------","-----------",
    //       "-----------","-----------");
    
    /* Event quantities. */
    printf("%-20s  % 6i       % 6i       % 6i       % 6i       % 6i     \n",
           "Observed events     ",
        DDCalc::Events(XENON),
        DDCalc::Events(LUX),
        DDCalc::Events(SCDMS),
        DDCalc::Events(SIMPLE),
        DDCalc::Events(MyDetector));
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Expected background ",
        DDCalc::Background(XENON),
        DDCalc::Background(LUX),
        DDCalc::Background(SCDMS),
        DDCalc::Background(SIMPLE),
        DDCalc::Background(MyDetector));
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Expected signal     ",
        DDCalc::Signal(XENON),
        DDCalc::Signal(LUX),
        DDCalc::Signal(SCDMS),
        DDCalc::Signal(SIMPLE),
        DDCalc::Signal(MyDetector));

    
    /* The log-likelihoods for the current WIMP; note these are _not_
       multiplied by -2.  The likelihood is calculated using a Poisson
       given the observed x.length() of events and expected signal +
       background. */
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Log-likelihood      ",
        DDCalc::LogLikelihood(XENON),
        DDCalc::LogLikelihood(LUX),
        DDCalc::LogLikelihood(SCDMS),
        DDCalc::LogLikelihood(SIMPLE),
        DDCalc::LogLikelihood(MyDetector));
    
    /* The logarithm of the p-value, calculated without background
       subtraction, using either the maximum gap statistic or a Poisson
       statistic, depending on how the detector was initialized.  Note
       that this is actually a conservative upper _bound_ on the p-value
       in the event of an unknown background and is useful for excluding
       WIMP parameters.  However, since it is not a true p-value, it
       should not be interpreted as being related to any particular
       likelihood. */
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Max gap log(p-value)",
        DDCalc::LogPValue(XENON),
        DDCalc::LogPValue(LUX),
        DDCalc::LogPValue(SCDMS),
        DDCalc::LogPValue(SIMPLE),
        DDCalc::LogPValue(MyDetector));
    
    /* Returns a factor x by which the current WIMP cross-sections must
       be multiplied (sigma -> x*sigma, applied to all four WIMP-nucleon
       cross-sections) to achieve the given p-value (specified by its
       logarithm).  Useful for finding the no-background-subtraction
       exclusion limits.  For example, if setWIMP_msigma(100.0,10.0,
       10.0,0.0,0.0) is called, then x*(10. pb) would be the SI
       cross-section at a WIMP mass of 100 GeV at which the experiment
       is excluded at the 90% CL (p=1-CL). */
    double lnp = log(0.1);  // default value for optional argument
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Max gap x for 90% CL",
        DDCalc::ScaleToPValue(XENON),
        DDCalc::ScaleToPValue(LUX),
        DDCalc::ScaleToPValue(SCDMS),
        DDCalc::ScaleToPValue(SIMPLE),
        DDCalc::ScaleToPValue(MyDetector));
    std::cout << "         * Factor x such that sigma->x*sigma gives desired p-value" << std::endl;

  }  // END INPUT LOOP <<<<<<<<<<<<<<<<<<<<<<<<<

  // Clean up all the objects
  DDCalc::FreeAll();

} 


// UTILITY FUNCTION DEFINITIONS ----------------------------------------

/* Write a description of how input parameters should be specified. */
void WriteDescription(const int type)
{
  std::cout << std::endl;
  std::cout << "Enter WIMP parameters below.  Only the first two are necessary." << std::endl;
  std::cout << "A blank line terminates input.  The parameters are:" << std::endl;
  //std::cout << "Enter WIMP parameters below.  The parameters are:" << std::endl;
  std::cout << std::endl;
  switch (type)
  {
    case TYPE_MFA:
      std::cout << "  M     WIMP mass [GeV]" << std::endl;
      std::cout << "  fp    Spin-independent WIMP-proton effective coupling [GeV^-2]" << std::endl;
      std::cout << "  fn    Spin-independent WIMP-neutron effective coupling [GeV^-2]" << std::endl;
      std::cout << "  ap    Spin-dependent WIMP-proton effective coupling [unitless]" << std::endl;
      std::cout << "  an    Spin-dependent WIMP-neutron effective coupling [unitless]" << std::endl;
      break;
  }
}


/* Read WIMP parameters (mass & couplings) from standard input. */
bool GetWIMPParams(const int type, double& M, double& xpSI, double& xnSI,
                   double& xpSD, double& xnSD)
{
  std::cout << std::endl;
  std::cout << "------------------------------------------------------------" << std::endl;
  switch (type)
  {
    case TYPE_MFA:
      std::cout << "Enter values <M fp fn ap an>:" << std::endl;
      break;
    }
  
  // Get input line for parsing
  std::string line;
  getline(std::cin, line);
  std::istringstream iss(line);

  // Parse input line
  if (!(iss >> M)) return false;
  if (!(iss >> xpSI)) return false;
  if (!(iss >> xnSI))
  {
    xnSI = xpSI; xpSD = 0.0; xnSD = 0.0;
    return true;
  }
  if (!(iss >> xpSD))
  {
    xpSD = 0.0; xnSD = 0.0;
    return true;
  }
  if (!(iss >> xnSD))
  {
    xnSD = xpSD;
    return true;
  }
  return true;
}


// END FILE ------------------------------------------------------------

 
