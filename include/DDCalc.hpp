/**********************************************************************
 *
 * C/C++ interface to various DDCalc routines.
 *
 *      C Savage       Nordita          2014
 *      A Scaffidi     U of Adelaide    2014
 *      M White        U of Adelaide    2014
 *      P Scott        Imperial College 2016
 *
 **********************************************************************/

#ifndef DDCalc_H
#define DDCalc_H

#include "DDExperiments.hpp"

// INTERNAL ROUTINES ---------------------------------------------------
// Here we just declare the symbols presented by the DDCalc library.
// Ignore this unless you want to import symbols yourself directly;
// otherwise just include this header and call the external routines
// defined further down.

extern "C"
{
    
  // Default initialisation (factory functions)
  int C_DDWIMP_ddcalc_initwimp();
  int C_DDHalo_ddcalc_inithalo();
  int C_DDExperiments_ddcalc_initdetector();

  // Halo setter
  void C_DDCalc_ddcalc_setshm(const int&, const double&, const double&, const double&,
                              const double&);
  
  // WIMP parameter setters and getters
  void C_DDCalc_ddcalc_setwimp_mfa(const int&, const double&,
                  const double&, const double&, const double&, const double&);
  void C_DDCalc_ddcalc_getwimp_mfa(const int&, double&,
                  double&, double&, double&, double&);
  void C_DDCalc_ddcalc_setwimp_mg(const int&, const double&,
                  const double&, const double&, const double&, const double&);
  void C_DDCalc_ddcalc_getwimp_mg(const int&, double&,
                  double&, double&, double&, double&);
  void C_DDCalc_ddcalc_setwimp_higgsportal(const int&, const double&,
                  const double&, const double&, const double&, const double&);
  void C_DDCalc_ddcalc_getwimp_higgsportal(const int&, double&,
                  double&, double&, double&, double&);
  
  // Detector parameter setter (minimum recoil energy to consider)
  void C_DDCalc_ddcalc_setdetectoremin(const int&, const double&);

  // Perform calculations
  void C_DDRates_ddcalc_calcrates(const int&, const int&, const int&);
  
  // Inspect results of calculations
  int C_DDRates_ddcalc_events(const int&);          // Number of events
  double C_DDRates_ddcalc_background(const int&);   // Expected backgrounds
  double C_DDRates_ddcalc_signal(const int&);       // Expected signal
  double C_DDStats_ddcalc_loglikelihood(const int&); // Log-likelihood
  double C_DDStats_ddcalc_scaletopvalue(const int&, const double&);
   // Factor x by which sigma -> x*sigma would yield given p-value (given as log(p))

  // Do memory cleanup
  void C_DDUtils_ddcalc_freewimps();
  void C_DDUtils_ddcalc_freehalos();
  void C_DDUtils_ddcalc_freedetectors();
  void C_DDUtils_ddcalc_freeall();

}


// EXTERNAL (USER) ROUTINES --------------------------------------------

namespace DDCalc
{

  //########## Default initialisation (factory functions) ##############
  
  // Create a WIMPStruct object internally in DDCalc with default values
  // and return an index for it. 
  int InitWIMP() { return C_DDWIMP_ddcalc_initwimp(); }

  // Create a HaloStruct object internally in DDCalc with default values
  // and return an index for it. 
  int InitHalo() { return C_DDHalo_ddcalc_inithalo(); }

  // Create a DetectorStruct object internally in DDCalc with default values
  // and return an index for it. 
  int InitDetector() { return C_DDExperiments_ddcalc_initdetector(); }


  //########## Halo setter #############################################
    
  // (Re-)Set the Standard Halo Model parameters:
  //    rho     Local dark matter density [GeV/cm^3]
  //    vrot    Local disk rotation speed [km/s]
  //    v0      Maxwell-Boltzmann most probable speed [km/s]
  //    vesc    Galactic escape speed [km/s]
  // This example uses the default values (and is thus optional).
  void SetSHM(const int HaloIndex, const double rho=0.4, const double vrot=235.0,
              const double v0=235.0, const double vesc=550.0)
  {
    C_DDCalc_ddcalc_setshm(HaloIndex,rho,vrot,v0,vesc);
  }

  //########## WIMP parameter setters and getters ######################
  
  // (Re-)Set the WIMP parameters.
  // There are three ways to specify the WIMP-nucleon couplings, with the
  // WIMP mass [GeV] always the first argument:
  //   * SetWIMP_mfa(m,fp,fn,ap,an)
  //     The standard couplings fp,fn [GeV^-2] & ap,an [unitless]
  //   * SetWIMP_mG(m,GpSI,GnSI,GpSD,GnSD)
  //     The effective 4 fermion vertex couplings GpSI,GnSI,GpSD,GnSD
  //     [GeV^-2], related by:
  //         GpSI = 2 fp        GpSD = 2\sqrt{2} G_F ap
  //         GnSI = 2 fn        GnSD = 2\sqrt{2} G_F an
  // In the above, 'p' is for proton, 'n' is for neutron, 'SI' is for
  // spin-independent, and 'SD' is for spin-dependent.
  void SetWIMP_mfa(const int WIMPIndex, const double m, const double fp, const double fn,
                   const double ap, const double an)
  {
    C_DDCalc_ddcalc_setwimp_mfa(WIMPIndex,m,fp,fn,ap,an);
  }

  void SetWIMP_mG(const int WIMPIndex, const double m, const double GpSI, const double GnSI,
                   const double GpSD, const double GnSD)
  {
    C_DDCalc_ddcalc_setwimp_mg(WIMPIndex,m,GpSI,GnSI,GpSD,GnSD);
  }

  void SetWIMP_Higgsportal(const int WIMPIndex, const double m, const double fsp, const double fsn,
                   const double app, const double apn)
  {
    C_DDCalc_ddcalc_setwimp_higgsportal(WIMPIndex,m,fsp,fsn,app,apn);
  }
    
  // Get the WIMP parameters with the same signatures and units as above.
  // The only difference is that WIMP-nucleon cross-sections are always
  // positive.
  void GetWIMP_mfa(const int WIMPIndex, double& m, double& fp, double& fn,
                          double& ap, double& an)
  {
    C_DDCalc_ddcalc_getwimp_mfa(WIMPIndex,m,fp,fn,ap,an);
  }

  void GetWIMP_mG(const int WIMPIndex, double& m, double& GpSI, double& GnSI,
                          double& GpSD, double& GnSD)
  {
    C_DDCalc_ddcalc_getwimp_mg(WIMPIndex,m,GpSI,GnSI,GpSD,GnSD);
  }

  void GetWIMP_Higgsportal(const int WIMPIndex, double& m, double& fsp, double& fsn,
                          double& app, double& apn)
  {
    C_DDCalc_ddcalc_getwimp_higgsportal(WIMPIndex,m,fsp,fsn,app,apn);
  }
  
  //########## Detector setter #########################################

  // Minimum recoil energy to consider (keV)
  void SetDetectorEmin(const int DetectorIndex, const double Emin)
  {
    C_DDCalc_ddcalc_setdetectoremin(DetectorIndex,Emin);
  }

  //########## Rate calculations #######################################
   
  // Perform the rate calculations necessary for the likelihoods.
  void CalcRates(const int DetectorIndex, const int WIMPIndex, const int HaloIndex)
  {
    C_DDRates_ddcalc_calcrates(DetectorIndex, WIMPIndex, HaloIndex);
  }
  
  //########## Inspect results of calculations #########################

  // Read off the observed events for the specified analysis.
  int Events(const int DetectorIndex)
  {
    return C_DDRates_ddcalc_events(DetectorIndex);
  }
  
  // Read off the expected background for the specified analysis.
  double Background(const int DetectorIndex)
  {
    return C_DDRates_ddcalc_background(DetectorIndex);
  }
  
  // Read off the expected signal for the specified analysis.
  double Signal(const int DetectorIndex)
  {
    return C_DDRates_ddcalc_signal(DetectorIndex);
  }
 
  
  // Read off the log-likelihoods for the specificed analysis; note these
  // are _not_ multiplied by -2.  The likelihood is calculated using a Poisson
  // given the observed number of events and expected signal + background.
  double LogLikelihood(const int DetectorIndex)
  {
    return C_DDStats_ddcalc_loglikelihood(DetectorIndex);
  }
  
  // Returns a factor x by which the WIMP cross-sections must
  // be multiplied (sigma -> x*sigma, applied to all four WIMP-nucleon
  // cross-sections) to achieve the given p-value (specified by its
  // logarithm).  Useful for finding the no-background-subtraction
  // exclusion limits.  For example, if setWIMP_msigma(100.0,10.0,10.0,
  // 0.0,0.0) is called, then x*(10. pb) would be the SI cross-section
  // at a WIMP mass of 100 GeV at which the experiment is excluded at
  // the 90% CL (p=1-CL).
  double ScaleToPValue(const int DetectorIndex, const double logp=-2.302585)
  {
    return C_DDStats_ddcalc_scaletopvalue(DetectorIndex, logp);
  }

  //########## Do memory cleanup #######################################

  // Delete all WIMPStructs
  void FreeWIMPs() { C_DDUtils_ddcalc_freewimps(); }

  // Delete all HaloStructs
  void FreeHalos() { C_DDUtils_ddcalc_freehalos(); }

  // Delete all DetectorStructs
  void FreeDetectors() { C_DDUtils_ddcalc_freedetectors(); }

  // Delete all WIMPStructs, HaloStructs and DetectorStructs
  void FreeAll() { C_DDUtils_ddcalc_freeall(); }
    
}

#endif  // DDCalc_H

