#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include "DDCalc.hpp"

int main()
{  
  
  double BGlogL,mDM,sigmatest,limit,logLlimit,mDMmin,mDMmax, sigmin, sigmax, sig, logLbest, mDMbest, sigbest, s, s1, s2;
  int mDMsteps, mDMi, sigsteps, sigi;

  int WIMP;
  int Halo;
  int Detector;  

// Define the smallest and largest DM mass to be considered, as well as the number of intermediate values
  mDMmin = 1.0;
  mDMmax = 1.0e4;
  mDMsteps = 80;

  sigmin = 1.0e-11;
  sigmax = 1.0e-7;
  sigsteps = 80;

  double logL[mDMsteps+1][sigsteps+1];

// Define the (arbitrary) reference cross section (in pb) for which event rates will be calculated
  sigmatest = 1.0e-10;

// Initialize WIMP and DM Halo object with default values
  Halo = DDCalc::InitHalo();
  WIMP = DDCalc::InitWIMP();	

// Optional: Change the parameters of the Standard Halo Model to correspond to standard conventions
  DDCalc::SetSHM(Halo, 0.3, 220.0, 220.0);

  Detector = DDCalc::XENON1T_2018_Init();

// *************************************************
// *** Routines for calculating exclusion limits ***
// *************************************************

// Here we assume that the experiment sees no excess, such that the best-fit point is identical to the background-only hypothesis.
// See below for the case of experiments seeing an excess.

// Step 1: Calculate the log likelihood for the background-only hypothesis
// This can be achieved by setting the DM couplings equal to zero (in which case the DM mass is irrelevant)
  DDCalc::SetWIMP_msigma(WIMP, mDMmin, 0., 0., 0., 0.); 
  DDCalc::CalcRates(Detector,WIMP,Halo);
  BGlogL = DDCalc::LogLikelihood(Detector);

// Step 2: Calculate the critical value of the log likelihood that corresponds to the exclusion limit
// Here we assume the asymptotic limit, such that -2 delta log L follows a 1/2 chi^2 distribution with 1 degree of freedom
// The one-sided upper bound at 90% confidence level is then given by -2 delta log L = 1.64
  limit = 1.64;
  logLlimit = BGlogL - (limit / 2.0);

  printf("%f",logLlimit);

// Step 3: Calculate the spin-independent exclusion limit (assuming equal couplings to protons and neutrons)
  printf("\n");
  printf("**************************** Spin-independent limit ****************************\n");
  printf("\n");
  printf("Assuming spin-independent scattering with equal couplings to protons and neutrons.\n");
  printf("\n");
  printf("    log_10(m_DM/GeV)  log_10(sigma_p/pb)\n");
  printf("\n");

  for( mDMi = 0; mDMi <= mDMsteps; mDMi = mDMi + 1 )
  {
    mDM = mDMmin * pow(mDMmax/mDMmin, mDMi/( (double) mDMsteps) );
    DDCalc::SetWIMP_msigma(WIMP, mDM, sigmatest, sigmatest, 0., 0.);
    DDCalc::CalcRates(Detector,WIMP,Halo);
    printf("          %10.4f          %10.4f\n",log(mDM)/log(10.0),log(DDCalc::ScaleToPValue(Detector,logLlimit)*sigmatest)/log(10.0));
  }


// Step 4: Calculate spin-dependent exclusion limit (assuming couplings only to protons)
  printf("\n");
  printf("***************************** Spin-dependent limit *****************************\n");
  printf("\n");
  printf("Assuming spin-dependent scattering with couplings to protons only.\n");
  printf("\n");
  printf("    log_10(m_DM/GeV)  log_10(sigma_p/pb)\n");
  printf("\n");

  for( mDMi = 0; mDMi <= mDMsteps; mDMi = mDMi + 1 )
  {
    mDM = mDMmin * pow(mDMmax/mDMmin, mDMi/( (double) mDMsteps) );
    DDCalc::SetWIMP_msigma(WIMP, mDM, 0., 0., sigmatest, 0.);
    DDCalc::CalcRates(Detector,WIMP,Halo);
    printf("          %10.4f          %10.4f\n",log(mDM)/log(10.0),log(DDCalc::ScaleToPValue(Detector,logLlimit)*sigmatest)/log(10.0));
  }


/*************************************************
 *** Routines for calculating best-fit regions ***
 *************************************************/

// Now we consider a more general case of an experiment seeing an excess,
// i.e. we are interested in constructing confidence intervals / regions.
// See arXiv:1407.6617 for details.

// Step 5: Feldman-Cousins
// If we consider only the total number of expected and observed events, we can employ the Feldman-Cousins method
// to construct a lower and an upper bound on the spin-independent scattering cross section as a function of DM mass

  s2=DDCalc::FeldmanCousinsUpper(-2.3,DDCalc::Events(Detector),DDCalc::Background(Detector));
  s1=DDCalc::FeldmanCousinsLower(-2.3,DDCalc::Events(Detector),DDCalc::Background(Detector));

  printf("\n");
  printf("**************************** Feldmann-Cousins bound ****************************\n");
  printf("\n");
  printf("Assuming spin-independent scattering with equal couplings to protons and neutrons.\n");
  printf("\n");
  if(s1 > 0)  
  {
    printf("      log_10(m_DM/GeV)  log_10(sigma_min/pb)  log_10(sigma_max/pb)\n");
    printf("\n");
  }
  else 
  {
    printf("Lower bound corresponds to zero cross section. Quoting upper bound only.\n");
    printf("\n");
    printf("      log_10(m_DM/GeV)  log_10(sigma_max/pb)\n");
    printf("\n");
  }

  for( mDMi = 0; mDMi <= mDMsteps; mDMi = mDMi + 1 )
  {
    mDM = mDMmin * pow(mDMmax/mDMmin, mDMi/( (double) mDMsteps) );
    DDCalc::SetWIMP_msigma(WIMP, mDM, sigmatest, sigmatest, 0., 0.);
    DDCalc::CalcRates(Detector,WIMP,Halo);
    s = DDCalc::Signal(Detector);
    if(s > 0)
    {
      if(s1 > 0)
      {
        printf("            %10.4f            %10.4f            %10.4f\n",log(mDM)/log(10.),log(s1/s*sigmatest)/log(10.),log(s2/s*sigmatest)/log(10.));
      }
      else
      {
        printf("            %10.4f            %10.4f\n",log(mDM)/log(10.),log(s2/s*sigmatest)/log(10.));
      }
    }
  }

// Step 6: 2d likelihood scan
// To include spectral information, we need to perform a full scan over the 2d parameter space
// We can then determine the combination of mass and cross section that give the best fit to the data and construct
// confidence regions around this point.
  logLbest = BGlogL;
  for( mDMi = 0; mDMi <= mDMsteps; mDMi = mDMi + 1 )
  {
    mDM = mDMmin * pow(mDMmax/mDMmin, mDMi/( (double) mDMsteps) );
    for( sigi = 0; sigi <= sigsteps; sigi = sigi + 1 )
    {
      sig = sigmin * pow(sigmax/sigmin, sigi/( (double) sigsteps) );
      DDCalc::SetWIMP_msigma(WIMP, mDM, sig, sig, 0., 0.);
      DDCalc::CalcRates(Detector,WIMP,Halo);
      logL[mDMi][sigi] = DDCalc::LogLikelihood(Detector);
      if(logL[mDMi][sigi] > logLbest)
      {
        logLbest = logL[mDMi][sigi];
        mDMbest = mDM;
        sigbest = sig;
      }
    }
  }

  printf("\n");
  printf("****************************** 2d likelihood scan ******************************\n");
  printf("\n");
  printf("Assuming spin-independent scattering with equal couplings to protons and neutrons.\n");
  printf("\n");
  printf(" Best fit point is m_DM = %7.1f GeV, sigma_p = %10.3e pb with log L = %7.3f\n",mDM,sigbest,logLbest);
  printf("\n");
  printf("    log_10(m_DM/GeV)  log_10(sigma_p/pb)               log L      -2 Delta log L\n");
  printf("\n");
  for( mDMi = 0; mDMi <= mDMsteps; mDMi = mDMi + 1 )
  {
    mDM = mDMmin * pow(mDMmax/mDMmin, mDMi/( (double) mDMsteps) );
    for( sigi = 0; sigi <= sigsteps; sigi = sigi + 1 )
    {
      sig = sigmin * pow(sigmax/sigmin, sigi/( (double) sigsteps) );
      printf("          %10.4f          %10.4f          %10.4f          %10.4f\n", log(mDM)/log(10.),log(sig)/log(10.),logL[mDMi][sigi],-2.0*(logL[mDMi][sigi]-logLbest));
    }
  }

}

