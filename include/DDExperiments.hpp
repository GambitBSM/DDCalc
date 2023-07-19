/**********************************************************************
 *
 * C/C++ interface to experiment-specific DDCalc routines.
 * 
 * Expand these if you want to add a new experiment.
 * You will also need to add the new experiment in
 * src/DDExperiments.f90 and src/analyses/<your_experiment_name>.f90
 *
 *      P Scott        Imperial College 2016
 *
 **********************************************************************/

#ifndef DDExperiments_H
#define DDExperiments_H

// INTERNAL ROUTINES ---------------------------------------------------

extern "C"
{
  
  // Experiment initialization (factory functions)
  int C_DDCalc_xenon100_2012_init();
  int C_DDCalc_lux_2013_init();
  int C_DDCalc_lux_2016_init();
  int C_DDCalc_pandax_2016_init();
  int C_DDCalc_pandax_2017_init();
  int C_DDCalc_xenon1t_2017_init();
  int C_DDCalc_xenon1t_2018_init();
  int C_DDCalc_xenon1t_2022_init();
  int C_DDCalc_lux_2015_init();
  int C_DDCalc_pico_2l_init();
  int C_DDCalc_pico_60_init();
  int C_DDCalc_pico_60_2017_init();
  int C_DDCalc_pico_60_2019_init();
  int C_DDCalc_supercdms_2014_init();
  int C_DDCalc_cdmslite_init();
  int C_DDCalc_simple_2014_init();
  int C_DDCalc_cresst_ii_init();
  int C_DDCalc_dummyexp_init();
  int C_DDCalc_lz_init();
  int C_DDCalc_lz_2022_init();
  int C_DDCalc_pandax_4t_init();
  int C_DDCalc_darwin_init();
  int C_DDCalc_darkside_20k_init();
  int C_DDCalc_darkside_50_init();
  int C_DDCalc_darkside_50_s2_init();
  int C_DDCalc_pico_500_init();

}
  

// EXTERNAL (USER) ROUTINES --------------------------------------------

namespace DDCalc
{

  // Initialize experiments for which likelihoods are to be calculated.
  // The argument indicates if extra calculations necessary for maximum gap
  // statistics should be performed (unnecessary for likelihoods).  The
  // flag is ignored for experimental analysis lacking the event energies
  // necessary for a maximum gap analysis.

  int XENON100_2012_Init()
  {
    C_DDCalc_xenon100_2012_init();
  }
  
  int LUX_2013_Init()
  {
    C_DDCalc_lux_2013_init();
  }
  
  int LUX_2016_Init()
  {
    C_DDCalc_lux_2016_init();
  }
  
  int PandaX_2016_Init()
  {
    C_DDCalc_pandax_2016_init();
  }

  int PandaX_2017_Init()
  {
    C_DDCalc_pandax_2017_init();
  }

  int XENON1T_2017_Init()
  {
    C_DDCalc_xenon1t_2017_init();
  }

  int XENON1T_2018_Init()
  {
    C_DDCalc_xenon1t_2018_init();
  }

  int XENON1T_2022_Init()
  {
    C_DDCalc_xenon1t_2022_init();
  }

  int LUX_2015_Init()
  {
    C_DDCalc_lux_2015_init();
  }

  int PICO_2L_Init()
  {
    C_DDCalc_pico_2l_init();
  }

  int PICO_60_Init()
  {
    C_DDCalc_pico_60_init();
  }

  int PICO_60_2017_Init()
  {
    C_DDCalc_pico_60_2017_init();
  }
 
  int PICO_60_2019_Init()
  {
    C_DDCalc_pico_60_2019_init();
  }
 
  int SuperCDMS_2014_Init()
  {
    C_DDCalc_supercdms_2014_init();
  }
  
  int SIMPLE_2014_Init()
  {
    C_DDCalc_simple_2014_init();
  }

  int CRESST_II_Init()
  {
    C_DDCalc_cresst_ii_init();
  }

  int LZ_2022_Init()
  {
    C_DDCalc_lz_2022_init();
  }

  int PandaX_4T_Init()
  {
    C_DDCalc_pandax_4t_init();
  }

  int DARWIN_Init()
  {
    C_DDCalc_darwin_init();
  }

  int DarkSide_20k_Init()
  {
    C_DDCalc_darkside_20k_init();
  }

  int DarkSide_50_Init()
  {
    C_DDCalc_darkside_50_init();
  }

  int DarkSide_50_S2_Init()
  {
    C_DDCalc_darkside_50_s2_init();
  }

  int PICO_500_Init()
  {
    C_DDCalc_pico_500_init();
  }
    
}

#endif  // DDCalc_H

