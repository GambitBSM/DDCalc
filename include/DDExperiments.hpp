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
  int C_DDCalc_xenon100_2012_init(const bool&);
  int C_DDCalc_lux_2013_init(const bool&);
  int C_DDCalc_lux_2016_init(const bool&);
  int C_DDCalc_pandax_2016_init(const bool&);
  int C_DDCalc_xenon1t_2017_init(const bool&);
  int C_DDCalc_lux_2015_init(const bool&);
  int C_DDCalc_pico_2l_init(const bool&);
  int C_DDCalc_pico_60_f_init(const bool&);
  int C_DDCalc_pico_60_i_init(const bool&);
  int C_DDCalc_pico_60_2017_init(const bool&);
  int C_DDCalc_supercdms_2014_init(const bool&);
  int C_DDCalc_simple_2014_init(const bool&);
  int C_DDCalc_darwin_ar_init(const bool&);
  int C_DDCalc_darwin_xe_init(const bool&);

}
  

// EXTERNAL (USER) ROUTINES --------------------------------------------

namespace DDCalc
{

  // Initialize experiments for which likelihoods are to be calculated.
  // The argument indicates if extra calculations necessary for maximum gap
  // statistics should be performed (unnecessary for likelihoods).  The
  // flag is ignored for experimental analysis lacking the event energies
  // necessary for a maximum gap analysis.

  int XENON100_2012_Init(const bool intervals=true)
  {
    C_DDCalc_xenon100_2012_init(intervals);
  }
  
  int LUX_2013_Init(const bool intervals=true)
  {
    C_DDCalc_lux_2013_init(intervals);
  }
  
  int LUX_2016_Init(const bool intervals=true)
  {
    C_DDCalc_lux_2016_init(intervals);
  }
  
  int PandaX_2016_Init(const bool intervals=true)
  {
    C_DDCalc_pandax_2016_init(intervals);
  }

  int Xenon1T_2017_Init(const bool intervals=true)
  {
    C_DDCalc_xenon1t_2017_init(intervals);
  }

  int LUX_2015_Init(const bool intervals=true)
  {
    C_DDCalc_lux_2015_init(intervals);
  }

  int PICO_2L_Init(const bool intervals=true)
  {
    C_DDCalc_pico_2l_init(intervals);
  }

  int PICO_60_F_Init(const bool intervals=true)
  {
    C_DDCalc_pico_60_f_init(intervals);
  }

  int PICO_60_I_Init(const bool intervals=true)
  {
    C_DDCalc_pico_60_i_init(intervals);
  }

  int PICO_60_Init(const bool intervals=true)
  {
    C_DDCalc_pico_60_init(intervals);
  }

  int PICO_60_2017_Init(const bool intervals=true)
  {
    C_DDCalc_pico_60_2017_init(intervals);
  }
  
  int SuperCDMS_2014_Init(const bool intervals=true)
  {
    C_DDCalc_supercdms_2014_init(intervals);
  }
  
  int SIMPLE_2014_Init(const bool intervals=true)
  {
    C_DDCalc_simple_2014_init(intervals);
  }
  
  int DARWIN_Ar_Init(const bool intervals=true)
  {
    C_DDCalc_darwin_ar_init(intervals);
  }
  
  int DARWIN_Xe_Init(const bool intervals=true)
  {
    C_DDCalc_darwin_xe_init(intervals);
  }
    
}

#endif  // DDCalc_H

