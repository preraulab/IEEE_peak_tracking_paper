This Peak Tracking Toolbox contains files for tracking parameterized peaks 
in the spectrogram. These include a set of variant filters, several model classes
used by the filters, functions for plotting the estimates, files to compare the 
variant filter performances on simulated examples, and files for analyzing a 
real data example. 

The primary user-level files are:
  ex_real_data.m -- script to run the real data example. Can be run as is.
  compareEKFs_jobs_wRandSSMs.m -- function to estimate/compare filters on simulated examples. Can be run with defaults as is. 
  iekf.m -- the recommended iekf filter function.
  plotMPSpectWithSlices.m -- function to plot observed spectrogram with estimate.
  plotMPStateEstimates.m -- function to plot state estimates.

A complete listing is given below.

NOTE: The StateSpaceMultiPeak class was built with the intention of 
      allowing general SDEs for the state evolution equation. However, this 
      was never fully realized. The in-class simulate function and the 
      external filter functions impose very specific forms of the state 
      evolution and noise covariances. The function can be modified 
      accordingly to handle more general situations. 

*******************
* Variant filters * 
*******************
iekf -- function for main IEKF with discrete switching of peak-combination
  and improved initial trajectory using random draws. A basic IEKF version 
  without draws is obtained by setting num_particles = 1. An EKF version (with 
  or without draws) is obtained by setting num_iters = 1.  

iekfW2ndDeriv -- similar to iekf, except the filter updates are treated 
  as a Newton-Raphson mode optimization that utilizes the second derivatives.

ekf2ndOrder -- function for a second-order Gaussian EKF.

iekfWPostMode -- same as iekf, except it selects the draw of the initial 
  reference trajectory as the approximate posterior filter mode.

NOTE: For computation speed, each filter currently overwrites the state 
      evolution equation, assuming x(t) = 0.9*Phi*x(t-1) + w(t), 
      with Phi(t1,t2) = I. They also assume constant noise covariances Q and R.
      The filters must be modified accordingly to handle other situations.

*****************
* Model Classes * 
*****************
StateSpaceMultiPeak -- state-space model class extending the MultiPeakModel class.
  It adds the variables needed for forming the state-space model. It also
  has a simulate function, which will simulate the model. 
  NOTE: The simulate function currently assumes an identity transition 
        matrix for the state and constant noise covariances.

MultiPeakModel -- class for a multi-peak observation. It has the individual 
  component peaks, functions to access their properties, functions to evaluate 
  the multi-peak observation and its derivatives at given parameter values, 
  and functions for plotting.
   
PeakModel -- class for individual peaks. It contains the peak type, 
  peak parameters (with their bounds), the peak function and its derivatives.

PeakObjParam -- class for individual parameters. It defines their type of 
  bound and contains the bound function and its derivatives and inverse. 

makePeakCombos -- function to determine maximum available peak combinations
  based on which are dynamic. Could be moved inside StateSpaceMultiPeak class.

makeComboTransitionMatr -- function to form the transition matrix of the 
  On/Off-peaks combination. Could be moved inside StateSpaceMultiPeak class.

****************************************************************
* Files to compare filter variatants applied to simulated data *
****************************************************************
compareEKFs_jobs_wRandSSMs -- Function to generate simulated examples, compute 
  the variant filter estimates for each, and compute a set of statistical 
  performance measures. The first four inputs--an output directory (odir), 
  simulation type (sim_type), number of peaks (num_peaks), and number of 
  simulations (num_sims) are input as strings to enable batch deployment on 
  a cluster. 
  sim_type determines whether the peak parameters follow hard-bound random 
  walks ('randWalk') or pseudo-deterministic patterns ('pseudoDeterm'). 
  num_peaks is the number of peaks on an exponential-decay background.
  Can be run for default (one peak on decay background 
  with hard-bound random walk parameters) as is.  
  NOTE: May hit memory issue at call to moranI in checkEstimates, in which 
        case comment out either function as desired.

simulatePseudoDeterm -- function to generate simulation with peak parameters 
  following pseudo-deterministic parameters. 

checkEstimates -- function that computes a variety of statistical measures 
  from the various filter estimates for each simulation.

boxQ -- function to comptute the Box Q statistic. It is modified to be able 
  to limit the extent of the two-dimensional correlations evaluated. The default
  is 10 time lags, all frequency bins. Called in checkEstimates using 
  4 time bins and +/-4 frequency bins.

moranI -- function to compute Moran's I statistic over. It is modified to be able 
  to limit the extent of the two-dimensional correlations evaluated. The default
  is all time lags, all frequency bins. Called in checkEstimates using 
  +/-3 time bins and +/-4 frequency bins.

myDateStr -- utility function used to generate data string for filename in compareEKFs_jobs_wRandSSMs.m

*******************************
* Files for real data example *
*******************************
sop_spect_mt.mat -- example real data file containing pre-computed spectrogram (spect)
  with time and frequency bins (stimes and sfreqs, respectively) and sampling rate (Fs). 

ex_real_data -- script to analyze data in sop_spec_mt.mat. Can be run as is.

exp_baseline_fit -- function to determine initial state values of exponential-decay 
  background function. Called in ex_real_data.

nanpow2db -- utility function to convert to spectrogram to dB when data contains NaNs.

*************************************************
* Files for plotting filter and state estimates *
*************************************************
plotMPSpectWithSlices -- function to plot original spectrogram, estimated spectrogram,
  residual estimates, and estimated On/Off-peak combinations. There is also 
  option to plot interactive time-slices of the spectrum, estimate, and individual peaks.

plotMPStateEstimates -- function to plot state estimates.

climscale -- utility function to adjust color scale in imagesc spectrogram plots.

