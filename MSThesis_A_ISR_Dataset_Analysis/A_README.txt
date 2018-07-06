MAIN A: Code to Read, Prepare, and Analyze the ISR Dataset: Ross_Test.m
This code reads, prepares, and analyzes the ISR dataset in a multitude of ways. It was developed on MATLAB R2016b.

Dataset can be found at: http://madrigal.haystack.mit.edu/cgi-bin/madrigal/madExperiment.cgi?exp=experiments/2013/mlh/06may13&displayLevel=1&expTitle=Rapid%20TID as mlh130506g.001 which should be downloaded as an HDF5 file. The newer file mlh130506g.002 can also be used by extracting the alternating code (97) data and modifying the import code a bit. The newer file was not back-ported for the plots used in this paper because the SNR and ion velocity datasets did not visibly change.

This research was supported under NSF Grant AGS 12-41407 to The Pennsylvania State University.

Supporting functions are also provided to allow the code to be run.

FUN A-1: Function to Convert Date to Julian Date: date2jd.m
This function converts the date to the Julian date for aligning the time correctly. It was coded by Peter J. Acklam. It was run on MATLAB R2016b.

FUN A-2: Function to Run Lomb-Scargle Spectral Analysis: Ross_Lombscargle_optimized.m
This function quickly runs a Lomb-Scargle spectral analysis on a given time series that does not need to be evenly spaced in time. It was written by Dr. Brett Shoelson, modified by Dr. Sumanta Sarkhel, and finally optimized by me. It was developed on MATLAB R2016b.
