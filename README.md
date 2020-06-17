## ICBM—Integrated Combined Baseline Modification: An Algorithm for Segmented Baseline Estimation

Accelerograms are the primary source for characterizing strong ground motion. It is therefore of paramount interest to have high‐quality recordings free from any nonphysical contamination. Frequently, accelerograms are affected by baseline jumps and drifts, either related to the instrument and/or a major earthquake. In this work, I propose a correction method for these undesired baseline drifts based on segmented linear least squares. The algorithm operates on the integrated waveforms and combines all three instrument components to estimate a model that modifies the baseline to be at zero continuously. The procedure consists of two steps: first a suite of models with variable numbers of discontinuities is derived for all three instrument components. During this process, the number of discontinuities is reduced in a parsimonious way, for example, two very close discontinuities are merged into a single one. In the second step, the optimal model is selected on the basis of the Bayesian information criterion. 

See https://doi.org/10.1785/0220190134 for a detailed description of the theory and application examples of ICBM.

# Software Prerequisites

ICBM is written in C++ using the linear algebra library armadillo.

http://arma.sourceforge.net/

This header file can be directly included in C++ code 

or 

can be made available in Python through numpy and a cython interface (under construction).   


# C++ Prerequisites:
	
armadillo (tested with v.9.880)
	
Check http://arma.sourceforge.net/download.html for installation on your system
	
# Python Prerequisites:
	
Python 3: numpy, cython
		
# Installation

Install the header file in your C++ library path or link it in the compiler. If you link to the headerfile by yourself, in g++ compilation looks like

$ g++ program.cpp -o program -larmadillo -I/path/to/the/file/icbm.h -O3 
   
In python you have to build the package for ICBM once.
This is covered in the makefile (coming soon).

# DEFCON: DEFinition of CONfiguration

ICBM requires a few parameters, set in DEFCON.

DEFCON is simply a textfile listing the configuration parameters for ICBM.
You can alter the order of entries and change the number of whitespaces or tabs. The basic structure

name value
	
must remain, and parameter names are not alterable.

You may not have to alter most of these values, however some need to be adjusted for every new data set.

list of parameters:

Parameter | Data type | Description | Default value
--|--|--|--
Nsmin | integer | minimum number of segments to be evaluated, must be larger than 0 | 1
Nsmax | integer | maximum number of segments to be evaluated, Nsmax >= Nsmin | 6
thrshld | double precision float | minimum velocity drift of a jump to be kept in the baseline model | 1e-5
nmin | integer | minimum number of iterations | 500

# Output

The output of ICBM contains the corrected accelerograms with baseline discontinuities removed.

In C++ the traces are stored in

std::vector<std::vector>> (M, std:vector<double>(N))

where M is the number of traces (usually 1 or 3,l currently only three traces are supported) and N is the number of samples per trace.

Optionally, the correction model can be written out. This is of particular interest, if the baseline jump associated with static near-field displacement should not be corrected for. Due to its sequential definition, the ICBM model will still correct for the remaining jumps.

CHECK that!
