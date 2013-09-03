3 September 2013

This is the APS code used in Daniel et al (arXiv:1205.2708), which we will 
hereafter refer to as reference [1].  Please cite that reference in any 
published work resulting from this code.

We present here the code with both the cartoon likelihood function from section
4 of reference [1].  This code is contained in the class udder_likelihood in the
files likelihoodinterface.cpp and likelihoodinterface.h

We also present the WMAP7 likelihood interface, represented by the class
wmap_likelihood and the additional files wmap_wrapper.F90 and
camb_wmap_wrapper.F90 (note that this wrapper was written for the April 2011
version of CAMB, and may require some re-writing to interface with more modern
versions of CAMB).

If you do not want to compile with CAMB, simply remove the two -D_WMAP7_ flags
from the compiler definitions in the Makefile.

To write your own likelihood function, you must write a new class that inherits
from the class likelihood_function in the file likelihoodinterface.h. 
You must declare a new instance of that likelihood function in aps_interface.cpp
and pass it to the constructor
likelihood(int,double*,double*,covariance_function*,likelihood_function*).  See
how the pointer lk is handled in aps_interface.cpp for an example.

Below, we will list the files contained in this software package and briefly
describe the contents of each.  Consult the individual files for more detailed
comments.

-------------------------
gaussian_process.h and gaussian_process.cpp

These contain the class gp which stores the data for APS and performs the
Gaussian process interpolation.  The also contain the functors for different
covariograms.  Currently, the code supports the neural network covariogram
(equation 4.29 of Rasmussen and Williams 2006) and the Gaussian covariogram
(equation 4.9 of Rasmussen and Williams 2006).

To incorporate a new covariogram, one must write a class which inherits from the
class covariance_function.  The covariogram is passed to gp by way of the
likelihood class.  This is handled in aps_interface.cpp (see how the pointer cv
is used).

-------------------------
-------------------------
eigen_wrapper.h and eigen_wrapper.cpp

These files contain linear algebra routines wrapped from LAPACK and BLAS to a
C++ interface.  Principally, 

invert_lapack(double **matin, double **min, int el, int verb) 

takes the el-by-el symmetric matrix **matin (the input matrix), inverts it, 
and stores the invers in **min. verb is just a verbosity flag.

-------------------------
-------------------------
goto_tools.h and goto_tools.o

These contain generally useful subroutines to sort lists, generate random 
numbers, etc.

-------------------------


--------------------------General notes------------------

Throughout the code you will see references to "nodes" and "node"-based
subroutines.  This corresponds to a feature that did not perform well and has
been disabled.  These routines should not be called.

If you have an problems or questions please feel free to 
email scottvalscott@gmail.com

If this code results in published results, please cite arXiv:1205.2708
