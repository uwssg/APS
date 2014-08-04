21 July 2014

This is the code implementing the Active Parameter Search algorithm as described
in REFERENCE.  If you find it at all useful in your research, please cite
REFERENCE.

This readme file will provide an overview for how the code works.  More specific
documentation describing the classes and functions involved can be found in the
.h header files.

(that ``more specific documentation'' has not actually been written yet, but
will be pushed onto the git repository in the next few weeks; I apologize for
the delay -- Scott)

///////////////////////REQUIREMENTS////////////////////////

In order to run APS, you will have to install the linear algebra
libraries LAPACK, BLAS, and ARPACK.  Edit the INCLUDE and LIBRARY
declarations in the Makefile accordingly.

///////////////////////BASIC EXAMPLE CODE//////////////////

Users interesed in a rough-and-ready example can find it in
aps_runner_ellipses.cpp.  To compile this code simply type

        make ellipse

This will run APS on a multi-modal likelihood function in which each mode is an
ellipse in parameter space.  Running

        ./ellipse 22 5 3

will run APS on a 5-dimensional parameter space with 3 modes (the 22 is the seed
for the pseuod random number generator used in APS).

        ./ellipse 34 6 2

will run APS on a 6-dimensional parameter space with 2 modes.

//////////////////////OUTPUT/////////////////

These runs will create files

        timingFiles/ellipse_d[dimensionality]_c[number of modes]_timing.sav
        outputFiles/ellipse_d[dimensionality]_c[number of modes]_output.sav

The output.sav is the output from APS.  Each row contains the point sampled in
parameter space, its chisquared value, the mu and sigma predicted for that
Gaussian process (if the point was sampled with the unextended APS; these are
set to -2.0 if the point was sampled by some other search method) and an integer
flag indicating how the point was chosen (0=unextended APS; 1=bisection or
simplex minimization; 2=steps 1B-5B in the paper)

The timing.sav file stores information about how APS is performing as it runs. 
Each row represents a moment during the run of APS.  The columns are as follows

        number of points stored by APS
        
        number of times chisquared was called (these are not necessarily the
        same; the code has safeguards so that, for example, it will not store
        points that return invalid chisquared values)
        
        the amount of clock time spent in the chisquared function
        
        the average amount of clock time per call to the chisqured function
        
        the total amount of clock time elapsed since the APS run began
        
        this total amount of time divided by the number of calls to chisquared
        (effectively, the amount of time per call to chisquared, including
        overhead introduced by APS)
        
        The number of chisquared calls devoted to non-simplex searches
        
        The amount  of clock time spent on non-simplex searches
        
        The number of chisquared calls devoted to simplex searches
        
        The amount of clock time spent on simplex searches
        
        The amount of clock time spent optimizing the Gaussian process's hyper
        parameters
        
        The amount of clock time spent redesignign the kd_tree that drives the
        Guassian process's nearest neighbor search
        
        The minimum discovered value of chisquared
        
        chi^2_lim
        
        The hyper volume associated with each identified low-chisquared region
        
        The number of times steps 1A-3A were called
        
        The number of times steps 1B-5B were called
        
        A diagnostic to make sure that the unit spheres in steps 1C-4C are being
        used

//////////////////RUNNING THE WMAP7 CODE

(Note: this code was written with a version CAMB downloaded on April 19, 2011
and version 4.1 of the WMAP likelihood code; there is no guarantee that it will
run with the most up-to-date version of either code; the user may have to modify
the following source code that is charged with interfacing APS with WMAP:

    camb_wrapper_wmap.F90
    wmap_wrapper.F90
    the class wmap_likelihood in chisq_aps4p5.cpp

)

If the user has successfully downloaded and compiled both CAMB and the WMAP
likelihood code and made the appropriate modifications to the Makefile, she can
run APS on the WMAP data by compiling

    make apsWMAP

and then running

    ./apsWMAP seed n_init n_end

where seed is a seed for the random number generator (-1 to use the system
clock), n_init is the number of initial random samples to make, and n_end is the
total number of calls to chisquared to make

The source code for apsWMAP can be found in aps_runner_wmap7.cpp

/////////////LIST OF SOURCE FILES AND BASIC DESCRIPTION OF THEIR CONTENTS

(recall that more detailed documentation will be given in the .h header files)

containers.cpp -- the code presented in this package manipulates 1- and
2-dimensional arrays using classes defined in this file.  These classes are
defined such that, if you ask for data that does not exist (e.g. for vv[6] when
vv only has elements 0 through 5), they will throw an exception and try to tell
you where it occurred.  The available classes are:

        array_1d<T> -- a 1-dimensional array of type T
        array_2d<T> -- a 2-dimensional array of type T in which all rows have 
                       the same number of columns
        asymm_array_2d<T> -- a 2-dimensional array in which each row can have
                             a different number of columns

The template variables T have been set so that these classes are available for
bouth doubles and ints.  Elements in these arrays are set like

        x.set(0,2.0) -- sets the 0th element of x to 2.0
        aa.set(1,2,3.0) -- sets the 2nd column of the 1st row of aa to 3.0 (note
                           that these arrays are zero-indexed)

The values of elements in these classes are retrieved like

        z = xx.get_data(8) -- sets z equal to the 8th element of xx
        w = aa.get_data(4,3) -- sets w equal to the 3rd column of the 4th row
                                of aa

see containers.h for more detailed documentation.

aps.cpp -- this defines the class aps, which actually performs the APS searches
described in the paper.  It requires the user to provide it with a chisquared
function, and a covariogram.  Once the class is instantiated, it must first be
initalized using the aps.initialize() method.  Searches can then be performed
using aps.search()

chisq.cpp -- this file defines the class chisquared which provides the basic
interface between chisquared functions and APS.  If a user wishes to write her
own chisquared function, it must inherit from the class chisquared.  The user's
daughter class will need to include an operator () which accepts an
array_1d<double>&, representing a point in parameter space and returning a
double, representing the value of chisquared at that point.

chisq_aps4p5.cpp -- this file defines some specific daught classes of the
chisquared class from chisq.cpp.  Specifically, it defines the wmap_likelihood
class, which is how APS runs the WMAP 7 year likelihood function.  If you want
to run this option, you will have to install CAMB and the WMAP likelihood code
on your machine and change the Makefile accordingly.

camb_wrapper_wmap.F90 -- fortran code that interfaces the wmap_likelihood class
from chisq_aps4p5.cpp with CAMB.

wmap_wrapper.F90 -- fortran code that interfaces the wmap_likelihood class form
chisq_aps4p5.cpp with the WMAP 7 likelihood code.

gaussian_process.cpp -- this file defines the class gp.cpp which performs the
Gaussian process inferences for APS.

kd.cpp -- this file defines the class kd_tree, which stores the points sampled
by APS in a KD-tree for easier nearest neighbor searching.

eigen_wrapper.cpp -- this file wraps the linear algebra routines of LAPACK,
BLAS, and ARPACK into C++.

goto_tools.cpp -- this file provides some generally interesting utility
functions.  Most notably, goto_tools.h defines our pseudo random number
generator.

aps_extractor.cpp -- defines the class aps_extractor which takes an APS output
file and outputs both the Frequentist and the Bayesian parameter constraints
implied thereby.

mcmc.cpp -- defines the class mcmc which runs a Monte Carlo Markov Chain

mcmc_extractor.cpp -- defines the class mcmc_extractor which takes the outputs
of a Monte Carlo Markov chain and draws independent samples from them.

kde.cpp -- defines the class kde which takes the independent samples from
mcmc_extractor and uses kernel density estimation to convert them into Bayesian
credible limits on parameter space.
