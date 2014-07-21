21 July 2014

This is the code implementing the Active Parameter Search algorithm as described
in REFERENCE.  If you find it at all useful in your research, please cite
REFERENCE.

This readme file will provide an overview for how the code works.  More specific
documentation describing the classes and functions involved can be found in the
.h header files.

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
        
        The amount of clock time spent redesignign the kd_tree that drives the
        Guassian process's nearest neighbor search
        
        The amount of clock time spent optimizing the Gaussian process's hyper
        parameters
        
        The minimum discovered value of chisquared
        
        chi^2_lim
        
        The hyper volume associated with each identified low-chisquared region
        
        The number of times steps 1A-3A were called
        
        The number of times steps 1B-5B were called
        
        A diagnostic to make sure that the unit spheres in steps 1C-4C are being
        used
