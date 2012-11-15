15 November 2012

------------Differences relative to older versions of APS----------------

This version of APS differs from the previously released version in the
following ways.

1) In the old version, the Kriging parameter was arbitrarily raised and lowered
so that a subjectively sufficient number of points were sampled far from the
high likelihood region.  In this version, the Kriging paramter is set by the
subroutine set_kp(int*) (in gaussian_process.cpp) by performing a Gaussian
process prediction on the points that have already been sampled.  The Kriging
parameter is set to a value such that 68% of these predictions fall within
1-sigma of the true values.  This is an expensive procedure, and is thus
performed only every tenth time the code goes to write its output.  If you would
like to adjust the frequency with which this happens, the call to set_kp(int*)
is in the subroutine add_pt(double*,double,int) in likelihoodinterface.cpp

2) In the old version, whenever a point was sampled within 10% of the target
confidence limit, the code interrupted its usual sampling routine and spent some
user-defined number of iterations using gradient descent to try to find points
that were within the target confidence limit.  The new version still uses
gradient descent, however, it is less monomaniacal.  Whenever the code finds a
point that is within some threshold (defined by the paramter grat in
likelihoodinterface.cpp; this can be set through the parameter file), it creates
a structure known as a ``gradient wanderer.'' This is an indepenent search which
will use gradient descent to try to find points within the confidene limit. 
There are only a finite number of gradient wanderers that the code will allow at
a time (set by the parameter gwroom in likelihoodinterface.cpp; this can be
controlled using the #Nw flag in the parameter file).  As the code is searching,
it alternates between doing a ``vanilla'' APS search (without gradient descent)
and doing gradient descent on each of its active wanderers.  The alternation is
such that the code spends around half of its time doing the vanilla
search.  Thus the code does not sacrifice its robust search of all parameter
space in favor of finding points within the confidence limit.  Note: when the
code sets the Kriging parameter (see point (1) above), it only considers points
that were sampled according to vanilla APS.

------------------Basic users' guide----------------------

To compile this version of aps, edit the Makefile so that it can find your
LAPACK and BLAS libraries, the libraries that will allow Fortran to talk to C++
and your C++ compiler and type

     make aps4p5

you will now have an executable aps4p5 which can be called with a parameter file
e.g.

     ./aps4p5 params_cartoon.sav

params_cartoon.sav contains various user-specified inputs which are detailed in
that file.  Note that all keywords must begin with a # in order for the code to
register them

The code that actually reads the parameter file is in aps_interface.cpp.  This
code instantiates an object of the class likelihood which does most of the work
of APS.  The class likelihood is stored within the files likelihoodinterface.h
and likelihoodinterface.cpp.  Some features which you may want to change are:

The subroutine call_likelihood(double*) actually calls your likelihood function
on the provided point in parameter space.  You will definitely want to change
this.

The subroutine write_pts() writes the output files (detailed below).  It
occasionally checks to make sure that the kd tree storing your sampled points
has been constructed correctly.  This is an expensive procedure.  If you want to
change how often it happens, look for the code check_tree(-1).  (The -1 is a
flag that tells check_tree(int) to check the whole tree, rather than just the
piece of the tree that descends from the point specified by int).  
The code that controls the kd tree is in kd.h and kd.cpp.  You should not need 
to change any of that code as it only exists to store the points sampled by 
APS and perform the nearest neighbor searches necessary for the Gaussian 
process prediction.

The class likelihood contains a member object of the class gp.  This class
controls the Gaussian process underlying APS.  The code for it exists in
gaussian_process.h and gaussian_process.cpp (which I have tried to comment
usefully).  Subroutines to be aware of:

covariogram(double*,double*,double*,int) -- this is the subroutine that encodes
the function underlying the covariogram.  You may want to replace it if you are
unhappy with the exp(-d^2/2) covariogram that we have chosen to use.  Note: if
you do replace it, the gradient descent code requires that this routine also be
able to return the gradient of the covariogram in the third double pointer.  See
the code for an example.

user_predict(double*,double*,int,double*) -- this is the routine that other code
calls to predict the value of chisquared using the Gaussian process.

user_predict_gradient(double*,double*,int) -- this is the routine that the
gradient wanderers call to predict the value of the gradient of chisquared using
the Gaussian process.

------------------Output files--------------

In the parameter file, there is a keyword #outputfile.  Running APS will produce
a file with the specified name which consists of a list of the sampled points in
parameter space, their chisquared values, a flag (which is zero if they were
sampled by vanilla MCMC an 1 if they were sampled by gradient descent), and the
value of the kriging parameter when the point was written.  These outputs are
generated by the subrouitne write_pts() in likelihoodinterface.cpp.  It is 
important that the kriging parameter always be the last column in this file 
(if you decide to change the output) so that interrupted APS runs can be resumed
by the subroutine resume(char*), also in likelihoodinterface.cpp.

APS will also produce a file named timingfile.sav which will help you gauge the
performance of APS.  This file contains (by column) 

 >the name of the output file
 
 >the number of points sampled
 
 > the number of gradient wanderers currently active
 
 >the number of points inside the confidence limit that have been found by the
   gradient wanderers
   
 >the number of points whose chisquared values were improved by the gradient 
   wanderers
   
 >the number of good points (inside the confidence limit) found
 
 >the average number of seconds (per sampled point) spent on the Gaussian 
  process in vanilla APS (if this is too large, you should reduce the number of
  candidates by changing #Nc in the parameter file)
  
 >the number of times vanilla APS evaluated candidates
 
 >the average number of seconds spent adding a point to the list of sampled
   points by vanilla APS (i.e. the length of each call to your likelihood 
   function)
   
 >the number of times a point was added by vanilla APS (this should be one less
    than the number of Gaussian process calls because of when write_pts() is
    called and when the counters are incremented)
 
 >the average number of seconds spent in a call to a gradient wanderer
 
 >the number of calls made to a gradient wanderer
 
 >the number of nearest neighbors used by the Gaussian process
 
 >the number of candidates considered by vanilla APS
 
 >a diagnostic of the kd tree (1 if everything is fine; 0 if not)
 
 >the current value of the Kriging parameter
 
 >the minimum value of chisquared found so far
 
Note: if a run gets interrupted (or stops naturally and you would like to
restart it), just call aps4p5 with the #resumefile flag set to the name of the
output file you would like to read in.  If you do not change #outputfile, the
code will not overwrite the old outputfile.  It will just append the new points
to the end.  If you do change #outputfile, the old outputs will get rewritten
into the new output file (except that the Kriging parametr column will be the
new value for the Kriging paramter).

If you do not change #end and the file you are resuming contains that many
points, no new points will be sampled.

#resumefile currently has no way of treating gradient wanderers, so you will
start over again from the case of no active wanderers

--------------------------General notes------------------

Throughout the code you will see references to "nodes" and "node"-based
subroutines.  This corresponds to a feature that did not perform well and has
been disabled.  These routines should not be called.

If you have an problems or questions please feel free to 
email scottvalscott@gmail.com

If this code results in published results, please cite arXiv:1205.2708
