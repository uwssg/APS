21 July 2014

This is the code implementing the Active Parameter Search algorithm as described
in REFERENCE.  If you find it at all useful in your research, please cite
REFERENCE.

This readme file will provide an overview for how the code works.  More specific
documentation describing the classes and functions involved can be found in the
.h header files.

Users interesed in a rough-and-ready example can find it in
aps_runner_ellipses.cpp.  To compile this code simply type

	make ellipse

This will run APS on a multi-modal likelihood function in which each mode is an
ellipse in parameter space.  Running

	./ellipse 22 5 3

will run APS on a 5-dimensional parameter space with 3 modes (the 22 is the seed
for the pseuod random number generator used in APS).

	./ellipse 34 6 2

will run APS ona  6-dimensional parameter space with 2 modes.


