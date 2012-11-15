#where is your lapack library
LAPACK =-L/Users/noldor/physics/lapack-3.1.1/ -llapack

#where is your blas library
BLAS =-L/Users/noldor/physics/BLAS/ -lblas

#flags so that c++ can call Fortran routines (i.e.lapack and blas)
FORTRAN = -L/opt/local/lib/g95/x86_64-apple-darwin10/4.2.4/ -lf95 \
-lgfortran

#flags to call CAMB if you are using the CMB case
CAMB = -I/Users/noldor/physics/CAMB_110419/ \
-L/Users/noldor/physics/CAMB_110419/ -lcamb

#flags to call CFITSIO if you are using the CMB case
CFITSIO = /Users/noldor/physics/cfitsio

#flags to call WMAP if you are using the CMB case
WMAP = -I/Users/noldor/physics/WMAP7likelihood/ \
-L/Users/noldor/physics/WMAP7likelihood -lwmap7

WMAPLIBS = -L$(CFITSIO) -lcfitsio -L/Users/noldor/physics/lapack-3.1.1 \
-llapack -L/Users/noldor/physics/BLAS -lblas -lpthread

#C++ compiler
gg = /usr/local/bin/g++ 

#Fortran compiler
ff = /usr/local/bin/gfortran 

goto_tools.o: goto_tools.h goto_tools.cpp
	$(gg) -c goto_tools.cpp $(LAPACK) $(BLAS) $(FORTRAN)

eigen_wrapper.o: goto_tools.o eigen_wrapper.cpp eigen_wrapper.h
	$(gg) -c eigen_wrapper.cpp goto_tools.o $(LAPACK) $(BLAS) \
	$(FORTRAN) -Wno-deprecated

camb_wrapper_wmap.o: camb_wrapper_wmap.F90
	$(ff) -c camb_wrapper_wmap.F90 $(CAMB)
	
wmap_wrapper.o: wmap_wrapper.F90
	$(ff) -c wmap_wrapper.F90 -I/Users/noldor/physics/WMAP7likelihood/ \
	$(WMAP) $(FORTRAN) $(WMAPLIBS) -lgfortran

kd.o: kd.cpp kd.h goto_tools.o
	$(gg) -c kd.cpp goto_tools.o $(LAPACK) $(BLAS) $(FORTRAN) \
         -Wno-deprecated

gaussian_process.o: gaussian_process.cpp gaussian_process.h kd.o goto_tools.o \
eigen_wrapper.o
	$(gg) -c gaussian_process.cpp goto_tools.o kd.o eigen_wrapper.o \
	$(LAPACK) $(FORTRAN) $(BLAS) $(ARPACK)

#use the compile command below if you want to use the CAMB and WMAP 7 
#likelihood function
#
#make sure that likelihood::call_likelihood() in likelihoodinterface.cpp
#uses the proper likelihood function
#
#likelihoodinterface.o: camb_wrapper_wmap.o likelihoodinterface.h \
likelihoodinterface.cpp kd.o goto_tools.o eigen_wrapper.o wmap_wrapper.o \
gaussian_process.o
#	$(gg) -D CMB -c likelihoodinterface.cpp goto_tools.o \
	eigen_wrapper.o camb_wrapper_wmap.o kd.o wmap_wrapper.o \
	gaussian_process.o \
	$(LAPACK) $(BLAS) $(FORTRAN) -Wno-deprecated

#use the compile command below if you only want to use the cartoon likelihood 
#function

likelihoodinterface.o: likelihoodinterface.h \
likelihoodinterface.cpp kd.o goto_tools.o eigen_wrapper.o gaussian_process.o
	$(gg) -c likelihoodinterface.cpp goto_tools.o \
	eigen_wrapper.o kd.o gaussian_process.o \
	$(LAPACK) $(BLAS) $(FORTRAN) -Wno-deprecated



aps4p5: aps_interface.cpp goto_tools.o eigen_wrapper.o \
likelihoodinterface.o kd.o gaussian_process.o
	$(gg) -o aps4p5 aps_interface.cpp goto_tools.o eigen_wrapper.o kd.o \
	 likelihoodinterface.o gaussian_process.o \
	 $(LAPACK) $(BLAS) $(FORTRAN) \
	 -Wno-deprecated
