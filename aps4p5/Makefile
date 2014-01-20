
#directory where you can find the libraries for the LAPACK linear
#algebra package
LAPACK_LIB = -L/Users/noldor/physics/lapack-3.1.1/

#directory where you can find the libraries for the BLAS linear algebra package
BLAS_LIB = -L/Users/noldor/physics/BLAS/

#directory where you can find the libraries that link fortran to c++
FORTRAN_LIB = -L/opt/local/lib/

############below are directories and libraries for the WMAP7 likelihood

CAMB_LIB = -L/Users/noldor/physics/CAMB_110419/camb/ -lcamb
CAMB_INC = -I/Users/noldor/physics/CAMB_110419/camb/

WMAP7_LIB = -L/Users/noldor/physics/WMAP7likelihood/ -lpthread -lwmap7
WMAP7_INC = -I/Users/noldor/physics/WMAP7likelihood/

CFITSIO = -L/Users/noldor/physics/cfitsio -lcfitsio

###########################

LIBRARIES = $(LAPACK_LIB) \
$(BLAS_LIB) \
$(FORTRAN_LIB) \
$(CAMB_LIB) \
$(WMAP7_LIB) \
$(CFITSIO) \
-lblas -llapack -lf95 -lgfortran 

INCLUDES = $(CAMB_INC) $(WMAP_INC)

#C++ compiler
gg = /usr/local/bin/g++ -Wno-write-strings -D_WMAP7_

#Fortran compiler
ff = /usr/local/bin/gfortran -D_WMAP7_

goto_tools.o: goto_tools.h goto_tools.cpp
	$(gg) -c goto_tools.cpp $(LIBRARIES) $(INCLUDES)

eigen_wrapper.o: goto_tools.o eigen_wrapper.cpp eigen_wrapper.h
	$(gg) -c eigen_wrapper.cpp goto_tools.o $(LIBARIES) \
	$(INCLUDES) -Wno-deprecated

camb_wrapper_wmap.o: camb_wrapper_wmap.F90
	$(ff) -c camb_wrapper_wmap.F90 $(LIBRARIES) $(INCLUDES)
	
wmap_wrapper.o: wmap_wrapper.F90
	$(ff) -c wmap_wrapper.F90 -I/Users/noldor/physics/WMAP7likelihood/ \
	$(LIBRARIES) $(INCLUDES) -lgfortran

kd.o: kd.cpp kd.h goto_tools.o
	$(gg) -c kd.cpp goto_tools.o $(LIBRARIES) $(INCLUDES) \
         -Wno-deprecated

gaussian_process.o: gaussian_process.cpp gaussian_process.h kd.o goto_tools.o \
eigen_wrapper.o
	$(gg) -c gaussian_process.cpp goto_tools.o kd.o eigen_wrapper.o \
	$(LIBRARIES) $(INCLUDES)

chisq.o: chisq.h chisq.cpp wmap_wrapper.o camb_wrapper_wmap.o goto_tools.o
	$(gg) -c chisq.cpp goto_tools.o wmap_wrapper.o camb_wrapper_wmap.o \
	$(LIBRARIES) $(INCLUDES)
	
exoplanet.o: chisq.o exoplanet.h exoplanet.cpp kd.o gaussian_process.o \
eigen_wrapper.o
	$(gg) -c exoplanet.cpp goto_tools.o wmap_wrapper.o camb_wrapper_wmap.o \
	chisq.o kd.o eigen_wrapper.o gaussian_process.o $(LIBRARIES) $(INCLUDES)

test_exoplanet: exoplanet.o test_exoplanet.cpp kd.o gaussian_process.o \
eigen_wrapper.o
	$(gg) -o test_exoplanet test_exoplanet.cpp exoplanet.o goto_tools.o chisq.o \
	kd.o eigen_wrapper.o gaussian_process.o \
	wmap_wrapper.o camb_wrapper_wmap.o \
	$(LIBRARIES) $(INCLUDES)

test_1_exoplanet: exoplanet.o test_1_exoplanet.cpp
	$(gg) -o test_1_exoplanet test_1_exoplanet.cpp exoplanet.o goto_tools.o chisq.o \
	wmap_wrapper.o camb_wrapper_wmap.o \
	$(LIBRARIES) $(INCLUDES)

likelihoodinterface.o: camb_wrapper_wmap.o likelihoodinterface.h \
likelihoodinterface.cpp kd.o goto_tools.o eigen_wrapper.o wmap_wrapper.o \
gaussian_process.o chisq.o
	$(gg) -D CMB -c likelihoodinterface.cpp goto_tools.o \
	eigen_wrapper.o camb_wrapper_wmap.o kd.o wmap_wrapper.o \
	gaussian_process.o chisq.o \
	$(LIBRARIES) $(INCLUDES) -Wno-deprecated

aps4p5: aps_interface.cpp goto_tools.o eigen_wrapper.o \
likelihoodinterface.o kd.o gaussian_process.o wmap_wrapper.o \
camb_wrapper_wmap.o chisq.o exoplanet.o
	$(gg) -o aps4p5 aps_interface.cpp goto_tools.o eigen_wrapper.o kd.o \
	 likelihoodinterface.o gaussian_process.o \
	 camb_wrapper_wmap.o wmap_wrapper.o chisq.o exoplanet.o \
	 $(LIBRARIES) $(INCLUDES) \
	 -Wno-deprecated

udder: udder_test.cpp likelihoodinterface.o chisq.o
	$(gg) -o udder udder_test.cpp goto_tools.o eigen_wrapper.o kd.o \
	likelihoodinterface.o gaussian_process.o chisq.o camb_wrapper_wmap.o \
	wmap_wrapper.o \
	$(LIBRARIES) $(INCLUDE) \
	-Wno-deprecated

clean:
	rm aps4p5 test_exoplanet udder *.o
