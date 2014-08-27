LAPACK_LIB =-L/Users/noldor/physics/lapack-3.1.1/ -llapack
BLAS_LIB =-L/Users/noldor/physics/BLAS/ -lblas

ARPACK_LIB = -L/Users/noldor/physics/ARPACK/ -larpack

FORTRAN_LIB = -L/opt/local/lib/ -lf95 -lgfortran

CAMB_LIB = -L/Users/noldor/physics/CAMB_110419/camb/ -lcamb
CAMB_INC = -I/Users/noldor/physics/CAMB_110419/camb/

CFITSIO_LIB = -L/Users/noldor/physics/cfitsio/ -lcfitsio

WMAP_LIB = -L/Users/noldor/physics/WMAP7likelihood/ -lwmap7 -lpthread
WMAP_INC = -I/Users/noldor/physics/WMAP7likelihood/

R_PATH = /Users/noldor/physics/lib/

LIBRARIES = $(LAPACK_LIB) $(BLAS_LIB) $(ARPACK_LIB) $(FORTRAN_LIB)

INCLUDE =

WMAP_LIBRARIES = $(LIBRARIES) $(WMAP_LIB) $(CAMB_LIB) $(CFITSIO_LIB)

WMAP_INCLUDE = $(INCLUDE) $(CAMB_INC) $(WMAP_INC)

#do not use these compilers with omp
#gg = g++ -Wno-write-strings -O3
#ff = gfortran -O3

#use these compilers if you are going to use omp
gg = /opt/local/bin/g++-mp-4.8 -Wno-write-strings -O3 -fopenmp -DUSE_OPENMP -g
ff = /opt/local/bin/gfortran-mp-4.8 -O3 -g

containers.o: containers.cpp containers.h
	$(gg) -c containers.cpp

goto_tools.o: goto_tools.h goto_tools.cpp containers.o
	$(gg) -c goto_tools.cpp containers.o $(LIBRARIES) $(INCLUDE)

test_containers: containers.o test_containers.cpp goto_tools.o
	$(gg) -o test_containers test_containers.cpp containers.o goto_tools.o

kd.o: kd.cpp kd.h goto_tools.o containers.o
	$(gg) -c kd.cpp goto_tools.o containers.o $(LIBRARIES) $(INCLUDE) \
         -Wno-deprecated

kde.o: kde.cpp kde.h
	$(gg) -c kde.cpp

test_kd: test_kd.cpp kd.o
	$(gg) -o test_kd test_kd.cpp containers.o goto_tools.o kd.o

chisq.o: goto_tools.o chisq.h chisq.cpp containers.o kd.o
	$(gg) -c chisq.cpp goto_tools.o containers.o kd.o $(LIBRARIES) $(INCLUDE)

chisq_aps4p5.o: chisq.o
	$(gg) -D_WMAP7_ -c chisq_aps4p5.cpp goto_tools.o containers.o chisq.o \
	$(LIBRARIES) $(INCLUDE)

camb_wrapper_wmap.o:
	$(ff) -D_WMAP7_ -c camb_wrapper_wmap.F90 $(WMAP_LIBRARIES) \
	$(WMAP_INCLUDE)

wmap_wrapper.o:
	$(ff) -D_WMAP7_ -c wmap_wrapper.F90 $(WMAP_LIBRARIES) $(WMAP_INCLUDE)

eigen_wrapper.o: goto_tools.o eigen_wrapper.cpp eigen_wrapper.h containers.o
	$(gg) -c eigen_wrapper.cpp goto_tools.o $(LIBRARIES) \
	$(INCLUDE) -Wno-deprecated

mcmc_extractor.o: mcmc_extractor.h mcmc_extractor.cpp goto_tools.o kd.o kde.o
	$(gg) -c mcmc_extractor.cpp

aps_extractor.o: aps_extractor.h aps_extractor.cpp goto_tools.o
	$(gg) -c aps_extractor.cpp

aps_extract: aps_extraction_runner.cpp aps_extractor.o kd.o
	$(gg) -o aps_extract aps_extraction_runner.cpp \
	containers.o goto_tools.o aps_extractor.o kd.o

aps_extract_ellipse: aps_extraction_ellipse.cpp aps_extractor.o kd.o
	$(gg) -o aps_extract_ellipse aps_extraction_ellipse.cpp \
	containers.o goto_tools.o aps_extractor.o kd.o

mcmc.o: chisq.o mcmc.h mcmc.cpp eigen_wrapper.o mcmc_extractor.o kd.o
	$(gg) -c mcmc.cpp \
	goto_tools.o containers.o chisq.o eigen_wrapper.o mcmc_extractor.o \
        kd.o kde.o \
	$(LIBRARIES) $(INCLUDE)

extract: mcmc_extraction_runner.cpp mcmc_extractor.o eigen_wrapper.o kd.o
	$(gg) -o extract mcmc_extraction_runner.cpp containers.o kd.o \
	goto_tools.o mcmc_extractor.o eigen_wrapper.o kde.o \
        $(LIBRARIES) $(INCLUDE)

test_eigen: test_eigen.cpp eigen_wrapper.o
	$(gg) -o test_eigen test_eigen.cpp goto_tools.o containers.o \
	eigen_wrapper.o $(LIBRARIES) $(INCLUDE) -Wno-deprecated

gaussian_process.o: gaussian_process.cpp gaussian_process.h kd.o goto_tools.o \
eigen_wrapper.o containers.o
	$(gg) -c gaussian_process.cpp goto_tools.o kd.o eigen_wrapper.o \
	containers.o \
	$(LIBRARIES) $(INCLUDE)

aps.o: aps.h aps.cpp kd.o goto_tools.o eigen_wrapper.o gaussian_process.o \
chisq.o containers.o
	$(gg) -c aps.cpp goto_tools.o \
	eigen_wrapper.o kd.o gaussian_process.o \
	chisq.o containers.o $(LIBRARIES) \
	$(INCLUDE) -Wno-deprecated

ellipse: aps_runner_ellipses.cpp aps.o chisq.o
	$(gg) -o ellipse aps_runner_ellipses.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o aps.o $(LIBRARIES) $(INCLUDE)

ellipseBrute: aps_ellipse_brute.cpp aps.o chisq.o
	$(gg) -o ellipseBrute aps_ellipse_brute.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o aps.o $(LIBRARIES) $(INCLUDE)



integrate: integrate_boundary.cpp kd.o eigen_wrapper.o gaussian_process.o chisq.o
	$(gg) -o integrate integrate_boundary.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o $(LIBRARIES) $(INCLUDE)

mcmc_test: mcmc.o mcmc_test.cpp
	$(gg) -o mcmc_test mcmc_test.cpp \
	goto_tools.o containers.o eigen_wrapper.o chisq.o mcmc.o \
	mcmc_extractor.o kd.o \
	$(LIBRARIES) $(INCLUDE)

mcmcBrute: mcmc.o mcmc_brute.cpp
	$(gg) -o mcmcBrute mcmc_brute.cpp \
	goto_tools.o containers.o eigen_wrapper.o chisq.o mcmc.o \
	mcmc_extractor.o kd.o \
	$(LIBRARIES) $(INCLUDE)

apsWMAP: aps_runner_wmap7.cpp aps.o chisq_aps4p5.o camb_wrapper_wmap.o \
wmap_wrapper.o
	$(gg) -D_WMAP7_ -o apsWMAP aps_runner_wmap7.cpp \
	goto_tools.o containers.o wmap_wrapper.o camb_wrapper_wmap.o \
	kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o chisq_aps4p5.o aps.o \
	$(WMAP_LIBRARIES) $(WMAP_INCLUDE)

all:
	make test_containers
	make test_kd
	make test_eigen
	make ellipse
	make integrate
	make mcmc_test
	make extract
	make aps_extract
	make aps_extract_ellipse
	make ellipseBrute
	make mcmcBrute
	make apsWMAP


clean:
	rm *.o test_containers test_kd test_eigen ellipse \
	integrate mcmc_test extract aps_extract aps_extract_ellipse \
	ellipseBrute mcmcBrute apsWMAP
