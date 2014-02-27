LAPACK_LIB =-L/Users/noldor/physics/lapack-3.1.1/ -llapack
BLAS_LIB =-L/Users/noldor/physics/BLAS/ -lblas

ARPACK_LIB = -L/Users/noldor/physics/ARPACK/ -larpack

FORTRAN_LIB = -L/opt/local/lib/ -lf95 -lgfortran

R_PATH = /Users/noldor/physics/lib/

LIBRARIES = $(LAPACK_LIB) $(BLAS_LIB) $(ARPACK_LIB) $(FORTRAN_LIB) \

INCLUDE = 

gg = g++ -Wno-write-strings -O3

ff = gfortran -O3


containers.o: containers.cpp containers.h
	$(gg) -c containers.cpp

goto_tools.o: goto_tools.h goto_tools.cpp containers.o
	$(gg) -c goto_tools.cpp containers.o $(LIBRARIES) $(INCLUDE)

test_containers: containers.o test_containers.cpp goto_tools.o
	$(gg) -o test_containers test_containers.cpp containers.o goto_tools.o

kd.o: kd.cpp kd.h goto_tools.o containers.o
	$(gg) -c kd.cpp goto_tools.o containers.o $(LIBRARIES) $(INCLUDE) \
         -Wno-deprecated

test_kd: test_kd.cpp kd.o
	$(gg) -o test_kd test_kd.cpp containers.o goto_tools.o kd.o

chisq_1311.o: goto_tools.o chisq_1311.h chisq_1311.cpp containers.o
	$(gg) -c chisq_1311.cpp goto_tools.o containers.o $(LIBRARIES) $(INCLUDE)

chisq_aps4p5.o: chisq_1311.o
	$(gg) -c chisq_aps4p5.cpp goto_tools.o containers.o chisq_1311.o \
	$(LIBRARIES) $(INCLUDE)

#sne.o: chisq_1311.o sne.h sne.cpp
#	$(gg) -c sne.cpp goto_tools.o chisq_1311.o $(LIBRARIES) $(INCLUDE)

eigen_wrapper.o: goto_tools.o eigen_wrapper.cpp eigen_wrapper.h containers.o
	$(gg) -c eigen_wrapper.cpp goto_tools.o $(LIBRARIES) \
	$(INCLUDE) -Wno-deprecated

test_eigen: test_eigen.cpp eigen_wrapper.o
	$(gg) -o test_eigen test_eigen.cpp goto_tools.o containers.o \
	eigen_wrapper.o $(LIBRARIES) $(INCLUDE) -Wno-deprecated

gaussian_process.o: gaussian_process.cpp gaussian_process.h kd.o goto_tools.o \
eigen_wrapper.o containers.o
	$(gg) -c gaussian_process.cpp goto_tools.o kd.o eigen_wrapper.o \
	containers.o \
	$(LIBRARIES) $(INCLUDE)

exoplanet.o: exoplanet.h exoplanet.cpp chisq_1311.o
	$(gg) -c exoplanet.cpp goto_tools.o containers.o chisq_1311.o \
	$(LIBRARIES) $(INCLUDE)

test_exoplanet: test_exoplanet.cpp exoplanet.o
	$(gg) -o test_exoplanet test_exoplanet.cpp \
	goto_tools.o containers.o chisq_1311.o exoplanet.o \
	$(LIBRARIES) $(INCLUDE)

test_1_exoplanet: test_1_exoplanet.cpp exoplanet.o
	$(gg) -o test_1_exoplanet test_1_exoplanet.cpp \
	goto_tools.o containers.o chisq_1311.o exoplanet.o \
	$(LIBRARIES) $(INCLUDE)

aps.o: aps.h aps.cpp kd.o goto_tools.o eigen_wrapper.o gaussian_process.o \
chisq_1311.o containers.o
	$(gg) -c aps.cpp goto_tools.o \
	eigen_wrapper.o kd.o gaussian_process.o \
	chisq_1311.o containers.o $(LIBRARIES) \
	$(INCLUDE) -Wno-deprecated

planet: aps_runner_planets.cpp aps.o exoplanet.o
	$(gg) -o planet aps_runner_planets.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq_1311.o exoplanet.o aps.o \
	$(LIBRARIES) $(INCLUDE)

planet2: aps_runner_2_planets.cpp aps.o exoplanet.o
	$(gg) -o planet2 aps_runner_2_planets.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq_1311.o exoplanet.o aps.o \
	$(LIBRARIES) $(INCLUDE)

analyze_planets: goto_tools.o containers.o analyze_planets.cpp
	$(gg) -o analyze_planets analyze_planets.cpp \
	goto_tools.o containers.o 

aps_test: aps_tester.cpp aps.o
	$(gg) -o aps_test aps_tester.cpp goto_tools.o kd.o \
	eigen_wrapper.o gaussian_process.o chisq_1311.o aps.o \
	$(LIBRARIES) $(INCLUDE)

aps_run: aps_runner.cpp aps.o
	$(gg) -o aps_run aps_runner.cpp goto_tools.o kd.o \
	eigen_wrapper.o gaussian_process.o chisq_1311.o aps.o \
	containers.o \
	$(LIBRARIES) $(INCLUDE)



apsForceGibbs: aps_interface.cpp goto_tools.o eigen_wrapper.o \
likelihoodinterface.o kd.o gaussian_process.o \
chisq_1311.o 
	$(gg) -o apsForceGibbs aps_interface.cpp goto_tools.o eigen_wrapper.o kd.o \
	 likelihoodinterface.o gaussian_process.o \
	 chisq_1311.o  $(LIBRARIES) \
	 $(INCLUDE) -Wno-deprecated


wp: write_params.cpp
	$(gg) -o wp write_params.cpp

#cartoon compilation statement
#aps2: aps_interface.cpp goto_tools.o eigen_wrapper.o gaussian_process.o \
likelihoodinterface.o kd.o 
#	$(gg) -o aps2 aps_interface.cpp goto_tools.o eigen_wrapper.o kd.o \
	 likelihoodinterface.o gaussian_process.o \
	 $(LAPACK) $(BLAS) $(FORTRAN) \
	 -Wno-deprecated

test: test_gp.cpp goto_tools.o eigen_wrapper.o \
likelihoodinterface.o kd.o camb_wrapper_wmap.o wmap_wrapper.o gaussian_process.o
	$(gg) -o test test_gp.cpp goto_tools.o eigen_wrapper.o kd.o \
	 likelihoodinterface.o camb_wrapper_wmap.o wmap_wrapper.o \
	 gaussian_process.o \
	 $(LAPACK) $(BLAS) $(FORTRAN) $(CAMB) $(WMAP) $(WMAPLIBS) \
	 -Wno-deprecated

clean:
	rm *.o aps_run aps_test test_containers test_kd
