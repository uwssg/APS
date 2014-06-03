LAPACK_LIB =-L/Users/noldor/physics/lapack-3.1.1/ -llapack
BLAS_LIB =-L/Users/noldor/physics/BLAS/ -lblas

ARPACK_LIB = -L/Users/noldor/physics/ARPACK/ -larpack

FORTRAN_LIB = -L/opt/local/lib/ -lf95 -lgfortran

R_PATH = /Users/noldor/physics/lib/

LIBRARIES = $(LAPACK_LIB) $(BLAS_LIB) $(ARPACK_LIB) $(FORTRAN_LIB) \

INCLUDE = 

#do not use these compilers with omp
#gg = g++ -Wno-write-strings -O3
#ff = gfortran -O3

#use these compilers if you are going to use omp
gg = /opt/local/bin/g++-mp-4.8 -Wno-write-strings -O3 -fopenmp -DUSE_OPENMP
ff = /opt/local/bin/gfortran-mp-4.8 -O3

containers.o: containers.cpp containers.h
	$(gg) -c containers.cpp

goto_tools.o: goto_tools.h goto_tools.cpp containers.o
	$(gg) -c goto_tools.cpp containers.o $(LIBRARIES) $(INCLUDE)

test_containers: containers.o test_containers.cpp goto_tools.o
	$(gg) -o test_containers test_containers.cpp containers.o goto_tools.o

hyper_cubes.o: hyper_cubes.h hyper_cubes.cpp goto_tools.o
	$(gg) -c hyper_cubes.cpp goto_tools.o containers.o

test_hyper_cubes: test_hyper_cubes.cpp hyper_cubes.o
	$(gg) -o test_hyper_cubes test_hyper_cubes.cpp \
	containers.o goto_tools.o hyper_cubes.o

kd.o: kd.cpp kd.h goto_tools.o containers.o
	$(gg) -c kd.cpp goto_tools.o containers.o $(LIBRARIES) $(INCLUDE) \
         -Wno-deprecated

test_kd: test_kd.cpp kd.o
	$(gg) -o test_kd test_kd.cpp containers.o goto_tools.o kd.o

chisq.o: goto_tools.o chisq.h chisq.cpp containers.o
	$(gg) -c chisq.cpp goto_tools.o containers.o $(LIBRARIES) $(INCLUDE)

chisq_aps4p5.o: chisq.o
	$(gg) -c chisq_aps4p5.cpp goto_tools.o containers.o chisq.o \
	$(LIBRARIES) $(INCLUDE)

#sne.o: chisq.o sne.h sne.cpp
#	$(gg) -c sne.cpp goto_tools.o chisq.o $(LIBRARIES) $(INCLUDE)

eigen_wrapper.o: goto_tools.o eigen_wrapper.cpp eigen_wrapper.h containers.o
	$(gg) -c eigen_wrapper.cpp goto_tools.o $(LIBRARIES) \
	$(INCLUDE) -Wno-deprecated

mcmc_extractor.o: mcmc_extractor.h mcmc_extractor.cpp goto_tools.o
	$(gg) -c mcmc_extractor.cpp

aps_extractor.o: aps_extractor.h aps_extractor.cpp goto_tools.o
	$(gg) -c aps_extractor.cpp

aps_extract: aps_extraction_runner.cpp aps_extractor.o
	$(gg) -o aps_extract aps_extraction_runner.cpp \
	containers.o goto_tools.o aps_extractor.o

mcmc.o: chisq.o mcmc.h mcmc.cpp eigen_wrapper.o mcmc_extractor.o
	$(gg) -c mcmc.cpp \
	goto_tools.o containers.o chisq.o eigen_wrapper.o mcmc_extractor.o \
	$(LIBRARIES) $(INCLUDE)

extract: mcmc_extraction_runner.cpp mcmc_extractor.o eigen_wrapper.o
	$(gg) -o extract mcmc_extraction_runner.cpp containers.o \
	goto_tools.o mcmc_extractor.o eigen_wrapper.o \
        $(LIBRARIES) $(INCLUDE)

gp_to_mcmc.o: gp_to_mcmc.cpp gp_to_mcmc.h chisq.o gaussian_process.o
	$(gg) -c gp_to_mcmc.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o $(LIBRARIES) $(INCLUDE)

test_eigen: test_eigen.cpp eigen_wrapper.o
	$(gg) -o test_eigen test_eigen.cpp goto_tools.o containers.o \
	eigen_wrapper.o $(LIBRARIES) $(INCLUDE) -Wno-deprecated

gaussian_process.o: gaussian_process.cpp gaussian_process.h kd.o goto_tools.o \
eigen_wrapper.o containers.o
	$(gg) -c gaussian_process.cpp goto_tools.o kd.o eigen_wrapper.o \
	containers.o \
	$(LIBRARIES) $(INCLUDE)

exoplanet.o: exoplanet.h exoplanet.cpp chisq.o
	$(gg) -c exoplanet.cpp goto_tools.o containers.o chisq.o \
	$(LIBRARIES) $(INCLUDE)

dummy_planet.o: dummy_planet.h dummy_planet.cpp chisq.o
	$(gg) -c dummy_planet.cpp goto_tools.o containers.o chisq.o \
	$(LIBRARIES) $(INCLUDE)

test_exoplanet: test_exoplanet.cpp exoplanet.o
	$(gg) -o test_exoplanet test_exoplanet.cpp \
	goto_tools.o containers.o chisq.o exoplanet.o \
	$(LIBRARIES) $(INCLUDE)

test_1_exoplanet: test_1_exoplanet.cpp exoplanet.o
	$(gg) -o test_1_exoplanet test_1_exoplanet.cpp \
	goto_tools.o containers.o chisq.o exoplanet.o \
	$(LIBRARIES) $(INCLUDE)

aps.o: aps.h aps.cpp kd.o goto_tools.o eigen_wrapper.o gaussian_process.o \
chisq.o containers.o
	$(gg) -c aps.cpp goto_tools.o \
	eigen_wrapper.o kd.o gaussian_process.o \
	chisq.o containers.o $(LIBRARIES) \
	$(INCLUDE) -Wno-deprecated

planet: aps_runner_planets.cpp aps.o exoplanet.o
	$(gg) -o planet aps_runner_planets.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o exoplanet.o aps.o \
	$(LIBRARIES) $(INCLUDE)

ellipse: aps_runner_ellipses.cpp aps.o chisq.o
	$(gg) -o ellipse aps_runner_ellipses.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o aps.o $(LIBRARIES) $(INCLUDE)

integrate: integrate_boundary.cpp aps.o chisq.o
	$(gg) -o integrate integrate_boundary.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o aps.o $(LIBRARIES) $(INCLUDE)

dummy: aps_runner_dummy_planets.cpp aps.o dummy_planet.o
	$(gg) -o dummy aps_runner_dummy_planets.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o dummy_planet.o aps.o \
	$(LIBRARIES) $(INCLUDE)

planet2: aps_runner_2_planets.cpp aps.o exoplanet.o
	$(gg) -o planet2 aps_runner_2_planets.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o exoplanet.o aps.o \
	$(LIBRARIES) $(INCLUDE)

planet5: aps_runner_all_planets.cpp aps.o exoplanet.o
	$(gg) -o planet5 aps_runner_all_planets.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o exoplanet.o aps.o \
	$(LIBRARIES) $(INCLUDE)

planetTest: aps_planet_test.cpp aps.o exoplanet.o
	$(gg) -o planetTest aps_planet_test.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o exoplanet.o aps.o \
	$(LIBRARIES) $(INCLUDE)

analyze_planets: goto_tools.o containers.o analyze_planets.cpp
	$(gg) -o analyze_planets analyze_planets.cpp \
	goto_tools.o containers.o 

aps_test: aps_tester.cpp aps.o
	$(gg) -o aps_test aps_tester.cpp goto_tools.o kd.o \
	eigen_wrapper.o gaussian_process.o chisq.o aps.o \
	$(LIBRARIES) $(INCLUDE)

aps_run: aps_runner.cpp aps.o
	$(gg) -o aps_run aps_runner.cpp goto_tools.o kd.o \
	eigen_wrapper.o gaussian_process.o chisq.o aps.o \
	containers.o \
	$(LIBRARIES) $(INCLUDE)

apsForceGibbs: aps_interface.cpp goto_tools.o eigen_wrapper.o \
likelihoodinterface.o kd.o gaussian_process.o \
chisq.o 
	$(gg) -o apsForceGibbs aps_interface.cpp goto_tools.o eigen_wrapper.o kd.o \
	 likelihoodinterface.o gaussian_process.o \
	 chisq.o  $(LIBRARIES) \
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

gp2mcmc: gp_to_mcmc_driver.cpp gp_to_mcmc.o mcmc.o
	$(gg) -o gp2mcmc gp_to_mcmc_driver.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o gp_to_mcmc.o mcmc.o \
	$(LIBRARIES) $(INCLUDE)

mcmc_test: mcmc.o mcmc_test.cpp
	$(gg) -o mcmc_test mcmc_test.cpp \
	goto_tools.o containers.o eigen_wrapper.o chisq.o mcmc.o \
	mcmc_extractor.o \
	$(LIBRARIES) $(INCLUDE)

validate: validate_gp2mcmc.cpp gp2mcmc
	$(gg) -o validate validate_gp2mcmc.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o gp_to_mcmc.o mcmc.o \
	$(LIBRARIES) $(INCLUDE)

posterior_cubes: posterior_from_hyper_cubes.cpp kd.o
	$(gg) -o posterior_cubes posterior_from_hyper_cubes.cpp \
	containers.o goto_tools.o kd.o 

posterior_spheres: posterior_from_hyper_spheres.cpp kd.o
	$(gg) -o posterior_spheres posterior_from_hyper_spheres.cpp \
	containers.o goto_tools.o kd.o

explore: test_exploration.cpp gaussian_process.o
	$(gg) -o explore test_exploration.cpp \
	goto_tools.o containers.o eigen_wrapper.o kd.o gaussian_process.o \
	$(LIBRARIES) $(INCLUDE)


clean:
	rm *.o aps_run aps_test test_containers test_kd planet dummy planet2 \
	planet3 planet4 planet5 mcmc_test gp_to_mcmc ellipse posterior_cubes \
	posterior_spheres
