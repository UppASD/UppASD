##############################################################################################
# CMAKE Configurations for UppASD
# (C) Anders Bergman 2019
#
# Acknowledgements: 
# https://github.com/SethMMorton/cmake_fortran_template  (for starting template)
# https://github.com/openmeeg/findmkl_cmake/blob/master/cmake/FindMKL.cmake (for FindMKL)
# https://github.com/libigl/eigen/blob/master/cmake/FindFFTW.cmake (for FindFFTW)
# 
##############################################################################################

target_sources(${UppASD_LIB} PRIVATE	gneb_driver.f90	mc_driver.f90	pt_driver.f90	sd_driver.f90	wl_driver.f90 sld_driver.f90 sx_driver.f90  ms_driver.f90 pyasd.f90 uppasd.f90)
###############################################################################################
## Add Sources in sub-directories
###############################################################################################
add_subdirectory(Clusters)
add_subdirectory(CoarseGraining)
add_subdirectory(Correlation)
add_subdirectory(Evolution)
add_subdirectory(Fields)
add_subdirectory(GNEB)
add_subdirectory(Hamiltonian)
add_subdirectory(Input)
add_subdirectory(KMC)
add_subdirectory(LSF)
add_subdirectory(Measurement)
add_subdirectory(MonteCarlo)
add_subdirectory(Parameters)
add_subdirectory(RNG)
add_subdirectory(SpinWaves)
add_subdirectory(System)
add_subdirectory(Tools)
add_subdirectory(Lattice)
add_subdirectory(Multiscale)
add_subdirectory(Third_party/Gridding)
##############################################################################################

