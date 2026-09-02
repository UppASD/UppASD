target_sources(${UppASD_LIB} PRIVATE reducedstencil.f90
dipolecommon.f90
dipolemanager.f90
energy.f90
ewaldmom.f90
hamiltoniandatatype.f90
hamiltoniandata.f90
hamiltonianactions.f90
hamiltonianactions_lite.f90
hamiltonianinit.f90
neighbourmap.f90
neighbourmap_sfc.f90
targetorder.f90
printhamiltonian.f90)

if(USE_FFTW)
   target_sources(${UppASD_LIB} PRIVATE cpu_fft_provider.f90 cpuconvolution.f90)
endif()
