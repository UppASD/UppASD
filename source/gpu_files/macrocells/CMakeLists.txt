target_sources(${UppASD_LIB} PRIVATE
        gpuMacroHamiltonianCalculations.hpp
        gpuMacroHamiltonianCalculations.cpp
        hamiltonianFactory.hpp
        hamiltonian.hpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
