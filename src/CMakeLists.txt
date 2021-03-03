cmake_minimum_required(VERSION 2.7)
project(scspin)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O2 -std=c++11 -march=native -Wall -fopenmp -pipe -DHEADERONLY=0")

#find_package(Eigen3 REQUIRED)
#find_package(RapidJSON REQUIRED)
#find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#endif()
#add_subdirectory(.)

include_directories(${EIGEN3_INCLUDE_DIR})

#set(SOURCE_FILES main.cpp SpinState.cpp SpinState.h MolecularState.cpp MolecularState.h SpinState.cpp SpinState.h SpinStateEvolver.cpp SpinStateEvolver.h Simulation.cpp Simulation.h RandomSampler.cpp RandomSampler.h KahanSum.cpp KahanSum.h Parameters.cpp Parameters.h SimulationOptions.cpp SimulationOptions.h MarkovStateParameters.cpp MarkovStateParameters.h sobol.h)
#add_executable(scspin ${SOURCE_FILES} SpinStateEvolver.cpp SpinStateEvolver.h Simulation.cpp exitSimulation.h KahanSum.cpp KahanSum.h Parameters.cpp Parameters.h sobol.cpp.o)

#set(SOURCE_FILES main.cpp SpinState.cpp SpinState.h MolecularState.cpp MolecularState.h SpinState.cpp SpinState.h SpinStateEvolver.cpp SpinStateEvolver.h Simulation.cpp Simulation.h RandomSampler.cpp RandomSampler.h KahanSum.cpp KahanSum.h Parameters.cpp Parameters.h SimulationOptions.cpp SimulationOptions.h MarkovStateParameters.cpp MarkovStateParameters.h )

set(SOURCE_FILES main.cpp SpinState.cpp SpinState.h MolecularState.cpp MolecularState.h SpinState.cpp SpinState.h SpinStateEvolver.cpp SpinStateEvolver.h Simulation.cpp Simulation.h RandomSampler.cpp RandomSampler.h KahanSum.cpp KahanSum.h Parameters.cpp Parameters.h SimulationOptions.cpp SimulationOptions.h MarkovStateParameters.cpp MarkovStateParameters.h sobol.cpp sobol.h)

add_executable(scspin ${SOURCE_FILES} SpinStateEvolver.cpp SpinStateEvolver.h Simulation.cpp Simulation.h KahanSum.cpp KahanSum.h Parameters.cpp Parameters.h sobol.h)
