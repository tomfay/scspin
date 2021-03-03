//
// Created by Tom on 04/12/2017.
//

#ifndef SEMICLASSICAL_SIMULATION_H
#define SEMICLASSICAL_SIMULATION_H

#define USE_SOBOL false
#define USE_NEW_VEC_LENS false
#define NUM_CHECKPOINTS 10


#include <fstream>
#include <iostream>
#include <string>
#include "SpinState.h"
#include "SpinStateEvolver.h"
#include "RandomSampler.h"
#include "KahanSum.h"
#include "omp.h"
#include "Parameters.h"
#include "MarkovStateParameters.h"
#include "MolecularState.h"

using namespace std ;

/** Class for running semi-classical spin dynamics simulations of radical pair reactions.
 *
 */
class Simulation
{
public:
  // Constructor for simulation object that constructs the propagator and spin state objects with the nput parameters
  Simulation( const ArrayXd nuc_couplings_1_in,
              const ArrayXd nuc_couplings_2_in,
              const double exchange_coupling_in,
              const double singlet_rate_in,
              const double triplet_rate_in,
              const Vector3d omega_1_in,
              const Vector3d omega_2_in );

  // Default constructor
  Simulation() ;

  // Runs a simulations with the input number of samples, number of time steps and time step size
  void runSimulation( int const num_samples, int const num_steps, double const time_step ) ;

  // Runs a difference simulation between two states run from the same initial configurations
  void runDifferenceSimulation( int const num_samples, int const num_steps, double const time_step, Parameters parameters_1 , Parameters parameters_2 );

  // Runs a difference simulation between two states run from the same initial configurations with hopping
  void runDiffSimWithHopping( int const num_samples, int const num_steps, double const time_step,
                              Parameters parameters_1 , Parameters parameters_2 ,
                              MarkovStateParameters hopping_params_1,
                              MarkovStateParameters hopping_params_2);

  // Runs a simulation with anisotropic hyperfine and electron-spin coupling interactions
  void runAnisoSimulation( int const num_samples,
                      int const num_steps,
                      double const time_step,
                      Parameters parameters_1);

  // Runs a simulation with anisotropic hyperfine and electron-spin coupling interactions with nuclear spin vectors frozen
  void runAnisoSwSimulation( int const num_samples,
                           int const num_steps,
                           double const time_step,
                           Parameters parameters_1);

  // Runs a simulations with the input number of samples, number of time steps and time step size
  void runSimulation( int const num_samples,
                 int const num_steps,
                 double const time_step,
                 Parameters parameters_1 );

  // Runs a simulation with the input number of samples, number of time steps and time step size with vector evolution only
  void runAnisoVectorOnlySimulation( int const num_samples,
                                                 int const num_steps,
                                                 double const time_step,
                                                 Parameters parameters_1 );

  // Runs a simulation with anisotropic hyperfine and electron-spin coupling interactions
  void runAnisoAltNormSimulation( int const num_samples,
                           int const num_steps,
                           double const time_step,
                           Parameters parameters_1);

  // Runs a simulation with anisotropic hyperfine and electron-spin coupling interactions
  void runAnisoNewLengthsSimulation( int const num_samples,
                                  int const num_steps,
                                  double const time_step,
                                  Parameters parameters_1);


private:

  // number of nuclear spins on each radical
  int num_nuc_1_ , num_nuc_2_ ;

  // The SC spin state that is used in simulations
  SpinState state_ ;

  // SC propagator for the spin state
  SpinStateEvolver propagator_ ;

  // Random sampler for generating initial states
  RandomSampler random_sampler_ ;

  // KahanSum object for accurately adding a large number of MC samples
  KahanSum kahan_sum_ ;

  // Runs a SC trajectory and returns some data from the trajectory, either performs with internal state/propagator
  // or an external state/propagator.
  ArrayXXd runTrajectory( int const num_steps , double const time_step) ;
  ArrayXXd runTrajectory( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator ) ;
  ArrayXXd runDifferenceTrajectory( int const num_steps, double const time_step,
                                    SpinState &state_1, SpinStateEvolver &propagator_1 ,
                                    SpinState &state_2, SpinStateEvolver &propagator_2 ) ;

  ArrayXXd runDiffTrajWithHopping( int const num_steps, double const time_step,
                                    SpinState &state_1, SpinStateEvolver &propagator_1 ,
                                    SpinState &state_2, SpinStateEvolver &propagator_2 ,
                                    MolecularState &molecular_state_1, MolecularState &molecular_state_2) ;
  ArrayXXd runAnisoTrajectory( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator ) ;

  ArrayXXd runAnisoSwTrajectory( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator );

  ArrayXXd runAnisoVectorOnlyTrajectory( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator );

  ArrayXXd runAnisoTrajectoryAltNorm( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator ) ;


  ArrayXXd runAnisoTrajectoryNewLengths( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator ) ;



  // Samples the initial spin state for the radical pair systems, either samples the private copy or an external one
  void sampleInitialState() ;
  void sampleInitialState( SpinState &state ) ;
  void sampleInitialState( SpinState &state, RandomSampler &random_sampler ) ;
  void sampleInitialState( SpinState &state, RandomSampler &random_sampler ,Array2d electron_vector_lengths);
  void sampleInitialStateUnitElectrons( SpinState &state, RandomSampler &random_sampler );

  // extracts data from two states
  ArrayXXd processTrajData( ArrayXXd &trajectory_data , Array2d &initial_data ) ;

  // Updates a trajectory data array with info from two states beign propagated simultaneously
  void updateDiffTrajData(SpinState &state_1 , SpinState &state_2, ArrayXXd &trajectory_data, int const i) ;

  // Writes data out to a file
  void writeDataOut( string const filename, ArrayXXd const &data_array , int const num_samples , double const time_step ) ;

  // Checkpointing - checks to see how far finished simulation is and if it is n * 5% finished then writes data out to a checkpoint file
  void checkpoint( int const &num_traj_completed, int const &num_samples, ArrayXXd const &data_array , double const &time_step) ;

  // calculate yields for a trajectory
  ArrayXd calculateYields( ArrayXXd const data , double const time_step, double const sing_rate, double const trip_rate);

  // integrates data with const spacing by trapezoidal rule
  double integrateTrapezoidal(ArrayXd const data_vals, double const spacing ) ;
};


#endif //SEMICLASSICAL_SIMULATION_H
