//
// Created by Tom on 04/12/2017.
//

#include "Simulation.h"

/** Constructor for the Simulation object.
 *
 * @param [in] nuc_couplings_1_in isotropic hyperfine couplings for radical 1
 * @param [in] nuc_couplings_2_in isotropic hyperfine couplings for radical 1
 * @param [in] exchange_coupling_in exchange coupling in the form H_12 = exchange_coupling * S_1.S_2
 * @param [in] singlet_rate_in singlet recombination rate
 * @param [in] triplet_rate_in triplet recombination rate
 * @param [in] omega_1_in zeeman frequency vector for radical 1
 * @param [in] omega_2_in zeeman frequency vector for radical 2
 */
Simulation::Simulation( const ArrayXd nuc_couplings_1_in,
                        const ArrayXd nuc_couplings_2_in,
                        const double exchange_coupling_in,
                        const double singlet_rate_in,
                        const double triplet_rate_in,
                        const Vector3d omega_1_in,
                        const Vector3d omega_2_in )
{
  // extracts the number of nuclei on each radical
  num_nuc_1_ = (int) nuc_couplings_1_in.size() ;
  num_nuc_2_ = (int) nuc_couplings_2_in.size() ;

  // constructs the SC propagator
  propagator_ = SpinStateEvolver(nuc_couplings_1_in, nuc_couplings_2_in , exchange_coupling_in ,
                                 singlet_rate_in, triplet_rate_in, omega_1_in, omega_2_in) ;

  // constructs the SC state with the correct number of nuclei on each radical.
  state_ = SpinState( num_nuc_1_ , num_nuc_2_) ;
}

/** Samples the inital state for an SC trajectory.
 *
 * Electron and nuclear spins are sampled from spherical distributions, T = S_1 x S_2 and unit = 1.
 */
void Simulation::sampleInitialState()
{
  // Sample the nuclear spins
  for (int k = 0 ;  k < num_nuc_1_ ; k++)
  {
    state_.nuc_spins_1_.col(k) = state_.nuc_spin_lengths_1_(k) * random_sampler_.sampleUnitSphere() ;
  }

  for (int k = 0 ;  k < num_nuc_2_ ; k++)
  {
    state_.nuc_spins_2_.col(k) = state_.nuc_spin_lengths_2_(k) * random_sampler_.sampleUnitSphere() ;
  }

  // Sample the electron spin variables S_i/2(0)
  state_.electron_spin_variables_.segment(0,3) = 0.5 * sqrt(0.75) * random_sampler_.sampleUnitSphere() ;
  state_.electron_spin_variables_.segment(3,3) = 0.5 * sqrt(0.75) * random_sampler_.sampleUnitSphere() ;
  Vector3d new_vec ;
  // sets (T_k1, T_k2, T_k3) to S_1k(0) * S_2(0),
  for (int k = 0 ; k < 3 ; k++)
  {
//    state_.setTensorRow( k , 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ) ;
//    state_.electron_spin_variables_.segment(6+(3*k),3) = 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ;
    new_vec = 2.0 * state_.electron_spin_variables_(k) * state_.getSpinVector(2) ;
    state_.setTensorRow( k , new_vec) ;
  }
  // sets the initial value of the unit operator/4 to 1/4
  state_.electron_spin_variables_(15) = 0.25 ;
}

/** Samples the inital state for an SC trajectory.
 *
 * @param [inout] state the state which is sampled.
 *
 * Electron and nuclear spins are sampled from spherical distributions, T = S_1 x S_2 and unit = 1.
 */
void Simulation::sampleInitialState( SpinState &state )
{
  // Sample the nuclear spins
  for (int k = 0 ;  k < state.num_nuc_spins_1_ ; k++)
  {
    state.nuc_spins_1_.col(k) = state.nuc_spin_lengths_1_(k) * random_sampler_.sampleUnitSphere() ;
  }

  for (int k = 0 ;  k < state.num_nuc_spins_2_ ; k++)
  {
    state.nuc_spins_2_.col(k) = state.nuc_spin_lengths_2_(k) * random_sampler_.sampleUnitSphere() ;
  }

  // Sample the electron spin variables S_i/2(0)
  state.electron_spin_variables_.segment(0,3) = 0.5 * sqrt(0.75) * random_sampler_.sampleUnitSphere() ;
  state.electron_spin_variables_.segment(3,3) = 0.5 * sqrt(0.75) * random_sampler_.sampleUnitSphere() ;
  Vector3d new_vec ;
  // sets (T_k1, T_k2, T_k3) to S_1k(0) * S_2(0),
  for (int k = 0 ; k < 3 ; k++)
  {
//    state_.setTensorRow( k , 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ) ;
//    state_.electron_spin_variables_.segment(6+(3*k),3) = 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ;
    new_vec = 2.0 * state.electron_spin_variables_(k) * state.getSpinVector(2) ;
    state.setTensorRow( k , new_vec) ;
  }
  // sets the initial value of the unit operator/4 to 1/4
  state.electron_spin_variables_(15) = 0.25 ;
}

/** Samples the inital state for an SC trajectory.
 *
 * @param [inout] state the state which is sampled.
 *
 * Electron and nuclear spins are sampled from spherical distributions, T = S_1 x S_2 and unit = 1.
 */
void Simulation::sampleInitialState( SpinState &state, RandomSampler &random_sampler )
{
  // Sample the nuclear spins
  for (int k = 0 ;  k < state.num_nuc_spins_1_ ; k++)
  {
    state.nuc_spins_1_.col(k) = state.nuc_spin_lengths_1_(k) * random_sampler.sampleUnitSphere() ;
  }

  for (int k = 0 ;  k < state.num_nuc_spins_2_ ; k++)
  {
    state.nuc_spins_2_.col(k) = state.nuc_spin_lengths_2_(k) * random_sampler.sampleUnitSphere() ;
  }

  // Sample the electron spin variables S_i(0)/2
  #if USE_SOBOL
    VectorXd electron_spins(6) ;
    random_sampler.sampleTwoSpheresSobol(electron_spins) ;
    #if USE_NEW_VEC_LENS
      state.electron_spin_variables_.segment(0,6) = 0.5 * 0.5 * electron_spins ;
    #else
      state.electron_spin_variables_.segment(0,6) = 0.5 * sqrt(0.75) * electron_spins ;
    #endif
  #else
    #if USE_NEW_VEC_LENS
      state.electron_spin_variables_.segment(0,3) = 0.5 * 0.5 * random_sampler.sampleUnitSphere() ;
      state.electron_spin_variables_.segment(3,3) = 0.5 * 0.5 * random_sampler.sampleUnitSphere() ;
    #else
      state.electron_spin_variables_.segment(0,3) = 0.5 * sqrt(0.75) * random_sampler.sampleUnitSphere() ;
      state.electron_spin_variables_.segment(3,3) = 0.5 * sqrt(0.75) * random_sampler.sampleUnitSphere() ;
    #endif
  #endif


  Vector3d new_vec ;
  // sets (T_k1, T_k2, T_k3) to S_1k(0) * S_2(0),
  for (int k = 0 ; k < 3 ; k++)
  {
//    state_.setTensorRow( k , 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ) ;
//    state_.electron_spin_variables_.segment(6+(3*k),3) = 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ;
    new_vec = 2.0 * state.electron_spin_variables_(k) * state.getSpinVector(2) ;
    state.setTensorRow( k , new_vec) ;
  }
  // sets the initial value of the unit operator/4 to 1/4
  state.electron_spin_variables_(15) = 0.25 ;
}


/** Samples the inital state for an SC trajectory.
 *
 * @param [inout] state the state which is sampled.
 *
 * Electron and nuclear spins are sampled from spherical distributions, T = S_1 x S_2 and unit = 1.
 */
void Simulation::sampleInitialState( SpinState &state, RandomSampler &random_sampler ,Array2d electron_vector_lengths)
{
  // Sample the nuclear spins
  for (int k = 0 ;  k < state.num_nuc_spins_1_ ; k++)
  {
    state.nuc_spins_1_.col(k) = state.nuc_spin_lengths_1_(k) * random_sampler.sampleUnitSphere() ;
  }

  for (int k = 0 ;  k < state.num_nuc_spins_2_ ; k++)
  {
    state.nuc_spins_2_.col(k) = state.nuc_spin_lengths_2_(k) * random_sampler.sampleUnitSphere() ;
  }
  double length_1 = electron_vector_lengths(0) ;
  double length_2 = electron_vector_lengths(1) ;
  // Sample the electron spin variables S_i(0)/2
  #if USE_SOBOL
  VectorXd electron_spins(6) ;
  random_sampler.sampleTwoSpheresSobol(electron_spins) ;
  state.electron_spin_variables_.segment(0,3) = 0.5 * length_1 * electron_spins.segment(0,3) ;
  state.electron_spin_variables_.segment(3,6) = 0.5 * length_2 * electron_spins.segment(3,6) ;
  #else
  state.electron_spin_variables_.segment(0,3) = 0.5 * length_1 * random_sampler.sampleUnitSphere() ;
  state.electron_spin_variables_.segment(3,6) = 0.5 * length_2 * random_sampler.sampleUnitSphere() ;
  #endif


  Vector3d new_vec ;
  // sets (T_k1, T_k2, T_k3) to S_1k(0) * S_2(0),
  for (int k = 0 ; k < 3 ; k++)
  {
//    state_.setTensorRow( k , 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ) ;
//    state_.electron_spin_variables_.segment(6+(3*k),3) = 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ;
    new_vec = 2.0 * state.electron_spin_variables_(k) * state.getSpinVector(2) ;
    state.setTensorRow( k , new_vec) ;
  }
  // sets the initial value of the unit operator/4 to 1/4
  state.electron_spin_variables_(15) = 0.25 ;
}

/** Samples the inital state for an SC trajectory.
 *
 * @param [inout] state the state which is sampled.
 *
 * Electron and nuclear spins are sampled from spherical distributions, T = S_1 x S_2 and unit = 1.
 */
void Simulation::sampleInitialStateUnitElectrons( SpinState &state, RandomSampler &random_sampler )
{
  // Sample the nuclear spins
  for (int k = 0 ;  k < state.num_nuc_spins_1_ ; k++)
  {
    state.nuc_spins_1_.col(k) = state.nuc_spin_lengths_1_(k) * random_sampler.sampleUnitSphere() ;
  }

  for (int k = 0 ;  k < state.num_nuc_spins_2_ ; k++)
  {
    state.nuc_spins_2_.col(k) = state.nuc_spin_lengths_2_(k) * random_sampler.sampleUnitSphere() ;
  }

  // Sample the electron spin variables S_i(0)/2
#if USE_SOBOL
  VectorXd electron_spins(6) ;
  random_sampler.sampleTwoSpheresSobol(electron_spins) ;
  state.electron_spin_variables_.segment(0,6) = 0.5 * electron_spins ;

#else
  state.electron_spin_variables_.segment(0,3) = 0.5 * random_sampler.sampleUnitSphere() ;
  state.electron_spin_variables_.segment(3,6) = 0.5 * random_sampler.sampleUnitSphere() ;
#endif


  Vector3d new_vec ;
  // sets (T_k1, T_k2, T_k3) to S_1k(0) * S_2(0),
  for (int k = 0 ; k < 3 ; k++)
  {
//    state_.setTensorRow( k , 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ) ;
//    state_.electron_spin_variables_.segment(6+(3*k),3) = 4.0 * state_.electron_spin_variables_(k) * state_.electron_spin_variables_.segment(3,3) ;
    new_vec = 2.0 * state.electron_spin_variables_(k) * state.getSpinVector(2) ;
    state.setTensorRow( k , new_vec) ;
  }
  // sets the initial value of the unit operator/4 to 1/4
  state.electron_spin_variables_(15) = 0.25 ;
}

/** Runs a semi-classical spin dynamics simulation of a radical pair system
 *
 * @param [in] num_samples number of monte carlo samples for sampling integral
 * @param [in] num_steps number of steps for time evolution
 * @param [in] time_step time step for the simulation
 *
 * Data is output to a file "out.dat"
 */
void Simulation::runSimulation( int const num_samples, int const num_steps, double const time_step )
{
  // set the time step and create propagator matrices in the propagator
  propagator_.setTimeStep( time_step ) ;
  propagator_.createEvolutionOperators() ;

  // data array for holding 1, P_S(0)P_S and P_S(0)P_T correlation functions
  ArrayXXd data_array = ArrayXXd::Zero(3, num_steps+1) ;


  // creates a kahan sum object to perform kahan sums
  KahanSum kahan_sum_shared(3,num_steps+1);

  RandomSampler random_sampler_shared ;
  SpinStateEvolver propagator_shared = propagator_ ;
  SpinState state_shared = state_ ;

  // Monte carlo sampling for integrals
  #pragma omp parallel shared(kahan_sum_shared, random_sampler_shared, data_array, state_shared, propagator_shared)
  {
    if (omp_get_thread_num() == 0)
    {
      cout << "Number of threads created: " << omp_get_num_threads() << endl ;
    }
    // private copies of the state_ and propagator_ variables are created to avoid things fucking up
    ArrayXXd trajectory_data = ArrayXXd::Zero(3, num_steps+1) ;
    SpinState state_private ;
    SpinStateEvolver propagator_private ;


    // copy the original variables into the private versions
    #pragma omp critical
    {
      state_private = state_shared ;
      propagator_private = propagator_shared ;
    };
    #pragma omp for
    for ( int k = 0; k < num_samples; k ++ )
    {
      // sample the initial state
      #pragma omp critical
      {
        sampleInitialState( state_private , random_sampler_shared) ;
      };

      // run a trajectory and output data to an arrat
      trajectory_data = runTrajectory( num_steps, time_step , state_private ,propagator_private );

      // add the array to the monte carlo average
      #pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data );
      };
//    data_array += trajectory_data ;
    }
  }
  // rescale the sum to given correct integrals
  data_array = ( 4.0 / ( (double) num_samples ) ) * data_array ;


  // output data to "out.dat"
  ofstream data_file ;
  data_file.open("out.dat") ;
  for (int i = 0 ; i < num_steps+1 ; i++)
  {
    data_file << ((double) i) * time_step << " "
              << data_array(0,i) << " "
              << data_array(1,i) << " "
              << data_array(2,i) << endl ;
  }
  data_file.close() ;

}

/** Runs a single SC trajectory outputting P_S(0) * (1(t), P_S(t), P_T(t))
 *
 * @param [in] num_steps number of times the propagator is applied, number of steps
 * @param [in] time_step size of time step for integration
 * @return array of P_S(0) * (1(t), P_S(t), P_T(t)) for each time step
 */
ArrayXXd Simulation::runTrajectory( int const num_steps, double const time_step )
{
  // set up a container for trajectory data
  ArrayXXd trajectory_data = ArrayXXd::Zero(3, num_steps+1) ;
  // this is 0.25 here because everything is multiplied by 4 later
  Array3d data_initial( state_.calculateSingletProbability() , state_.calculateSingletProbability() , state_.calculateSingletProbability() ) ;

  // run dynamics storing trajectory data
  for (int i = 0 ; i < num_steps ; i++ )
  {
    trajectory_data(0,i) = state_.calculateSurvivalProbability() ;
    trajectory_data(1,i) = state_.calculateSingletProbability()  ;
    trajectory_data(2,i) = state_.calculateTripletProbability() ;
    propagator_.evolveState(state_) ;
  }


  trajectory_data(0,num_steps) = state_.calculateSurvivalProbability() ;
  trajectory_data(1,num_steps) = state_.calculateSingletProbability()  ;
  trajectory_data(2,num_steps) = state_.calculateTripletProbability() ;

  // multiply stuff by the initial variables to get correct integrand
  for (int i = 0 ; i < num_steps+1 ; i++)
  {

    trajectory_data(0,i) = data_initial(0) * trajectory_data(0,i) ;
    trajectory_data(1,i) = data_initial(1) * trajectory_data(1,i) ;
    trajectory_data(2,i) = data_initial(2) * trajectory_data(2,i) ;
  }
  return trajectory_data ;
}

/** Runs a single SC trajectory outputting P_S(0) * (1(t), P_S(t), P_T(t))
 *
 * @param [in] num_steps number of times the propagator is applied, number of steps
 * @param [in] time_step size of time step for integration
 * @param [inout] state the state to be propagated
 * @param [in] propagator the propagator for the dynamics
 * @return array of P_S(0) * (1(t), P_S(t), P_T(t)) for each time step
 */
ArrayXXd Simulation::runTrajectory( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator )
{
  // set up a container for trajectory data
  ArrayXXd trajectory_data = ArrayXXd::Zero(4, num_steps+1) ;
  // this is 0.25 here because everything is multiplied by 4 later
  ArrayXd data_initial = ArrayXd::Zero(4) ;
  data_initial(0) = state.calculateSingletProbability() ;
  data_initial(1) = data_initial(0) ;
  data_initial(2) = state.calculateTripletProbability() ;
  data_initial(3) = data_initial(2) ;

  // run dynamics storing trajectory data
  for (int i = 0 ; i < num_steps ; i++ )
  {
    trajectory_data(0,i) = state.calculateSingletProbability() ;
    trajectory_data(1,i) = state.calculateTripletProbability()  ;
    trajectory_data(2,i) = trajectory_data(0,i) ;
    trajectory_data(3,i) = trajectory_data(1,i) ;
    propagator.evolveState(state) ;
  }

  trajectory_data(0,num_steps) = state.calculateSingletProbability() ;
  trajectory_data(1,num_steps) = state.calculateTripletProbability() ;
  trajectory_data(2,num_steps) = trajectory_data(0,num_steps) ;
  trajectory_data(3,num_steps) = trajectory_data(1,num_steps) ;

  // multiply stuff by the initial variables to get correct integrand
  for (int i = 0 ; i < num_steps+1 ; i++)
  {

    trajectory_data(0,i) = 4.0*data_initial(0) * trajectory_data(0,i) ;
    trajectory_data(1,i) = 4.0*data_initial(1) * trajectory_data(1,i) ;
    trajectory_data(2,i) = 4.0*data_initial(2) * trajectory_data(2,i) ;
    trajectory_data(3,i) = 4.0*data_initial(3) * trajectory_data(3,i) ;
  }
  return trajectory_data ;
}

/** Runs two simultaneous simulations with different parameters and calculates difference signals
 *
 * @param num_samples
 * @param num_steps
 * @param time_step
 * @param parameters_1
 * @param parameters_2
 */
void Simulation::runDifferenceSimulation( int const num_samples,
                                          int const num_steps,
                                          double const time_step,
                                          Parameters parameters_1,
                                          Parameters parameters_2 )
{

  // creates the propagators and empty state objects
  SpinStateEvolver propagator_1(parameters_1) ;
  SpinStateEvolver propagator_2(parameters_2) ;
  SpinState state_1( parameters_1.getNumNuc()(0) , parameters_1.getNumNuc()(1) ) ;
  SpinState state_2( parameters_2.getNumNuc()(0) , parameters_2.getNumNuc()(1) ) ;
  state_1.setNucSpinLengths(parameters_1.getNucVecLengths(0),parameters_1.getNucVecLengths(1)) ;
  state_2.setNucSpinLengths(parameters_2.getNucVecLengths(0),parameters_2.getNucVecLengths(1)) ;

  // set the time step and create propagator matrices in the propagator
  propagator_1.setTimeStep( time_step ) ;
  propagator_1.createEvolutionOperators() ;
  propagator_2.setTimeStep( time_step ) ;
  propagator_2.createEvolutionOperators() ;

  // data array containing simulation data
  // Contains < P_x(0) P_y(t) > and < (P_x(0) P_y(t))^2 > for error calculations with x,y = S,T
  int num_data_vals = 7 * 5  ;
  ArrayXXd data_array = ArrayXXd::Zero( num_data_vals , num_steps+1 ) ;
  KahanSum kahan_sum_shared(num_data_vals,num_steps+1);
  int num_traj_complete = 0;
  RandomSampler random_sampler_shared ;


#pragma omp parallel \
  shared(kahan_sum_shared, data_array, state_1, propagator_1, state_2, propagator_2, random_sampler_shared, num_traj_complete)
  {
    if (omp_get_thread_num() == 0)
    {
      cout << "Number of threads created: " << omp_get_num_threads() << endl ;
    }

    // create private copies of the propagators and states
    SpinStateEvolver propagator_1_private;
    SpinStateEvolver propagator_2_private;
    SpinState state_2_private;
    SpinState state_1_private;

    #pragma omp critical
    {
      propagator_1_private = propagator_1;
      propagator_2_private = propagator_2;
//      propagator_1_private.printInfo() ;
//      propagator_2_private.printInfo() ;
      state_1_private = state_1;
      state_2_private = state_2;
//      state_1_private.printInfo() ;
//      state_2_private.printInfo() ;
    };

    #pragma omp barrier

    // creates a private copy of an array to store data from each trajectory
    ArrayXXd trajectory_data_private = ArrayXXd::Zero( num_data_vals, num_steps + 1 );


    // loops over samples
    if (omp_get_thread_num() == 0)
    {
      cout << "Starting sampling..."<< endl ;
    }
    #pragma omp for
    for ( int i = 0; i < num_samples; i ++ ) {
      // Sample initial state on state 1 and copy into state 2 - we assume the number of nuclei are the same on both.
      #pragma omp critical
      {
        sampleInitialState( state_1_private , random_sampler_shared);
      };
      state_2_private = state_1_private;

      // runs a difference trajectory and outputs the required data
      trajectory_data_private = runDifferenceTrajectory( num_steps, time_step, state_1_private, propagator_1_private,
                                                         state_2_private, propagator_2_private );

      // Adds data to the sum
      #pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
//        data_array = data_array + trajectory_data_private ;
        num_traj_complete = num_traj_complete + 1;
        checkpoint( num_traj_complete, num_samples, data_array, time_step );
      };
    }
  }

  // rescale the sum to given correct integrals
  data_array = ( 4.0 / ( (double) num_samples ) ) * data_array ;

  ofstream output_file ;
  output_file.open("diffout.dat") ;
  for (int k = 0 ; k < num_steps+1 ; k++)
  {
    output_file << time_step * k << " " << data_array.col(k).transpose() << endl ;
  }
  output_file.close() ;

}

/** Runs a difference trajectory - simultaneously propagates two states evolving with different propagators
 *
 * @param [in] num_steps number of time steps for evolution
 * @param [in] time_step time step for evolution
 * @param [in,out] state_1 state 1 to be propagated
 * @param [in] propagator_1 propagator for state 1
 * @param [in,out] state_2 state 2 to be propagated
 * @param [in] propagator_2 propagator for state 2
 * @return
 */
ArrayXXd Simulation::runDifferenceTrajectory( int const num_steps,
                                              double const time_step,
                                              SpinState &state_1,
                                              SpinStateEvolver &propagator_1,
                                              SpinState &state_2,
                                              SpinStateEvolver &propagator_2 )
{
  // An array of the initial
  Array2d data_initial( state_1.calculateSingletProbability(), state_1.calculateTripletProbability() ) ;

  // Array of trajectory
  int num_output_data = 7 ;

  // Trajectory data
  ArrayXXd trajectory_data = ArrayXXd::Zero( num_output_data , num_steps+1 ) ;

  // run dynamics storing trajectory data
  for (int i = 0 ; i < num_steps ; i++ )
  {
    // updates the trajectory data array
    updateDiffTrajData( state_1, state_2, trajectory_data, i ) ;

    // propagates the two states
    propagator_1.evolveState(state_1) ;
    propagator_2.evolveState(state_2) ;
  }

  // updates the trajectory data array with data from the final step
  updateDiffTrajData( state_1, state_2, trajectory_data, num_steps ) ;

  // returns the processed data with everything needed to calculate errors
  ArrayXXd output_data = processTrajData( trajectory_data , data_initial ) ;
  return output_data ;
}

/** Processes data from a trajectory with all data needed to estimate statistical errors
 *
 * @param trajectory_data data from the trajectory
 * @param initial_data initial data in the form (P_S(0), P_T(0))
 * @return
 */
ArrayXXd Simulation::processTrajData( ArrayXXd &trajectory_data, Array2d &initial_data )
{
  int num_traj_data = trajectory_data.rows();
  int traj_length = trajectory_data.cols() ;
  int num_output_data = num_traj_data * 5 ;

  ArrayXXd output_data = ArrayXXd::Zero(num_output_data , traj_length ) ;
  ArrayXd data_row =  ArrayXd::Zero(traj_length) ;
  for (int k = 0 ; k < num_traj_data ; k++ )
  {
    data_row = trajectory_data.row(k) ;
    output_data.row( (5 * k) ) = initial_data(0) * data_row ;
    output_data.row( (5 * k + 1) ) = initial_data(1) * data_row ;

    data_row = data_row * data_row ;
    output_data.row( (5 * k + 2) ) = (initial_data(0) * initial_data(0)) * data_row ;
    output_data.row( (5 * k + 3) ) = (initial_data(0) * initial_data(1)) * data_row ;
    output_data.row( (5 * k + 4) ) = (initial_data(1) * initial_data(1)) * data_row ;
  }
  return output_data;
}

/** Default constructor for Simulation class - does nothing.
 *
 */
Simulation::Simulation()
{

}

/** Updates the ith column of the trajectory data vector
 *
 * @param [in] state_1 SpinState for state 1
 * @param [in] state_2 SpinState for state 2
 * @param [inout] trajectory_data array of data for the diff trajectory
 * @param [in] i step to be updated/written
 */
void Simulation::updateDiffTrajData( SpinState &state_1, SpinState &state_2, ArrayXXd &trajectory_data, const int i )
{
  // writes (unit_1(t) -unit_2(t), unit_1(t), unit_2(t), P_S1(t), P_S2(t), P_T1(t), P_T2(t)) to the column of the array
  trajectory_data(0,i) = state_1.calculateSurvivalProbability() - state_2.calculateSurvivalProbability() ;
  trajectory_data(1,i) = state_1.calculateSurvivalProbability() ;
  trajectory_data(2,i) = state_2.calculateSurvivalProbability() ;
  trajectory_data(3,i) = state_1.calculateSingletProbability() ;
  trajectory_data(4,i) = state_2.calculateSingletProbability() ;
  trajectory_data(5,i) = state_1.calculateTripletProbability() ;
  trajectory_data(6,i) = state_2.calculateTripletProbability() ;
}

/** writes data out to a file
 *
 * @param [in] filename
 * @param [in] data_array
 * @param [in] num_samples
 */
void Simulation::writeDataOut( string const filename, ArrayXXd const &data_array, int const num_samples, double const time_step )
{
  // output filestream
  ofstream fileout ;
  fileout.open( filename ) ;

  //  create a copy of the data array scaled for the number of samples
  ArrayXXd data_array_out( data_array.rows(), data_array.cols()) ;
  data_array_out = (4.0 / ( (double) num_samples ) ) * data_array;

  // write out each column to a file
  for (int k = 0 ;  k < data_array_out.cols() ; k++)
  {
    fileout << ((double) k ) * time_step << " " << data_array_out.col(k).transpose() << endl ;
  }

  // close the output file
  fileout.close() ;
}

/** Checks to see how complete simulation is - if n*5% complete then outputs checkpoint file
 *
 * @param num_traj_completed
 * @param num_samples
 * @param data_array
 * @param time_step
 */
void Simulation::checkpoint( int const &num_traj_completed,
                             int const &num_samples,
                             ArrayXXd const &data_array,
                             double const &time_step )
{
  // checks to see if n* 5% finished
  if ((num_traj_completed > 0)&&(num_traj_completed%(num_samples/NUM_CHECKPOINTS)==0))
  {
    cout << (100)*num_traj_completed/(num_samples) << "% complete." << endl ;

    string filename = "checkpoint-"+to_string( (NUM_CHECKPOINTS*num_traj_completed) / num_samples )+".dat" ;

    writeDataOut( filename , data_array , num_traj_completed , time_step) ;
  }
}

void Simulation::runDiffSimWithHopping( int const num_samples,
                                        int const num_steps,
                                        double const time_step,
                                        Parameters parameters_1,
                                        Parameters parameters_2,
                                        MarkovStateParameters hopping_params_1,
                                        MarkovStateParameters hopping_params_2 )
{

  // creates the propagators and empty state objects
  SpinStateEvolver propagator_1(parameters_1) ;
  SpinStateEvolver propagator_2(parameters_2) ;
  SpinState state_1( parameters_1.getNumNuc()(0) , parameters_1.getNumNuc()(1) ) ;
  SpinState state_2( parameters_2.getNumNuc()(0) , parameters_2.getNumNuc()(1) ) ;
  state_1.setNucSpinLengths(parameters_1.getNucVecLengths(0),parameters_1.getNucVecLengths(1)) ;
  state_2.setNucSpinLengths(parameters_2.getNucVecLengths(0),parameters_2.getNucVecLengths(1)) ;

  // set the time step and create propagator matrices in the propagator
  propagator_1.setTimeStep( time_step ) ;
  propagator_1.createEvolutionOperators() ;
  propagator_2.setTimeStep( time_step ) ;
  propagator_2.createEvolutionOperators() ;

  // create the markov state objects
  MolecularState molecular_state_1( hopping_params_1 ) ;
  MolecularState molecular_state_2( hopping_params_2 ) ;


  // set the time steps for the molecular state hopping
  molecular_state_1.setTimStepFromTotal( time_step ) ;
  molecular_state_2.setTimStepFromTotal( time_step ) ;

  // data array containing simulation data
  // Contains < P_x(0) P_y(t) > and < (P_x(0) P_y(t))^2 > for error calculations with x,y = S,T
  int num_data_vals = 7 * 5  ;
  ArrayXXd data_array = ArrayXXd::Zero( num_data_vals , num_steps+1 ) ;
  KahanSum kahan_sum_shared(num_data_vals,num_steps+1);
  int num_traj_complete = 0;
  RandomSampler random_sampler_shared ;


#pragma omp parallel \
  shared(kahan_sum_shared, data_array, state_1, propagator_1, state_2, propagator_2, molecular_state_1, molecular_state_2, random_sampler_shared, num_traj_complete)
  {
    if (omp_get_thread_num() == 0)
    {
      cout << "Number of threads created: " << omp_get_num_threads() << endl ;
    }

    #pragma omp barrier

    // create private copies of the propagators and states
    SpinStateEvolver propagator_1_private;
    SpinStateEvolver propagator_2_private;
    SpinState state_2_private;
    SpinState state_1_private;
    MolecularState molecular_state_1_private ;
    MolecularState molecular_state_2_private ;

    // copy the propagator and state objects into thread private copies
    #pragma omp critical
    {
      propagator_1_private = propagator_1;
      propagator_2_private = propagator_2;
      state_1_private = state_1;
      state_2_private = state_2;
      molecular_state_1_private = molecular_state_1 ;
      molecular_state_2_private = molecular_state_2 ;
//      molecular_state_1_private.printInfo() ;
//      molecular_state_2_private.printInfo() ;
    };

    // reset the seed for the molecular state objects
    molecular_state_1_private.resetSeed() ;
    molecular_state_2_private.resetSeed() ;

    #pragma omp barrier

    // creates a private copy of an array to store data from each trajectory
    ArrayXXd trajectory_data_private ;
    #pragma omp critical
    {
      trajectory_data_private = ArrayXXd::Zero( num_data_vals, num_steps + 1 );
    };
    #pragma omp barrier
    // loops over samples
    #pragma omp for
    for ( int i = 0; i < num_samples; i++ )
    {
      // Sample initial state on state 1 and copy into state 2 - we assume the number of nuclei are the same on both.
      #pragma omp critical
      {
        sampleInitialState( state_1_private , random_sampler_shared );
      };
      state_2_private = state_1_private;

      // Sample the initial markov states - currently assumes initial state is uniformly distributed
      molecular_state_1.sampleInitialStatesUniformly() ;
      molecular_state_2.sampleInitialStatesUniformly() ;

      // runs a difference trajectory and outputs the required data
      trajectory_data_private = runDiffTrajWithHopping( num_steps, time_step, state_1_private, propagator_1_private,
                                                        state_2_private, propagator_2_private ,
                                                        molecular_state_1_private, molecular_state_2_private );

      // Adds data to the sum
      #pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
        num_traj_complete = num_traj_complete + 1;
        checkpoint( num_traj_complete, num_samples, data_array, time_step );
      };

    }
  }

  // rescale the sum to given correct integrals
  data_array = ( 4.0 / ( (double) num_samples ) ) * data_array ;

  ofstream output_file ;
  output_file.open("hoppingout.dat") ;
  for (int k = 0 ; k < num_steps+1 ; k++)
  {
    output_file << time_step * k << " " << data_array.col(k).transpose() << endl ;
  }
  output_file.close() ;

}

ArrayXXd Simulation::runDiffTrajWithHopping( int const num_steps,
                                             double const time_step,
                                             SpinState &state_1,
                                             SpinStateEvolver &propagator_1,
                                             SpinState &state_2,
                                             SpinStateEvolver &propagator_2,
                                             MolecularState &molecular_state_1,
                                             MolecularState &molecular_state_2 )
{
  // An array of the initial
  Array2d data_initial( state_1.calculateSingletProbability(), state_1.calculateTripletProbability() ) ;

  // Array of trajectory
  int num_output_data = 7 ;

  // Trajectory data
  ArrayXXd trajectory_data = ArrayXXd::Zero( num_output_data , num_steps+1 ) ;


  // An array of the couplings for each radical
  ArrayXd couplings_1 = propagator_1.getNucCouplings1();
  ArrayXd couplings_2 = propagator_1.getNucCouplings2();

  // run dynamics storing trajectory data
  for (int i = 0 ; i < num_steps ; i++ )
  {
    // updates the trajectory data array
    updateDiffTrajData( state_1, state_2, trajectory_data, i ) ;

    // update the hyperfine couplings for each radical
    molecular_state_1.updateCouplings( couplings_1 ) ;
    molecular_state_2.updateCouplings( couplings_2 ) ;
    propagator_1.setNucCouplings1( couplings_1 ) ;
    propagator_1.setNucCouplings2( couplings_2 ) ;
    propagator_2.setNucCouplings1( couplings_1 ) ;
    propagator_2.setNucCouplings2( couplings_2 ) ;

    // propagates the two states
    propagator_1.evolveState(state_1) ;
    propagator_2.evolveState(state_2) ;
  }

  // updates the trajectory data array with data from the final step
  updateDiffTrajData( state_1, state_2, trajectory_data, num_steps ) ;

  // returns the processed data with everything needed to calculate errors
  ArrayXXd output_data = processTrajData( trajectory_data , data_initial ) ;
  return output_data ;
}

/** Runs a simulation with different anisotropic couplings
 *
 * @param num_samples
 * @param num_steps
 * @param time_step
 * @param parameters_1
 */
void Simulation::runAnisoSimulation( int const num_samples,
                                          int const num_steps,
                                          double const time_step,
                                          Parameters parameters_1 )
{

  // creates the propagators and empty state objects
  SpinStateEvolver propagator_1(parameters_1) ;

  SpinState state_1( parameters_1.getNumNuc()(0) , parameters_1.getNumNuc()(1) ) ;
  state_1.setNucSpinLengths(parameters_1.getNucVecLengths(0),parameters_1.getNucVecLengths(1)) ;


  // set the time step and create propagator matrices in the propagator
  propagator_1.setTimeStep( time_step ) ;
  propagator_1.createEvolutionOperators() ;

  // data array containing simulation data
  // Contains < P_x(0) P_y(t) > and < (P_x(0) P_y(t))^2 > for error calculations with x,y = S,T
  int num_data_vals = 4  ;
  ArrayXXd data_array = ArrayXXd::Zero( num_data_vals , num_steps+1 ) ;
  ArrayXd yields_data = ArrayXd::Zero( 8 ) ;
  ArrayXd yields = ArrayXd::Zero( 4 ) ;
  KahanSum kahan_sum_shared(num_data_vals,num_steps+1);
  int num_traj_complete = 0;
  RandomSampler random_sampler_shared ;


#pragma omp parallel \
  shared(kahan_sum_shared, data_array, state_1, propagator_1, random_sampler_shared, num_traj_complete, yields, yields_data)
  {
    if (omp_get_thread_num() == 0)
    {
      cout << "Number of threads created: " << omp_get_num_threads() << endl ;
    }

    // create private copies of the propagators and states
    SpinStateEvolver propagator_1_private;
    SpinState state_1_private;

#pragma omp critical
    {
      propagator_1_private = propagator_1;
      state_1_private = state_1;
    };

#pragma omp barrier

    // creates a private copy of an array to store data from each trajectory
    ArrayXXd trajectory_data_private = ArrayXXd::Zero( num_data_vals, num_steps + 1 );


    // loops over samples
    if (omp_get_thread_num() == 0)
    {
      cout << "Starting sampling..."<< endl ;
    }
#pragma omp for
    for ( int i = 0; i < num_samples; i ++ ) {
      // Sample initial state on state 1 and copy into state 2 - we assume the number of nuclei are the same on both.
#pragma omp critical
      {
        sampleInitialState( state_1_private , random_sampler_shared);
      };

      // runs a difference trajectory and outputs the required data
      trajectory_data_private = runAnisoTrajectory( num_steps, time_step, state_1_private, propagator_1_private);

      // Adds data to the sum
#pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
//        data_array = data_array + trajectory_data_private ;
        num_traj_complete = num_traj_complete + 1;
        yields =   calculateYields(trajectory_data_private, propagator_1.getTimeStep(), propagator_1.getSingletRate(), propagator_1.getTripletRate()) ;
        yields_data.segment(0,4) = yields_data.segment(0,4)+ yields ;
        yields_data.segment(4,4) = yields_data.segment(4,4) +  yields * yields ;
        checkpoint( num_traj_complete, num_samples, data_array, time_step );
      };
    }
  }

  // rescale the sum to given correct integrals
  data_array = ( 1.0 / ( (double) num_samples ) ) * data_array ;

  ofstream output_file ;
  output_file.open("anisosimout.dat") ;
  output_file.precision(10) ;
  for (int k = 0 ; k < num_steps+1 ; k++)
  {
    output_file << scientific << time_step * k << " " << data_array.col(k).transpose() << endl ;
  }
  output_file.close() ;

  output_file.open("anisoyieldsout.dat") ;
  output_file.precision(10) ;
  yields_data = yields_data / ((double) num_samples) ;
  output_file << scientific << yields_data << endl ;
  output_file.close() ;

}

/** Runs a simulation with different anisotropic couplings
 *
 * @param num_samples
 * @param num_steps
 * @param time_step
 * @param parameters_1
 */
void Simulation::runSimulation( int const num_samples,
                                     int const num_steps,
                                     double const time_step,
                                     Parameters parameters_1 )
{

  // creates the propagators and empty state objects
  SpinStateEvolver propagator_1(parameters_1) ;
  SpinState state_1( parameters_1.getNumNuc()(0) , parameters_1.getNumNuc()(1) ) ;
  state_1.setNucSpinLengths(parameters_1.getNucVecLengths(0), parameters_1.getNucVecLengths(1)) ;

  // set the time step and create propagator matrices in the propagator
  propagator_1.setTimeStep( time_step ) ;
  propagator_1.createEvolutionOperators() ;

  // data array containing simulation data
  // Contains < P_x(0) P_y(t) > and < (P_x(0) P_y(t))^2 > for error calculations with x,y = S,T
  int num_data_vals = 4  ;
  ArrayXXd data_array = ArrayXXd::Zero( num_data_vals , num_steps+1 ) ;
  ArrayXd yields_data = ArrayXd::Zero( 8 ) ;
  ArrayXd yields = ArrayXd::Zero(4);
  KahanSum kahan_sum_shared(num_data_vals,num_steps+1);
  int num_traj_complete = 0;
  RandomSampler random_sampler_shared ;


#pragma omp parallel \
  shared(kahan_sum_shared, data_array, state_1, propagator_1, random_sampler_shared, num_traj_complete, yields, yields_data)
  {
    if (omp_get_thread_num() == 0)
    {
      cout << "Number of threads created: " << omp_get_num_threads() << endl ;
    }

    // create private copies of the propagators and states
    SpinStateEvolver propagator_1_private;
    SpinState state_1_private;

#pragma omp critical
    {
      propagator_1_private = propagator_1;
//      propagator_1_private.printInfo() ;
//      propagator_2_private.printInfo() ;
      state_1_private = state_1;
//      state_1_private.printInfo() ;
//      state_2_private.printInfo() ;
    };

#pragma omp barrier

    // creates a private copy of an array to store data from each trajectory
    ArrayXXd trajectory_data_private = ArrayXXd::Zero( num_data_vals, num_steps + 1 );


    // loops over samples
    if (omp_get_thread_num() == 0)
    {
      cout << "Starting sampling..."<< endl ;
    }
#pragma omp for
    for ( int i = 0; i < num_samples; i ++ ) {
      // Sample initial state on state 1 and copy into state 2 - we assume the number of nuclei are the same on both.
#pragma omp critical
      {
        sampleInitialState( state_1_private , random_sampler_shared);
      };

      // runs a difference trajectory and outputs the required data
      trajectory_data_private = runTrajectory( num_steps, time_step, state_1_private, propagator_1_private);

      // Adds data to the sum
#pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
//        data_array = data_array + trajectory_data_private ;
        num_traj_complete = num_traj_complete + 1;
        yields =   calculateYields(trajectory_data_private, propagator_1.getTimeStep(), propagator_1.getSingletRate(), propagator_1.getTripletRate()) ;
        yields_data.segment(0,4) = yields_data.segment(0,4)+ yields ;
        yields_data.segment(4,4) = yields_data.segment(4,4) +  yields * yields ;
        checkpoint( num_traj_complete, num_samples, data_array, time_step );
      };
    }
  }

  // rescale the sum to given correct integrals
  data_array = ( 1.0 / ( (double) num_samples ) ) * data_array ;

  ofstream output_file ;
  output_file.open("simout.dat") ;
  output_file.precision(10) ;
  for (int k = 0 ; k < num_steps+1 ; k++)
  {
    output_file << scientific << time_step * k << " " << data_array.col(k).transpose() << endl ;
  }
  output_file.close() ;
  output_file.open("yieldsout.dat") ;
  output_file.precision(10) ;
  yields_data = yields_data / ((double) num_samples) ;
  output_file << scientific << yields_data << endl ;
  output_file.close() ;

}

/** Runs a single SC trajectory outputting P_S(0) * (1(t), P_S(t), P_T(t))
 *
 * @param [in] num_steps number of times the propagator is applied, number of steps
 * @param [in] time_step size of time step for integration
 * @param [inout] state the state to be propagated
 * @param [in] propagator the propagator for the dynamics
 * @return array of P_S(0) * (1(t), P_S(t), P_T(t)) for each time step
 */
ArrayXXd Simulation::runAnisoTrajectory( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator )
{
  // set up a container for trajectory data
  ArrayXXd trajectory_data = ArrayXXd::Zero(4, num_steps+1) ;
  // this is 0.25 here because everything is multiplied by 4 later
  ArrayXd data_initial = ArrayXd::Zero(4) ;
  #if USE_NEW_VEC_LENS
    data_initial(0) = 9.0*state.calculateSingletProbability()-2.0 ;
    data_initial(1) = data_initial(0) ;
    data_initial(2) = 9.0*state.calculateTripletProbability()-6.0;
    data_initial(3) = data_initial(2) ;
  #else
    data_initial(0) = state.calculateSingletProbability() ;
    data_initial(1) = data_initial(0) ;
    data_initial(2) = state.calculateTripletProbability() ;
    data_initial(3) = data_initial(2) ;
  #endif

  // run dynamics storing trajectory data
  for (int i = 0 ; i < num_steps ; i++ )
  {
    trajectory_data(0,i) = state.calculateSingletProbability() ;
    trajectory_data(1,i) = state.calculateTripletProbability()  ;
    trajectory_data(2,i) = trajectory_data(0,i) ;
    trajectory_data(3,i) = trajectory_data(1,i) ;
    propagator.evolveStateAniso(state) ;
  }

  trajectory_data(0,num_steps) = state.calculateSingletProbability() ;
  trajectory_data(1,num_steps) = state.calculateTripletProbability()  ;
  trajectory_data(2,num_steps) = trajectory_data(0,num_steps) ;
  trajectory_data(3,num_steps) = trajectory_data(1,num_steps) ;

  // multiply stuff by the initial variables to get correct integrand
  for (int i = 0 ; i < num_steps+1 ; i++)
  {

    trajectory_data(0,i) = 4.0*data_initial(0) * trajectory_data(0,i) ;
    trajectory_data(1,i) = 4.0*data_initial(1) * trajectory_data(1,i) ;
    trajectory_data(2,i) = 4.0*data_initial(2) * trajectory_data(2,i) ;
    trajectory_data(3,i) = 4.0*data_initial(3) * trajectory_data(3,i) ;
  }
  return trajectory_data ;
}

ArrayXd Simulation::calculateYields( ArrayXXd const data , double const time_step, double const sing_rate, double const trip_rate)
{
  // array of [phi_SS, phi_TS, phi_ST , phi_TT ]
  ArrayXd yields = ArrayXd::Zero(4) ;
  for (int i = 0 ; i <4 ; i++)
  {
    yields(i) = integrateTrapezoidal( data.row(i), time_step) ;
  }
  yields(0) = sing_rate * yields(0) ;
  yields(1) = trip_rate * yields(1) ;
  yields(2) = (1.0/3.0)*sing_rate * yields(2) ;
  yields(3) = (1.0/3.0)*trip_rate * yields(3) ;
  return yields;
}

double Simulation::integrateTrapezoidal( ArrayXd const data_vals, double const spacing )
{
  double integral_val = 0.0 ;
  int num_data = (int) data_vals.size() ;
  integral_val = integral_val + spacing * (0.5*data_vals(0) + 0.5*data_vals(num_data-1) + data_vals.segment(1,num_data-2).sum()) ;
  return integral_val;
}

/** Runs a simulation with different anisotropic couplings
 *
 * @param num_samples
 * @param num_steps
 * @param time_step
 * @param parameters_1
 */
void Simulation::runAnisoSwSimulation( int const num_samples,
                                     int const num_steps,
                                     double const time_step,
                                     Parameters parameters_1 )
{

  // creates the propagators and empty state objects
  SpinStateEvolver propagator_1(parameters_1) ;

  SpinState state_1( parameters_1.getNumNuc()(0) , parameters_1.getNumNuc()(1) ) ;
  state_1.setNucSpinLengths(parameters_1.getNucVecLengths(0),parameters_1.getNucVecLengths(1)) ;


  // set the time step and create propagator matrices in the propagator
  propagator_1.setTimeStep( time_step ) ;
  propagator_1.createEvolutionOperators() ;

  // data array containing simulation data
  // Contains < P_x(0) P_y(t) > and < (P_x(0) P_y(t))^2 > for error calculations with x,y = S,T
  int num_data_vals = 4  ;
  ArrayXXd data_array = ArrayXXd::Zero( num_data_vals , num_steps+1 ) ;
  ArrayXd yields_data = ArrayXd::Zero( 8 ) ;
  ArrayXd yields = ArrayXd::Zero( 4 ) ;
  KahanSum kahan_sum_shared(num_data_vals,num_steps+1);
  int num_traj_complete = 0;
  RandomSampler random_sampler_shared ;


#pragma omp parallel \
  shared(kahan_sum_shared, data_array, state_1, propagator_1, random_sampler_shared, num_traj_complete, yields, yields_data)
  {
    if (omp_get_thread_num() == 0)
    {
      cout << "Number of threads created: " << omp_get_num_threads() << endl ;
    }

    // create private copies of the propagators and states
    SpinStateEvolver propagator_1_private;
    SpinState state_1_private;

#pragma omp critical
    {
      propagator_1_private = propagator_1;
      state_1_private = state_1;
    };

#pragma omp barrier

    // creates a private copy of an array to store data from each trajectory
    ArrayXXd trajectory_data_private = ArrayXXd::Zero( num_data_vals, num_steps + 1 );


    // loops over samples
    if (omp_get_thread_num() == 0)
    {
      cout << "Starting sampling..."<< endl ;
    }
#pragma omp for
    for ( int i = 0; i < num_samples; i ++ ) {
      // Sample initial state on state 1 and copy into state 2 - we assume the number of nuclei are the same on both.
#pragma omp critical
      {
        sampleInitialState( state_1_private , random_sampler_shared);
      };

      // runs a difference trajectory and outputs the required data
      trajectory_data_private = runAnisoSwTrajectory( num_steps, time_step, state_1_private, propagator_1_private);

      // Adds data to the sum
#pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
//        data_array = data_array + trajectory_data_private ;
        num_traj_complete = num_traj_complete + 1;
        yields =   calculateYields(trajectory_data_private, propagator_1.getTimeStep(), propagator_1.getSingletRate(), propagator_1.getTripletRate()) ;
        yields_data.segment(0,4) = yields_data.segment(0,4)+ yields ;
        yields_data.segment(4,4) = yields_data.segment(4,4) +  yields * yields ;
        checkpoint( num_traj_complete, num_samples, data_array, time_step );
      };
    }
  }

  // rescale the sum to given correct integrals
  data_array = ( 1.0 / ( (double) num_samples ) ) * data_array ;

  ofstream output_file ;
  output_file.open("anisoswsimout.dat") ;
  output_file.precision(10) ;
  for (int k = 0 ; k < num_steps+1 ; k++)
  {
    output_file << scientific << time_step * k << " " << data_array.col(k).transpose() << endl ;
  }
  output_file.close() ;

  output_file.open("anisoswyieldsout.dat") ;
  output_file.precision(10) ;
  yields_data = yields_data / ((double) num_samples) ;
  output_file <<scientific << yields_data << endl ;
  output_file.close() ;

}

/** Runs a single SC trajectory outputting P_S(0) * (1(t), P_S(t), P_T(t))
 *
 * @param [in] num_steps number of times the propagator is applied, number of steps
 * @param [in] time_step size of time step for integration
 * @param [inout] state the state to be propagated
 * @param [in] propagator the propagator for the dynamics
 * @return array of P_S(0) * (1(t), P_S(t), P_T(t)) for each time step
 */
ArrayXXd Simulation::runAnisoSwTrajectory( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator )
{
  // set up a container for trajectory data
  ArrayXXd trajectory_data = ArrayXXd::Zero(4, num_steps+1) ;
  // this is 0.25 here because everything is multiplied by 4 later
  ArrayXd data_initial = ArrayXd::Zero(4) ;
  data_initial(0) = state.calculateSingletProbability() ;
  data_initial(1) = data_initial(0) ;
  data_initial(2) = state.calculateTripletProbability() ;
  data_initial(3) = data_initial(2) ;

  propagator.updateOmegaEffsAniso(state) ;
  // run dynamics storing trajectory data
  for (int i = 0 ; i < num_steps ; i++ )
  {
    trajectory_data(0,i) = state.calculateSingletProbability() ;
    trajectory_data(1,i) = state.calculateTripletProbability()  ;
    trajectory_data(2,i) = trajectory_data(0,i) ;
    trajectory_data(3,i) = trajectory_data(1,i) ;
    propagator.evolveStateAnisoSw(state) ;
  }

  trajectory_data(0,num_steps) = state.calculateSingletProbability() ;
  trajectory_data(1,num_steps) = state.calculateTripletProbability()  ;
  trajectory_data(2,num_steps) = trajectory_data(0,num_steps) ;
  trajectory_data(3,num_steps) = trajectory_data(1,num_steps) ;

  // multiply stuff by the initial variables to get correct integrand
  for (int i = 0 ; i < num_steps+1 ; i++)
  {

    trajectory_data(0,i) = 4.0*data_initial(0) * trajectory_data(0,i) ;
    trajectory_data(1,i) = 4.0*data_initial(1) * trajectory_data(1,i) ;
    trajectory_data(2,i) = 4.0*data_initial(2) * trajectory_data(2,i) ;
    trajectory_data(3,i) = 4.0*data_initial(3) * trajectory_data(3,i) ;
  }
  return trajectory_data ;
}

/** Runs a simulation with different anisotropic couplings
 *
 * @param num_samples
 * @param num_steps
 * @param time_step
 * @param parameters_1
 */
void Simulation::runAnisoVectorOnlySimulation( int const num_samples,
                                     int const num_steps,
                                     double const time_step,
                                     Parameters parameters_1 )
{

  // creates the propagators and empty state objects
  SpinStateEvolver propagator_1(parameters_1) ;

  SpinState state_1( parameters_1.getNumNuc()(0) , parameters_1.getNumNuc()(1) ) ;
  state_1.setNucSpinLengths(parameters_1.getNucVecLengths(0),parameters_1.getNucVecLengths(1)) ;


  // set the time step and create propagator matrices in the propagator
  propagator_1.setTimeStep( time_step ) ;
  propagator_1.createEvolutionOperators() ;

  // data array containing simulation data
  // Contains < P_x(0) P_y(t) > and < (P_x(0) P_y(t))^2 > for error calculations with x,y = S,T
  int num_data_vals = 4  ;
  ArrayXXd data_array = ArrayXXd::Zero( num_data_vals , num_steps+1 ) ;
  ArrayXd yields_data = ArrayXd::Zero( 8 ) ;
  ArrayXd yields = ArrayXd::Zero( 4 ) ;
  KahanSum kahan_sum_shared(num_data_vals,num_steps+1);
  int num_traj_complete = 0;
  RandomSampler random_sampler_shared ;


#pragma omp parallel \
  shared(kahan_sum_shared, data_array, state_1, propagator_1, random_sampler_shared, num_traj_complete, yields, yields_data)
  {
    if (omp_get_thread_num() == 0)
    {
      cout << "Number of threads created: " << omp_get_num_threads() << endl ;
    }

    // create private copies of the propagators and states
    SpinStateEvolver propagator_1_private;
    SpinState state_1_private;

#pragma omp critical
    {
      propagator_1_private = propagator_1;
      state_1_private = state_1;
    };

#pragma omp barrier

    // creates a private copy of an array to store data from each trajectory
    ArrayXXd trajectory_data_private = ArrayXXd::Zero( num_data_vals, num_steps + 1 );


    // loops over samples
    if (omp_get_thread_num() == 0)
    {
      cout << "Starting sampling..."<< endl ;
    }
#pragma omp for
    for ( int i = 0; i < num_samples; i ++ ) {
      // Sample initial state on state 1 and copy into state 2 - we assume the number of nuclei are the same on both.
#pragma omp critical
      {
        sampleInitialState( state_1_private , random_sampler_shared);
      };

      // runs a difference trajectory and outputs the required data
      trajectory_data_private = runAnisoVectorOnlyTrajectory( num_steps, time_step, state_1_private, propagator_1_private);

      // Adds data to the sum
#pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
//        data_array = data_array + trajectory_data_private ;
        num_traj_complete = num_traj_complete + 1;
        yields =   calculateYields(trajectory_data_private, propagator_1.getTimeStep(), propagator_1.getSingletRate(), propagator_1.getTripletRate()) ;
        yields_data.segment(0,4) = yields_data.segment(0,4)+ yields ;
        yields_data.segment(4,4) = yields_data.segment(4,4) +  yields * yields ;
        checkpoint( num_traj_complete, num_samples, data_array, time_step );
      };
    }
  }

  // rescale the sum to given correct integrals
  data_array = ( 1.0 / ( (double) num_samples ) ) * data_array ;

  ofstream output_file ;
  output_file.open("anisosimveconlyout.dat") ;
  output_file.precision(10) ;
  for (int k = 0 ; k < num_steps+1 ; k++)
  {
    output_file << scientific << time_step * k << " " << data_array.col(k).transpose() << endl ;
  }
  output_file.close() ;

  output_file.open("anisoyieldsveconlyout.dat") ;
  output_file.precision(10) ;
  yields_data = yields_data / ((double) num_samples) ;
  output_file <<scientific << yields_data << endl ;
  output_file.close() ;

}

/** Runs a single SC trajectory outputting P_S(0) * (1(t), P_S(t), P_T(t))
 *
 * @param [in] num_steps number of times the propagator is applied, number of steps
 * @param [in] time_step size of time step for integration
 * @param [inout] state the state to be propagated
 * @param [in] propagator the propagator for the dynamics
 * @return array of P_S(0) * (1(t), P_S(t), P_T(t)) for each time step
 */
ArrayXXd Simulation::runAnisoVectorOnlyTrajectory( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator )
{
  // set up a container for trajectory data
  ArrayXXd trajectory_data = ArrayXXd::Zero(4, num_steps+1) ;
  // this is 0.25 here because everything is multiplied by 4 later
  ArrayXd data_initial = ArrayXd::Zero(4) ;
  data_initial(0) = state.calculateSingletProbabilityVectorOnly() ;
  data_initial(1) = data_initial(0) ;
  data_initial(2) = state.calculateTripletProbabilityVectorOnly() ;
  data_initial(3) = data_initial(2) ;

  // run dynamics storing trajectory data
  for (int i = 0 ; i < num_steps ; i++ )
  {
    trajectory_data(0,i) = state.calculateSingletProbabilityVectorOnly() ;
    trajectory_data(1,i) = state.calculateTripletProbabilityVectorOnly()  ;
    trajectory_data(2,i) = trajectory_data(0,i) ;
    trajectory_data(3,i) = trajectory_data(1,i) ;
    propagator.evolveStateAnisoVectors(state) ;
  }

  trajectory_data(0,num_steps) = state.calculateSingletProbabilityVectorOnly() ;
  trajectory_data(1,num_steps) = state.calculateTripletProbabilityVectorOnly()  ;
  trajectory_data(2,num_steps) = trajectory_data(0,num_steps) ;
  trajectory_data(3,num_steps) = trajectory_data(1,num_steps) ;

  // multiply stuff by the initial variables to get correct integrand
  double k_bar = 0.25*propagator.getSingletRate()+0.75*propagator.getTripletRate() ;
  double dt = propagator.getTimeStep() ;
  double decay_factor = exp( -dt * k_bar ) ;
  double population = 1.0 ;
  for (int i = 0 ; i < num_steps+1 ; i++)
  {

    trajectory_data(0,i) = population * 4.0*data_initial(0) * trajectory_data(0,i) ;
    trajectory_data(1,i) = population * 4.0*data_initial(1) * trajectory_data(1,i) ;
    trajectory_data(2,i) = population * 4.0*data_initial(2) * trajectory_data(2,i) ;
    trajectory_data(3,i) = population * 4.0*data_initial(3) * trajectory_data(3,i) ;
    population = decay_factor * population ;
  }


  return trajectory_data ;
}

/** Runs a simulation with different anisotropic couplings
 *
 * @param num_samples
 * @param num_steps
 * @param time_step
 * @param parameters_1
 */
void Simulation::runAnisoAltNormSimulation( int const num_samples,
                                     int const num_steps,
                                     double const time_step,
                                     Parameters parameters_1 )
{

  // creates the propagators and empty state objects
  SpinStateEvolver propagator_1(parameters_1) ;
  SpinState state_1( parameters_1.getNumNuc()(0) , parameters_1.getNumNuc()(1) ) ;
  state_1.setNucSpinLengths(parameters_1.getNucVecLengths(0),parameters_1.getNucVecLengths(1)) ;
//  cout << state_1.nuc_spin_lengths_1_ << endl << state_1.nuc_spin_lengths_2_ <<endl ;


  // set the time step and create propagator matrices in the propagator
  propagator_1.setTimeStep( time_step ) ;
  propagator_1.createEvolutionOperators() ;

  // data array containing simulation data
  // Contains < P_x(0) P_y(t) > and < (P_x(0) P_y(t))^2 > for error calculations with x,y = S,T
  int num_data_vals = 4  ;
  ArrayXXd data_array = ArrayXXd::Zero( num_data_vals , num_steps+1 ) ;
  ArrayXd yields_data = ArrayXd::Zero( 8 ) ;
  ArrayXd yields = ArrayXd::Zero( 4 ) ;
  KahanSum kahan_sum_shared(num_data_vals,num_steps+1);
  int num_traj_complete = 0;
  RandomSampler random_sampler_shared ;


#pragma omp parallel \
  shared(kahan_sum_shared, data_array, state_1, propagator_1, random_sampler_shared, num_traj_complete, yields, yields_data)
  {
    if (omp_get_thread_num() == 0)
    {
      cout << "Number of threads created: " << omp_get_num_threads() << endl ;
    }

    // create private copies of the propagators and states
    SpinStateEvolver propagator_1_private;
    SpinState state_1_private;

#pragma omp critical
    {
      propagator_1_private = propagator_1;
      state_1_private = state_1;
    };

#pragma omp barrier

    // creates a private copy of an array to store data from each trajectory
    ArrayXXd trajectory_data_private = ArrayXXd::Zero( num_data_vals, num_steps + 1 );


    // loops over samples
    if (omp_get_thread_num() == 0)
    {
      cout << "Starting sampling..."<< endl ;
    }
#pragma omp for
    for ( int i = 0; i < num_samples; i ++ ) {
      // Sample initial state on state 1 and copy into state 2 - we assume the number of nuclei are the same on both.
#pragma omp critical
      {
        sampleInitialState( state_1_private , random_sampler_shared);
      };

      // runs a difference trajectory and outputs the required data
      trajectory_data_private = runAnisoTrajectoryAltNorm( num_steps, time_step, state_1_private, propagator_1_private);

      // Adds data to the sum
#pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
//        data_array = data_array + trajectory_data_private ;
        num_traj_complete = num_traj_complete + 1;
        yields =   calculateYields(trajectory_data_private, propagator_1.getTimeStep(), propagator_1.getSingletRate(), propagator_1.getTripletRate()) ;
        yields_data.segment(0,4) = yields_data.segment(0,4)+ yields ;
        yields_data.segment(4,4) = yields_data.segment(4,4) +  yields * yields ;
        checkpoint( num_traj_complete, num_samples, data_array, time_step );
      };
    }
  }

  // rescale the sum to given correct integrals
  data_array = ( 1.0 / ( (double) num_samples ) ) * data_array ;

  ofstream output_file ;
  output_file.open("anisosimaltnormout.dat") ;
  output_file.precision(10) ;
  for (int k = 0 ; k < num_steps+1 ; k++)
  {
    output_file << scientific << time_step * k << " " << data_array.col(k).transpose() << endl ;
  }
  output_file.close() ;

  output_file.open("anisoyieldsaltnormout.dat") ;
  output_file.precision(10) ;
  yields_data = yields_data / ((double) num_samples) ;
  output_file <<scientific << yields_data << endl ;
  output_file.close() ;

}

/** Runs a single SC trajectory outputting P_S(0) * (1(t), P_S(t), P_T(t))
 *
 * @param [in] num_steps number of times the propagator is applied, number of steps
 * @param [in] time_step size of time step for integration
 * @param [inout] state the state to be propagated
 * @param [in] propagator the propagator for the dynamics
 * @return array of P_S(0) * (1(t), P_S(t), P_T(t)) for each time step
 */
ArrayXXd Simulation::runAnisoTrajectoryAltNorm( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator )
{
  // set up a container for trajectory data
  ArrayXXd trajectory_data = ArrayXXd::Zero(4, num_steps+1) ;
  // this is 0.25 here because everything is multiplied by 4 later
  ArrayXd data_initial = ArrayXd::Zero(4) ;
  #if USE_NEW_VEC_LENS
    data_initial(0) = 9.0*state.calculateSingletProbability()-2.0 ;
    data_initial(1) = data_initial(0) ;
    data_initial(2) = 9.0*state.calculateTripletProbability()-6.0;
    data_initial(3) = data_initial(2) ;
  #else
    data_initial(0) = state.calculateSingletProbability() ;
    data_initial(1) = data_initial(0) ;
    data_initial(2) = state.calculateTripletProbability() ;
    data_initial(3) = data_initial(2) ;
  #endif

  // run dynamics storing trajectory data
  for (int i = 0 ; i < num_steps ; i++ )
  {
    trajectory_data(0,i) = state.calculateSingletProbability() ;
    trajectory_data(1,i) = state.calculateTripletProbability()  ;
    trajectory_data(2,i) = trajectory_data(0,i) ;
    trajectory_data(3,i) = trajectory_data(1,i) ;
    propagator.evolveStateAnisoAltNorm(state) ;
  }

  trajectory_data(0,num_steps) = state.calculateSingletProbability() ;
  trajectory_data(1,num_steps) = state.calculateTripletProbability()  ;
  trajectory_data(2,num_steps) = trajectory_data(0,num_steps) ;
  trajectory_data(3,num_steps) = trajectory_data(1,num_steps) ;

  // multiply stuff by the initial variables to get correct integrand
  for (int i = 0 ; i < num_steps+1 ; i++)
  {

    trajectory_data(0,i) = 4.0*data_initial(0) * trajectory_data(0,i) ;
    trajectory_data(1,i) = 4.0*data_initial(1) * trajectory_data(1,i) ;
    trajectory_data(2,i) = 4.0*data_initial(2) * trajectory_data(2,i) ;
    trajectory_data(3,i) = 4.0*data_initial(3) * trajectory_data(3,i) ;
  }
  return trajectory_data ;
}

void Simulation::runAnisoNewLengthsSimulation( int const num_samples,
                                            int const num_steps,
                                            double const time_step,
                                            Parameters parameters_1 )
{

  // creates the propagators and empty state objects
  SpinStateEvolver propagator_1(parameters_1) ;

  SpinState state_1( parameters_1.getNumNuc()(0) , parameters_1.getNumNuc()(1) ) ;
  state_1.setNucSpinLengths(parameters_1.getNucVecLengths(0),parameters_1.getNucVecLengths(1)) ;


  // set the time step and create propagator matrices in the propagator
  propagator_1.setTimeStep( time_step ) ;
  propagator_1.createEvolutionOperators() ;

  // data array containing simulation data
  // Contains < P_x(0) P_y(t) > and < (P_x(0) P_y(t))^2 > for error calculations with x,y = S,T
  int num_data_vals = 4  ;
  ArrayXXd data_array = ArrayXXd::Zero( num_data_vals , num_steps+1 ) ;
  ArrayXd yields_data = ArrayXd::Zero( 8 ) ;
  ArrayXd yields = ArrayXd::Zero( 4 ) ;
  KahanSum kahan_sum_shared(num_data_vals,num_steps+1);
  int num_traj_complete = 0;
  RandomSampler random_sampler_shared ;
  Array2d ones_array = ArrayXd::Ones( 2 ) ;



#pragma omp parallel \
  shared(kahan_sum_shared, data_array, state_1, propagator_1, random_sampler_shared, num_traj_complete, yields, yields_data, ones_array)
  {
    if (omp_get_thread_num() == 0)
    {
      cout << "Number of threads created: " << omp_get_num_threads() << endl ;
    }

    // create private copies of the propagators and states
    SpinStateEvolver propagator_1_private;
    SpinState state_1_private;
    SpinState state_1_initial ;
    Vector3d new_vec ;
    double weight ;

#pragma omp critical
    {
      propagator_1_private = propagator_1;
      state_1_private = state_1;
      state_1_initial = state_1 ;
    };

#pragma omp barrier

    // creates a private copy of an array to store data from each trajectory
    ArrayXXd trajectory_data_private = ArrayXXd::Zero( num_data_vals, num_steps + 1 );

    // loops over samples
    if (omp_get_thread_num() == 0)
    {
      cout << "Starting sampling..."<< endl ;
    }

    //
#pragma omp for
    for ( int i = 0; i < num_samples; i ++ ) {


#pragma omp critical
      {
        sampleInitialStateUnitElectrons( state_1_initial , random_sampler_shared);
      };

      // first run with S_1 = S_2 = 0
      state_1_private = state_1_initial;
      state_1_private.electron_spin_variables_.segment(0,3) = 0.0* state_1_initial.electron_spin_variables_.segment(0,3) ;
      state_1_private.electron_spin_variables_.segment(3,6) = 0.0* state_1_initial.electron_spin_variables_.segment(3,6) ;
      for (int k = 0 ; k < 3 ; k++)
      {
        new_vec = 2.0 * state_1_private.electron_spin_variables_(k) * state_1_private.getSpinVector(2) ;
        state_1_private.setTensorRow( k , new_vec) ;
      }
      // run the trajectory
      trajectory_data_private = runAnisoTrajectoryNewLengths( num_steps, time_step, state_1_private, propagator_1_private);
      // multiply trajectory by appropriate weights
      weight = ((-0.8)*(-0.8)) ;
      trajectory_data_private = weight * trajectory_data_private ;

      // Adds data to the sum
      #pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
        //        data_array = data_array + trajectory_data_private ;
        yields =   calculateYields(trajectory_data_private, propagator_1.getTimeStep(), propagator_1.getSingletRate(), propagator_1.getTripletRate()) ;
        yields_data.segment(0,4) = yields_data.segment(0,4)+ yields ;
        yields_data.segment(4,4) = yields_data.segment(4,4) +  yields * yields ;

      };

      //  run with S_1 = sqrt(5/12), S_2 = 0
      state_1_private = state_1_initial;
      state_1_private.electron_spin_variables_.segment(0,3) = sqrt(5.0/12.0)* state_1_initial.electron_spin_variables_.segment(0,3) ;
      state_1_private.electron_spin_variables_.segment(3,6) = 0.0* state_1_initial.electron_spin_variables_.segment(3,6) ;
      for (int k = 0 ; k < 3 ; k++)
      {
        new_vec = 2.0 * state_1_private.electron_spin_variables_(k) * state_1_private.getSpinVector(2) ;
        state_1_private.setTensorRow( k , new_vec) ;
      }
      // run the trajectory
      trajectory_data_private = runAnisoTrajectoryNewLengths( num_steps, time_step, state_1_private, propagator_1_private);
      // multiply trajectory by appropriate weights
      trajectory_data_private = ((-0.8)*(1.8))*trajectory_data_private ;

      // Adds data to the sum
      #pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
        //        data_array = data_array + trajectory_data_private ;
        yields =   calculateYields(trajectory_data_private, propagator_1.getTimeStep(), propagator_1.getSingletRate(), propagator_1.getTripletRate()) ;
        yields_data.segment(0,4) = yields_data.segment(0,4)+ yields ;
        yields_data.segment(4,4) = yields_data.segment(4,4) +  yields * yields ;

      };

      //  run with S_1 = 0, S_2 = sqrt(5/12)
      state_1_private = state_1_initial;
      state_1_private.electron_spin_variables_.segment(0,3) = 0.0* state_1_initial.electron_spin_variables_.segment(0,3) ;
      state_1_private.electron_spin_variables_.segment(3,6) = sqrt(5.0/12.0)* state_1_initial.electron_spin_variables_.segment(3,6) ;
      for (int k = 0 ; k < 3 ; k++)
      {
        new_vec = 2.0 * state_1_private.electron_spin_variables_(k) * state_1_private.getSpinVector(2) ;
        state_1_private.setTensorRow( k , new_vec) ;
      }
      // run the trajectory
      trajectory_data_private = runAnisoTrajectoryNewLengths( num_steps, time_step, state_1_private, propagator_1_private);
      // multiply trajectory by appropriate weights
      trajectory_data_private = ((-0.8)*(1.8))*trajectory_data_private ;

      // Adds data to the sum
      #pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
        //        data_array = data_array + trajectory_data_private ;
        yields =   calculateYields(trajectory_data_private, propagator_1.getTimeStep(), propagator_1.getSingletRate(), propagator_1.getTripletRate()) ;
        yields_data.segment(0,4) = yields_data.segment(0,4)+ yields ;
        yields_data.segment(4,4) = yields_data.segment(4,4) +  yields * yields ;

      };

      //  run with S_1 = sqrt(5/12), S_2 = sqrt(5/12)
      state_1_private = state_1_initial;
      state_1_private.electron_spin_variables_.segment(0,3) = sqrt(5.0/12.0)* state_1_initial.electron_spin_variables_.segment(0,3) ;
      state_1_private.electron_spin_variables_.segment(3,6) = sqrt(5.0/12.0)* state_1_initial.electron_spin_variables_.segment(3,6) ;
      for (int k = 0 ; k < 3 ; k++)
      {
        new_vec = 2.0 * state_1_private.electron_spin_variables_(k) * state_1_private.getSpinVector(2) ;
        state_1_private.setTensorRow( k , new_vec) ;
      }
      // run the trajectory
      trajectory_data_private = runAnisoTrajectoryNewLengths( num_steps, time_step, state_1_private, propagator_1_private);
      // multiply trajectory by appropriate weights
      trajectory_data_private = ((1.8)*(1.8)) * trajectory_data_private ;

      // Adds data to the sum
      #pragma omp critical
      {
        kahan_sum_shared.addToSum( data_array, trajectory_data_private );
        //        data_array = data_array + trajectory_data_private ;
        num_traj_complete = num_traj_complete + 1;
        yields =   calculateYields(trajectory_data_private, propagator_1.getTimeStep(), propagator_1.getSingletRate(), propagator_1.getTripletRate()) ;
        yields_data.segment(0,4) = yields_data.segment(0,4)+ yields ;
        yields_data.segment(4,4) = yields_data.segment(4,4) +  yields * yields ;
        checkpoint( num_traj_complete, num_samples, data_array, time_step );
      };


    }
  }

  // rescale the sum to given correct integrals
  data_array = ( 1.0 / ( (double) num_samples ) ) * data_array ;

  ofstream output_file ;
  output_file.open("anisosimnewlengthsout.dat") ;
  output_file.precision(10) ;
  for (int k = 0 ; k < num_steps+1 ; k++)
  {
    output_file << scientific << time_step * k << " " << data_array.col(k).transpose() << endl ;
  }
  output_file.close() ;

  output_file.open("anisoyieldsnewlengthsout.dat") ;
  output_file.precision(10) ;
  yields_data = yields_data / ((double) num_samples) ;
  output_file <<scientific << yields_data << endl ;
  output_file.close() ;

}

/** Runs a single SC trajectory outputting P_S(0) * (1(t), P_S(t), P_T(t))
 *
 * @param [in] num_steps number of times the propagator is applied, number of steps
 * @param [in] time_step size of time step for integration
 * @param [inout] state the state to be propagated
 * @param [in] propagator the propagator for the dynamics
 * @return array of P_S(0) * (1(t), P_S(t), P_T(t)) for each time step
 */
ArrayXXd Simulation::runAnisoTrajectoryNewLengths( int const num_steps, double const time_step, SpinState &state, SpinStateEvolver &propagator )
{
  // set up a container for trajectory data
  ArrayXXd trajectory_data = ArrayXXd::Zero(4, num_steps+1) ;
  // this is 0.25 here because everything is multiplied by 4 later
  ArrayXd data_initial = ArrayXd::Zero(4) ;
  data_initial(0) = state.calculateSingletProbability() ;
  data_initial(1) = data_initial(0) ;
  data_initial(2) = state.calculateTripletProbability() ;
  data_initial(3) = data_initial(2) ;

  // run dynamics storing trajectory data
  for (int i = 0 ; i < num_steps ; i++ )
  {
    trajectory_data(0,i) = state.calculateSingletProbability() ;
    trajectory_data(1,i) = state.calculateTripletProbability()  ;
    trajectory_data(2,i) = trajectory_data(0,i) ;
    trajectory_data(3,i) = trajectory_data(1,i) ;
    propagator.evolveStateAnisoAltNorm(state) ;
  }

  trajectory_data(0,num_steps) = state.calculateSingletProbability() ;
  trajectory_data(1,num_steps) = state.calculateTripletProbability()  ;
  trajectory_data(2,num_steps) = trajectory_data(0,num_steps) ;
  trajectory_data(3,num_steps) = trajectory_data(1,num_steps) ;

  // multiply stuff by the initial variables to get correct integrand
  for (int i = 0 ; i < num_steps+1 ; i++)
  {

    trajectory_data(0,i) = 4.0*data_initial(0) * trajectory_data(0,i) ;
    trajectory_data(1,i) = 4.0*data_initial(1) * trajectory_data(1,i) ;
    trajectory_data(2,i) = 4.0*data_initial(2) * trajectory_data(2,i) ;
    trajectory_data(3,i) = 4.0*data_initial(3) * trajectory_data(3,i) ;
  }
  return trajectory_data ;
}
