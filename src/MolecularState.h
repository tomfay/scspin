//
// Created by Tom on 04/12/2017.
//

#ifndef SEMICLASSICAL_MOLECULARSTATE_H
#define SEMICLASSICAL_MOLECULARSTATE_H

#include <iostream>
#include <vector>

#ifdef HEADERONLY
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#else
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#endif

#include <random>
#include <fstream>

#include "RandomSampler.h"
#include "MarkovStateParameters.h"

using namespace std ;
using namespace Eigen ;

class MolecularState
{
public:

  // constructor for the molecular state without an hopping rate information
  MolecularState( const vector<ArrayXi> &site_indices_, const vector<ArrayXXd> &site_couplings_ );

  // default constructor
  MolecularState() ;

  // sets the hopping rates from an input std::vector of Eigen::MatrixXd
  void setHoppingRates( vector<MatrixXd> const &hopping_rates_in ) ;

  // print the markov state information
  void printMarkovStates() ;

  // sets the time_step and creates the transition matrices
  void setTimeStep( double const time_step_in ) ;

  // sets the time step based on num_steps set in constructor and total time for each step
  void setTimStepFromTotal( double const total_time );

  //runs a series of hopping attempts
  void runMarkovChain( int const num_steps ) ;

  // runs a Markov chain model to update coupling constants
  void updateCouplings( int const num_steps , ArrayXd &couplings ) ;
  void updateCouplings( ArrayXd &couplings ) ;

  // constructor the markov state parameters from MarkovStateParameters object
  MolecularState( const MarkovStateParameters &params );

  // reset the random sampler seed
  void resetSeed() ;

  // prints info
  void printInfo() ;

  // samples the initial states uniformly
  void sampleInitialStatesUniformly() ;

private:

  // the number of steps used in the dynamics
  int num_steps_ ;

  // the number of groups of independent sites undergoing hops e.g. methyl groups
  int num_groups_ ;

  // the number of states that each site can hop between e.g. 3 for a methyl group
  ArrayXi array_num_states_ ;

  // the number of sites (hfccs) hopped between in each group
  ArrayXi array_num_sites_ ;

  // the indices of the sites that are hopped between for each independent group of sites
  vector<ArrayXi> site_indices_ ;

  // site couplings for each
  vector<ArrayXXd> site_couplings_ ;

  // matrix of hopping rates between different sites for each group of independent sites
  vector<MatrixXd> hopping_rates_ ;

  // time_step for hopping dynamics
  double time_step_ ;

  // transition probability matrices T = exp( K * dt )
  vector<MatrixXd> transition_matrices_ ;

  // cumulative probabilities
  vector<MatrixXd> cumulative_probabilities_ ;

  // average couplings for the sites in different groups over a markov chain run
  vector<ArrayXd> average_couplings_ ;

  // average populations/probabilities for being in a given markov state for each group of states
  vector<ArrayXd> state_probabilities_ ;

  // creates the transitions matrices
  void constructTransitionMatrices() ;

  // creates cumulative probability matrices for checking if a hop occurs
  void constructCumulativeProbs() ;

  // the states of all the independent groups
  ArrayXi current_states_ ;

  // distributions for each group of states
  vector<uniform_int_distribution<int>> int_distributions_ ;

  // generator for random number generation
  mt19937 generator_ ;



  // random sampler for generating hopping probabilities
  RandomSampler random_sampler_ ;

  // attempt hopping for each group of states for the current state
  void attemptHops() ;

  // random number
  double random_number_ ;

  // counter for the markov states over a run
  vector<ArrayXi> state_counts_ ;

  // updates the state counts
  void updateStateCounts() ;

  // reset state counts
  void resetStateCounts() ;

  // calculates the average couplings based on state counts for a given group
  void updateAverageCouplings( int const i_group) ;
  ArrayXd calculateStateProbabilities( const int i_group ) ;
};



#endif //SEMICLASSICAL_MOLECULARSTATE_H
