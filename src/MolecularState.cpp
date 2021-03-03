//
// Created by Tom on 04/12/2017.
//

#include "MolecularState.h"

/** Constructor for the Molecular State class taking in the sites that are altered in the hopping and the values that
 * they change to.
 *
 * @param [in] site_indices_in a std::vector of Eigen::ArrayXi where each array contains the indices of the sites (nuclear
 *  spin indices) that are altered in each markov state for each independent markov state model (group).
 *  e.g. site_indices_in[1] = (i_11, i_12, ... , i_1n1)
 *       site_indices_in[2] = (i_21, i_22, ... , i_2n2)
 *       ...
 *       site_indices_in[N] = (i_N1, i_N2, ... , i_NnN)
 *
 * @param [in] site_couplings_in a std::vector of Eigen::ArrayXXd where each element corresponds to the different values of
 *  the hfccs for each site in each state. The columns of the arrays are the hfccs for each site defined in site_indices_in
 *  group k has Mk Markov states and nk sites
 *  e.g. site_couplings[k] = ( a(i_k1)_1 , a(i_k1)_2, a(i_k1)_3 ,... a(i_k1)_Mk )
 *                           ( a(i_k2)_1 , a(i_k2)_2, a(i_k2)_3 ,... a(i_k2)_Mk )
 *                           ...
 *                           ( a(i_knk)_1 , a(i_knk)_2, a(i_knk)_3 ,... a(i_knk)_Mk )
 *  here a(i)_n is the coupling constant for site i when the system is in state n.
 */
MolecularState::MolecularState( const vector<ArrayXi> &site_indices_in, const vector<ArrayXXd> &site_couplings_in )
{
  // extract the number of independent groups
  num_groups_ = (int) site_indices_in.size() ;
  // set the sizes of the arrays of the number of sites and states in each group
  array_num_sites_.resize(num_groups_) ;
  array_num_states_.resize(num_groups_) ;

  // resize the sizes of the site indices and site couplings arrays
  site_indices_.resize(num_groups_) ;
  site_couplings_.resize(num_groups_) ;
  hopping_rates_.resize(num_groups_) ;
  transition_matrices_.resize(num_groups_) ;
  cumulative_probabilities_.resize(num_groups_) ;
  state_counts_.resize(num_groups_);
  state_probabilities_.resize(num_groups_) ;
  average_couplings_.resize(num_groups_) ;

  // copy the input stuff to the member functions
  site_indices_ = site_indices_in ;
  site_couplings_ = site_couplings_in ;

  for (int i = 0 ; i < num_groups_ ; i++)
  {
    // extracts the number of sites for a given group
    array_num_sites_(i) = (int) site_couplings_[i].rows() ;

    // extracts the number of states for a given group of independent states
    array_num_states_(i) = (int) site_couplings_[i].cols() ;

    // fills the hopping rates matrix and transition matrix for a group of independent sites with zeros
    hopping_rates_[i] = ArrayXXd::Zero(array_num_states_(i), array_num_states_(i)) ;
    transition_matrices_[i] = ArrayXXd::Zero(array_num_states_(i), array_num_states_(i)) ;
    cumulative_probabilities_[i] = ArrayXXd::Zero(array_num_states_(i), array_num_states_(i)) ;
    state_counts_[i] = ArrayXi::Zero(array_num_states_(i)) ;
    state_probabilities_[i] = ArrayXd::Zero(array_num_states_(i)) ;
    average_couplings_[i] = ArrayXd::Zero(array_num_sites_(i)) ;
  }

  time_step_ = 0.0 ;

  // setting up the random number stuff randomly seeds the generator
  random_device seed_generator ;
  generator_.seed(seed_generator()) ;

  // sets up the distributions
  int_distributions_.resize(num_groups_) ;
  for (int i = 0 ; i< num_groups_ ; i++)
  {
    int_distributions_[i] = uniform_int_distribution<int>(0,array_num_states_(i)-1) ;
  }

  // sets up the array of current states
  current_states_ = ArrayXi::Zero(num_groups_) ;
}

/** Prints information about each of the Markov state groups (a group is a set of sites that hop independently of the
 * rest of the molecule), this includes the number of sites and number of states in each group and the site indices for
 * each group and the couplings for each site in each state.
 */
void MolecularState::printMarkovStates()
{
  for (int i = 0 ; i< num_groups_ ; i++ )
  {
    cout << "Group: " << i << endl ;
    cout << "Number of sites:  " << array_num_sites_(i) << endl ;
    cout << "Number of states: " << array_num_states_(i) << endl ;
    cout << "Site indices: " << site_indices_[i].transpose() << endl ;
    cout << "Site Couplings: " << endl ;
    cout << site_couplings_[i] << endl << endl ;
  }
}

/** Sets the hopping rates between markov states for each group of independent markov state
 *
 * @param [in] hopping_rates_in std::vector of Eigen::MatrixXd where each element is a hopping rate matrix between
 *  different states in a group of independent sites.
 */
void MolecularState::setHoppingRates( vector<MatrixXd> const &hopping_rates_in )
{
  hopping_rates_ = hopping_rates_in ;
}

/** Creates the transition matrices for the given time step
 *
 *  The transition matrix T(dt) is given by
 *  T(dt) = exp( K * dt )
 *
 *  The columns of this matrix correspond to T_ij = P( i <- j ), the probability of transitioning from i to j
 *  in a time dt.
 */
void MolecularState::constructTransitionMatrices()
{
  for (int i = 0 ; i < num_groups_ ; i++)
  {
    transition_matrices_[i] = (time_step_ * hopping_rates_[i]).exp() ;
  }
}

/** Sets the time step for the hopping dynamics and constructs the corresponding transition matrix
 *
 * @param [in] time_step_in the time step for the hopping dynamics.
 */
void MolecularState::setTimeStep( double const time_step_in )
{
  time_step_ = time_step_in ;
  constructTransitionMatrices() ;
  constructCumulativeProbs() ;
}

/** Samples the initial states for the system uniformly for each independent group of states.
 *
 */
void MolecularState::sampleInitialStatesUniformly()
{
  // samples each of the states uniformly
  for ( int i = 0 ; i < num_groups_ ; i++ )
  {
    current_states_(i) = int_distributions_[i](generator_) ;
  }
}

/**
 *
 */
void MolecularState::attemptHops()
{
  int i,k, state_index ;
  for ( i = 0 ; i < num_groups_ ; i++ )
  {
    random_number_ = random_sampler_.sampleUniformVariable(0.0,1.0) ;
    k = 0 ;
//    while ( cumulative_probabilities_[i]( k , current_states_(i) ) < random_number_  )
//    {
//      k++ ;
//    }
//    current_states_(i) =  k ;
    state_index = 0 ;
    for (k = 0 ;  k < array_num_states_(i) ; k++)
    {
      state_index = k ;
      if (random_number_ < cumulative_probabilities_[i]( k , current_states_(i) ) )
      {
        break ;
      }
    }
    current_states_(i) = state_index ;
  }
}

/** Constructs cumulative probability matrices for hop checking
 *
 */
void MolecularState::constructCumulativeProbs()
{
  // loops over all groups of independent
  for (int i = 0 ; i < num_groups_ ; i++)
  {
    // loops over each row
    cumulative_probabilities_[i].row(0) = transition_matrices_[i].row(0) ;
    for (int j = 0 ; j < array_num_states_(i)-1 ; j++)
    {
      cumulative_probabilities_[i].row(j+1) = cumulative_probabilities_[i].row(j) + transition_matrices_[i].row(j+1);
    }
  }
}

/** Test for markov chain stuff
 *
 * @param num_steps
 */
void MolecularState::runMarkovChain( int const num_steps )
{
//  ArrayXXi markov_chain_states = ArrayXXi::Zero(num_groups_ , num_steps+1) ;
//  sampleInitialStatesUniformly() ;
  resetStateCounts() ;
  for (int k = 0 ; k < num_steps ; k++)
  {
//    markov_chain_states.col(k) = current_states_ ;
    updateStateCounts() ;
    attemptHops() ;
  }

}

/** Updates the counts for states over a given run.
 *
 */
void MolecularState::updateStateCounts()
{
  // loops over each group
  for (int i = 0 ; i < num_groups_ ; i++)
  {
    // adds to the counter for the current state
    state_counts_[i]( current_states_(i) )++ ;
  }
}

/** Resets all of the state counts to zero.
 *
 */
void MolecularState::resetStateCounts()
{
  // loops over each group
  for (int i = 0 ; i < num_groups_ ; i++)
  {
    // sets all the counts to zero
    state_counts_[i].fill(0) ;
  }
}

/** Updates isotropic hfccs based on a markov chain run over a number of steps
 *
 * @param [in] num_steps
 * @param [out] couplings
 */
void MolecularState::updateCouplings( int const num_steps, ArrayXd &couplings )
{
  // runs the Markov chain
  runMarkovChain( num_steps ) ;

  // updates the average couplings over the run for each group
  for (int i = 0 ; i < num_groups_ ;  i++)
  {
    updateAverageCouplings( i ) ;
    // updates the couplings in the input array of couplings
    for (int k = 0 ; k < array_num_sites_(i) ; k++)
    {
      couplings(site_indices_[i](k)) = average_couplings_[i](k) ;
    }
  }
}

/** Calculates the average probabilities of being in a given state for a given group over a run
 *
 * @param [in] i_group group to be updated
 * @return average probabilities of being in a state for the given group from counts over a markov chain run
 */
ArrayXd MolecularState::calculateStateProbabilities( const int i_group )
{
  return ( 1.0 / ( (double) state_counts_[i_group].sum() ) )
    * (state_counts_[i_group].cast<double>()) ;
}

/**
 *
 * @param [in] i_group index of the independent group of markov states for which
 */
void MolecularState::updateAverageCouplings( int const i_group )
{
  state_probabilities_[i_group] = calculateStateProbabilities( i_group ) ;

  average_couplings_[i_group].fill(0.0) ;

  for (int k = 0 ; k < array_num_states_(i_group) ; k++)
  {
    average_couplings_[i_group] = average_couplings_[i_group] + state_probabilities_[i_group](k) * site_couplings_[i_group].col(k)  ;
  }
}

MolecularState::MolecularState( const MarkovStateParameters &params )
{
  // extract the number of independent groups
  num_groups_ = (int) params.getSiteIndices().size() ;
  // set the sizes of the arrays of the number of sites and states in each group
  array_num_sites_.resize(num_groups_) ;
  array_num_states_.resize(num_groups_) ;

  // resize the sizes of the site indices and site couplings arrays
  site_indices_.resize(num_groups_) ;
  site_couplings_.resize(num_groups_) ;
  hopping_rates_.resize(num_groups_) ;
  transition_matrices_.resize(num_groups_) ;
  cumulative_probabilities_.resize(num_groups_) ;
  state_counts_.resize(num_groups_);
  state_probabilities_.resize(num_groups_) ;
  average_couplings_.resize(num_groups_) ;

  // copy the input stuff to the member functions
  site_indices_ = params.getSiteIndices() ;
  site_couplings_ = params.getSiteCouplings() ;

  for (int i = 0 ; i < num_groups_ ; i++)
  {
    // extracts the number of sites for a given group
    array_num_sites_(i) = (int) site_couplings_[i].rows() ;

    // extracts the number of states for a given group of independent states
    array_num_states_(i) = (int) site_couplings_[i].cols() ;

    // fills the hopping rates matrix and transition matrix for a group of independent sites with zeros
    hopping_rates_[i] = ArrayXXd::Zero(array_num_states_(i), array_num_states_(i)) ;
    transition_matrices_[i] = ArrayXXd::Zero(array_num_states_(i), array_num_states_(i)) ;
    cumulative_probabilities_[i] = ArrayXXd::Zero(array_num_states_(i), array_num_states_(i)) ;
    state_counts_[i] = ArrayXi::Zero(array_num_states_(i)) ;
    state_probabilities_[i] = ArrayXd::Zero(array_num_states_(i)) ;
    average_couplings_[i] = ArrayXd::Zero(array_num_sites_(i)) ;
  }

  time_step_ = 0.0 ;

  // setting up the random number stuff randomly seeds the generator
  random_device seed_generator ;
  generator_.seed(seed_generator()) ;

  // sets up the distributions
  int_distributions_.resize(num_groups_) ;
  for (int i = 0 ; i< num_groups_ ; i++)
  {
    int_distributions_[i] = uniform_int_distribution<int>(0,array_num_states_(i)-1) ;
  }

  // sets up the array of current states
  current_states_ = ArrayXi::Zero(num_groups_) ;

  // sets the hopping rates from the parameters object
  setHoppingRates( params.getHoppingRates() ) ;

  // sets the number of steps for hopping from the parameters object
  num_steps_ = params.getNumHoppingSteps() ;

}

void MolecularState::setTimStepFromTotal( double const total_time )
{
  if (num_steps_!=0)
  {
    time_step_ = total_time / ( (double) num_steps_ ) ;
    constructTransitionMatrices() ;
    constructCumulativeProbs() ;
  }
  else
  {
    time_step_ = 0.0 ;
    constructTransitionMatrices() ;
    constructCumulativeProbs() ;
  }
}

MolecularState::MolecularState()
{
  num_groups_ = 0 ;
  num_steps_ = 0 ;
}

void MolecularState::resetSeed()
{
  random_sampler_.resetSeed() ;
}

void MolecularState::updateCouplings( ArrayXd &couplings )
{
  // runs the Markov chain
  runMarkovChain( num_steps_ ) ;

  // updates the average couplings over the run for each group
  for (int i = 0 ; i < num_groups_ ;  i++)
  {
    updateAverageCouplings( i ) ;
    // updates the couplings in the input array of couplings
    for (int k = 0 ; k < array_num_sites_(i) ; k++)
    {
      couplings(site_indices_[i](k)) = average_couplings_[i](k) ;
    }
  }
}

void MolecularState::printInfo()
{
  cout << "Num groups: " << num_groups_ << endl ;
  cout << "Num steps: " << num_steps_ << endl ;
  cout << "Num sites in each group: " << endl  << array_num_sites_.transpose() << endl ;
  cout << "Num states in each group: " << endl << array_num_states_.transpose() << endl ;
  cout << "Current states: " << current_states_.transpose() << endl ;
}




