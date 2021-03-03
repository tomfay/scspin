//
// Created by Tom on 04/12/2017.
//

#include "SpinStateEvolver.h"

/** SpinStateEvolver constructor - takes in all parameters other than the time step
 *
 * @param [in] nuc_couplings_1_in nuclear hyperfine couplings for radical 1
 * @param [in] nuc_couplings_2_in nuclear hyperfine couplings for radical 2
 * @param [in] exchange_coupling_in the exchange coupling - assumes H_12 = exchange_coupling * S_1.S_2
 * @param [in] singlet_rate_in k_S - singlet reaction rate
 * @param [in] triplet_rate_in k_T - triplet reaction rate
 * @param [in] omega_1_in vector defining zeeman frequency for radical 1
 * @param [in] omega_2_in vector defining zeeman frequency for radical 2
 *
 * creates an object that propagates the semi classical spin dynamics for a given set of spin parameters
 */
SpinStateEvolver::SpinStateEvolver( const ArrayXd nuc_couplings_1_in,
                                    const ArrayXd nuc_couplings_2_in,
                                    const double exchange_coupling_in,
                                    const double singlet_rate_in,
                                    const double triplet_rate_in,
                                    const Vector3d omega_1_in,
                                    const Vector3d omega_2_in )
{
  // the nuclear couplings are set here
  nuc_couplings_1_ = nuc_couplings_1_in ;
  nuc_couplings_2_ = nuc_couplings_2_in ;

  // the number of nuclear spins on each radicla is extracted
  num_nuc_spins_1_ = (int) nuc_couplings_1_.size() ;
  num_nuc_spins_2_ = (int) nuc_couplings_2_.size() ;

  // the zeeman frequency vectors are set here
  omega_1_ = omega_1_in ;
  omega_2_ = omega_2_in ;

  // sets the exchange coupling
  exchange_coupling_ = exchange_coupling_in ;

  // sets the rate constants and assigns the values of k_bar and delta_k
  singlet_rate_ = singlet_rate_in ;
  triplet_rate_ = triplet_rate_in ;
  k_bar = 0.25 * singlet_rate_ + 0.75 * triplet_rate_ ;
  delta_k = 0.25 * (singlet_rate_ - triplet_rate_ ) ;
  dephasing_rate_ = 0.0 ;

  // sets anisotropic to false
  is_aniso_ = false ;

  return ;
}

/** SpinStateEvolver constructor - takes in all parameters other than the time step
 *
 * @param [in] nuc_couplings_1_in nuclear hyperfine couplings for radical 1
 * @param [in] nuc_couplings_2_in nuclear hyperfine couplings for radical 2
 * @param [in] exchange_coupling_in the exchange coupling - assumes H_12 = exchange_coupling * S_1.S_2
 * @param [in] singlet_rate_in k_S - singlet reaction rate
 * @param [in] triplet_rate_in k_T - triplet reaction rate
 * @param [in] omega_1_in vector defining zeeman frequency for radical 1
 * @param [in] omega_2_in vector defining zeeman frequency for radical 2
 *
 * creates an object that propagates the semi classical spin dynamics for a given set of spin parameters
 */
SpinStateEvolver::SpinStateEvolver( const ArrayXXd aniso_nuc_couplings_1_in,
                                    const ArrayXXd aniso_nuc_couplings_2_in,
                                    const Array33d coupling_tensor_in,
                                    const double singlet_rate_in,
                                    const double triplet_rate_in,
                                    const Vector3d omega_1_in,
                                    const Vector3d omega_2_in )
{
  // the nuclear couplings are set here
  aniso_nuc_couplings_1_ = aniso_nuc_couplings_1_in ;
  aniso_nuc_couplings_2_ = aniso_nuc_couplings_2_in ;

  // the number of nuclear spins on each radicla is extracted
  num_nuc_spins_1_ = (int) nuc_couplings_1_.size() ;
  num_nuc_spins_2_ = (int) nuc_couplings_2_.size() ;

  // the zeeman frequency vectors are set here
  omega_1_ = omega_1_in ;
  omega_2_ = omega_2_in ;

  // sets the exchange coupling
  aniso_e_coupling_tensor_ = coupling_tensor_in ;

  // sets the rate constants and assigns the values of k_bar and delta_k
  singlet_rate_ = singlet_rate_in ;
  triplet_rate_ = triplet_rate_in ;
  k_bar = 0.25 * singlet_rate_ + 0.75 * triplet_rate_ ;
  delta_k = 0.25 * (singlet_rate_ - triplet_rate_ ) ;
  dephasing_rate_ = 0.0 ;

  // sets to anisotropic
  is_aniso_ = true ;

  return ;
}

/** Constructs a recombination matrix such that d/dt(v) = ... + K * v
 *
 * @return the matrix that acts on the electron spin variables to produce recombination
 */
MatrixXd SpinStateEvolver::constructRecombinationMatrix()
{
  MatrixXd recombination_matrix = MatrixXd::Zero(16,16) ;

  // diagonal element are all k_bar
  recombination_matrix.diagonal() = -k_bar * VectorXd::Ones(16) ; // set the diagonal

  // The spin vector component parts off diagonal blocks
  // d/dt S_1i = ... - k_bar S_1i + delta_k S_2i
  recombination_matrix.block(3,0,3,3).diagonal() = delta_k * VectorXd::Ones(3) ; // set off diagonal
  recombination_matrix.block(0,3,3,3).diagonal() = delta_k * VectorXd::Ones(3) ; // set off diagonal

  int n, m, i, j ;
  // the off diagonal components of the diagonal part of the  tensor and unit operator
  // d/dt T_ii = ... -k_bar T_ii - SUM_(j=/=i) delta_k T_jj + delta_k (unit /4 )
  // d/dt (unit / 4) = -k_bar (unit/4) + delta_k SUM_j T_jj
  for (i = 0 ; i < 3 ; i++)
  {
      // The T_ii to unit/4 couplings
      n = (3*i) + i + 6 ; // T_ii index
      recombination_matrix(n,15) = delta_k ;
      recombination_matrix(15,n) = delta_k ;
      // The T_ii to T_jj couplings
      for (j = i+1 ; j < 3 ; j++)
      {
        m = (3*j) + j + 6 ;
        recombination_matrix(n,m) = -delta_k ;
        recombination_matrix(m,n) = -delta_k ;
      }
  }

  // the off diagonal components of the tensor
  // d/dt T_ij = - k_bar T_ij + delta_k T_ji
  for (i = 0 ; i < 3 ; i++)
  {
    for ( j = (i+1) ; j < 3 ; j++ )
    {
      n = (3*i) + j + 6 ; // T_ij index
      m = (3*j) + i + 6 ; // T_ji index
      recombination_matrix(n,m) = delta_k ;
      recombination_matrix(m,n) = delta_k ;
    }
  }
//  cout << recombination_matrix << endl;
  return recombination_matrix ;
}

/** Constructs a coupling matrix for the electron spin variables such that d/dt (v) = ... + C * v
 *
 * @param [in] coupling_tensor the coupling tensors that contains dipolar and exchange terms such that
 *  H_12 = S_1 . coupling_tensor . S_2
 * @return the coupling matrix that acts on the electron spin variables
 */
MatrixXd SpinStateEvolver::constructCouplingMatrix( Matrix3d const coupling_tensor )
{
  // See Lachlan's "derivation.pdf" document for details
  MatrixXd coupling_matrix = MatrixXd::Zero(16,16) ;
  int i,n ;
  // first the S_1 to T coupling
  for (i = 0 ; i < 3 ; i++)
  {
    // S_1x part
    n = tensorIndex(1 , i) ; // T_2i couplings
    coupling_matrix(n,0) = 0.5 * coupling_tensor(2,i) ;
    coupling_matrix(0,n) = - 0.5 * coupling_tensor(2,i) ;

    n = tensorIndex(2, i) ; // T_3i couplings
    coupling_matrix(n,0) = -0.5 * coupling_tensor(1,i) ;
    coupling_matrix(0,n) = 0.5 * coupling_tensor(1,i) ;

    // S_1y part
    n = tensorIndex(0,i) ; // T_1i
    coupling_matrix(n,1) = -0.5 * coupling_tensor(2,i) ;
    coupling_matrix(1,n) = 0.5 * coupling_tensor(2,i) ;

    n = tensorIndex(2,i) ; // T_3i
    coupling_matrix(n,1) = 0.5 * coupling_tensor(0,i) ;
    coupling_matrix(1,n) = -0.5 * coupling_tensor(0,i) ;

    // S_1z part
    n = tensorIndex(0,i) ; // T_1i
    coupling_matrix(n,2) = 0.5 * coupling_tensor(1,i) ;
    coupling_matrix(2,n) = -0.5 * coupling_tensor(1,i) ;

    n = tensorIndex(1,i) ; // T_2i
    coupling_matrix(n,2) = -0.5 * coupling_tensor(0,i) ;
    coupling_matrix(2,n) = 0.5 * coupling_tensor(0,i) ;
  }

  // then the S_2 to T coupling
  for (i = 0 ; i < 3 ; i++)
  {
    // S_2x part
    n = tensorIndex( i , 1 ) ; // T_i2 couplings
    coupling_matrix(n,3) = 0.5 * coupling_tensor(i,2) ;
    coupling_matrix(3,n) = - 0.5 * coupling_tensor(i,2) ;

    n = tensorIndex( i , 2 ) ; // T_i3 couplings
    coupling_matrix(n,3) = -0.5 * coupling_tensor(i,1) ;
    coupling_matrix(3,n) = 0.5 * coupling_tensor(i,1) ;

    // S_2y part
    n = tensorIndex( i , 0 ) ; // T_i1 couplings
    coupling_matrix(n,4) = -0.5 * coupling_tensor(i,2) ;
    coupling_matrix(4,n) = 0.5 * coupling_tensor(i,2) ;

    n = tensorIndex( i , 2 ) ; // T_i3 couplings
    coupling_matrix(n,4) = 0.5 * coupling_tensor(i,0) ;
    coupling_matrix(4,n) = -0.5 * coupling_tensor(i,0) ;

    // S_2z part
    n = tensorIndex( i , 0 ) ; // T_i1 couplings
    coupling_matrix(n,5) = 0.5 * coupling_tensor(i,1) ;
    coupling_matrix(5,n) = -0.5 * coupling_tensor(i,1) ;

    n = tensorIndex( i , 1 ) ; // T_i2 couplings
    coupling_matrix(n,5) = -0.5 * coupling_tensor(i,0) ;
    coupling_matrix(5,n) = 0.5 * coupling_tensor(i,0) ;

  }

  return coupling_matrix ;
}

/** Outputs a vector rotated by a given rotation vector
 *
 * @param [in] rotation_vector vector about which the rotation is performs, direction defines axis and
 *  size defines rotation angle
 * @param [in] vector_in the vector to be rotated
 * @return the rotated vector
 *
 * Rotation is performed analytically using the rodrigues formula.
 */
Vector3d SpinStateEvolver::rotateVector( Vector3d const &rotation_vector , Vector3d const &vector_in )
{
  #if USE_CAYLEY
    double angle_sq = 0.25*(rotation_vector.squaredNorm()) ;
    return (vector_in + (2.0/(1.0+angle_sq))
      *(0.5*rotation_vector.cross(vector_in)
      + 0.25*rotation_vector *(rotation_vector.dot(vector_in))
      - angle_sq * vector_in) ) ;
  #else
    // rotation angle is the norm of the rotation vector
    double angle = rotation_vector.norm();
    // if angle < some tolerance no rotation is performed
    if (abs(angle)<TOLERANCE)
    {
      return vector_in ;
    }
    // generates the normalised rotation vector about which the rotation is performed
    Vector3d normed_rotation_vector =  (1.0 / angle) * rotation_vector  ;
    // outputes the rotated vector using Rodrigues formula
    return (cos(angle) * vector_in
      + (sin(angle)) * normed_rotation_vector.cross(vector_in)
     + ( 1.0 - cos(angle) ) * ( normed_rotation_vector.dot( vector_in ) ) * normed_rotation_vector) ;
  #endif
}

/** Returns the index of the electron spin variable corresponding to tensor component T_ij
 *
 * @param [in] i index for spin 1 component
 * @param [in] j index for spin 2 component
 * @return index in the electron spin variable vector corresponding the the T_ij element
 */
int SpinStateEvolver::tensorIndex( int const i, int const j )
{
  return (3 * i) + j + 6 ;
}

/** Sets the time step for the propagator
 *
 * @param [in] time_step_in the time step for the simulation
 */
void SpinStateEvolver::setTimeStep( const double time_step_in )
{
  time_step_ = time_step_in ;
}

/** Constructs the electron spin coupling tensor - currently only uses echange coupling
 *
 * @return the 3x3 electron spin coupling tensor
 *
 * Returns the electron spin coupling tensor assuming only exchange coupling is present.
 * This is defined such that H_12 = S_1 . coupling_tensor . S_2 = exchange coupling S_1 . S_2
 */
Matrix3d SpinStateEvolver::constructCouplingTensor()
{
//  Matrix3d coupling_tensor ;
//  coupling_tensor <<   -0.9999847 , -0.7369246 , 0.5112104,
//                                0.7369246 ,-0.0826997 , 0.0655341,
//                                0.5112104 , 0.0655341 , -0.5620818 ;
//  return coupling_tensor ;
  if (is_aniso_)
  {
    return aniso_e_coupling_tensor_ ;
  }
  else
  {
    return exchange_coupling_ * Matrix3d::Identity();
  }
}

/** Constructs the propagator for the electron spin variables corresponding to the recombination process
 *
 * @return recombination propagator for dt/2 this is exp( dt/2 K )
 */
MatrixXd SpinStateEvolver::constructRecombEvolutionOperator()
{
  MatrixXd recomb_matrix = constructRecombinationMatrix() ;
  MatrixXd dephasing_matrix = constructDephasingMatrix() ;
  MatrixXd relaxation_matrix = constructRelaxationMatrix() ;
  MatrixXd matrix = (recomb_matrix) + (dephasing_matrix) + (relaxation_matrix);

  #if USE_CAYLEY
    return ( ( (MatrixXd::Identity(16,16)) - matrix * (0.25*time_step_)).inverse()
    * ((MatrixXd::Identity(16,16)) + matrix * (0.25*time_step_)) ) ;
  #else
    return (0.5 * time_step_ * (matrix) ).exp() ;
  #endif
}

/** Constructs the propagator for the electron spin coupling terms
 *
 * @return propagator for dt/2 for electron spin coupling exp(dt/2 C)
 *
 * Currently assumes only exchange coupling is present.
 */
MatrixXd SpinStateEvolver::constructCouplingEvolutionOperator()
{
  MatrixXd matrix = constructCouplingMatrix( constructCouplingTensor() ) ;
  #if USE_CAYLEY
    return ( (MatrixXd::Identity(16,16) - matrix * 0.25*time_step_).inverse()
      * ((MatrixXd::Identity(16,16)) + matrix * 0.25*time_step_) ) ;
  #else
    return (0.5 * time_step_ *  matrix).exp() ;
  #endif

}

/** Creates all the propagators for the electron spin variables. This includes the coupling and recombination components.
 * Also creates the matrices for split operator propagation: exp(dt/2 K) * exp( dt/2 C) and exp(dt/2 C) * exp( dt/2 K) .
 *
 */
void SpinStateEvolver::createEvolutionOperators()
{
  recomb_evolution_operator_ = constructRecombEvolutionOperator() ;
  coupling_evolution_operator_ = constructCouplingEvolutionOperator() ;
  evol_mat_step_2_ = recomb_evolution_operator_ * coupling_evolution_operator_ ;
  evol_mat_step_1_ = coupling_evolution_operator_ * recomb_evolution_operator_  ;
}

/** Evolves an input state for a full time step dt using a split operator/analytic rotation method.
 *
 * @param [in,out] state The state which is evolves from t = t_0 to t = t_0 + dt .
 */
void SpinStateEvolver::evolveState( SpinState &state )
{
  // Apply the recombination evolution operator for dt/2
//  state.electron_spin_variables_ = recomb_evolution_operator_ * state.electron_spin_variables_ ;
  // Apply the coupling evolution operator for dt/2
//  state.electron_spin_variables_ = coupling_evolution_operator_ * state.electron_spin_variables_ ;

  // Applies the combined recombination and coupling propgators for dt/2
  state.electron_spin_variables_ = evol_mat_step_1_ * state.electron_spin_variables_ ;
  // Evolve the electron spins in the hyperfine field of the nuclei at t for dt/2
  evolveElectronSpins( state );

  // Evolve the nuclear spins around the electron spins at t+dt/2 for dt
  evolveNuclearSpins( state );

  // Evolve the electron spins the the hyperfine field of the nuclear at t+dt for dt/2
  evolveElectronSpins( state );

  // Apply the coupling evolution operator for dt/2
//  state.electron_spin_variables_ = coupling_evolution_operator_ * state.electron_spin_variables_ ;
  // Apply the recombination evolution operator for dt/2
//  state.electron_spin_variables_ = recomb_evolution_operator_ * state.electron_spin_variables_ ;

  // Applies the combined coupling and recombination propgators for dt/2
  state.electron_spin_variables_ = evol_mat_step_2_ * state.electron_spin_variables_ ;
}

/** Updates the instantaneous effective zeeman frequencies for the electron spins.
 *
 * @param state state for which the instantaneous electron zeeman terms are calculated
 *
 * Updates the instantaneous effective zeeman frequency vectors for the electron spins. This is the sum of the external
 * zeeman term and the instantaneous hyperfine field.
 */
void SpinStateEvolver::updateOmegaEffs( SpinState &state )
{
  // static zeeman terms
  omega_eff_1_ = omega_1_ ;
  omega_eff_2_ = omega_2_ ;
  int i ;
  // adds the instantaneous hyperfine fields to both effective zeeman frequencies
  for (i = 0 ; i < num_nuc_spins_1_ ; i++)
  {
    omega_eff_1_ = omega_eff_1_ + nuc_couplings_1_(i) * state.nuc_spins_1_.col(i) ;
  }

  for (i = 0 ; i < num_nuc_spins_2_ ; i++)
  {
    omega_eff_2_ = omega_eff_2_ + nuc_couplings_2_(i) * state.nuc_spins_2_.col(i) ;
  }

}

/** Evolves the electron spin variables in the zeeman and instantaneous hyperfine fields for dt/2.
 *
 * @param state state for which the electron spin variables are evolved.
 */
void SpinStateEvolver::evolveElectronSpins( SpinState &state )
{
  // Update the omega_effs with the current hyperfine fields
  updateOmegaEffs( state ) ;

  // create the rotation vector for spin 1
  rotation_vector_ =   ( 0.5 * time_step_ ) * omega_eff_1_  ;

  // rotate the S_1 vector
  new_vector_ = rotateVector( rotation_vector_ , state.electron_spin_variables_.segment(0,3) ) ;
  state.electron_spin_variables_.segment(0,3) = new_vector_ ;

  // Rotate the T components
  for (int i = 0 ; i <3 ; i++)
  {
    new_vector_ = rotateVector( rotation_vector_ , state.getTensorCol(i) ) ;
    state.setTensorCol( i , new_vector_ );
  }

  // create the rotation vector for spin 1
  rotation_vector_ =   ( 0.5 * time_step_ ) * omega_eff_2_  ;

  // rotate the S_2 vector
  new_vector_ = rotateVector( rotation_vector_ , state.electron_spin_variables_.segment(3,3) ) ;
  state.electron_spin_variables_.segment(3,3) = new_vector_ ;

  // Rotate the T components
  for (int i = 0 ; i <3 ; i++)
  {
    new_vector_ = rotateVector( rotation_vector_ , state.getTensorRow(i) ) ;
    state.setTensorRow( i , new_vector_ );
  }
}

/** Evolves the nuclear spins around the instantaneous electron spin vectors for dt.
 *
 * @param state state for which the nuclear spins are evolved
 */
void SpinStateEvolver::evolveNuclearSpins( SpinState &state )
{
  int k ;
  for (k = 0 ;  k < num_nuc_spins_1_ ; k++)
  {
    // calculates the rotation vector for nuc spin 1,k
    rotation_vector_ = ( nuc_couplings_1_(k) * time_step_ ) * state.getNormedSpinVector(1) ;
    state.nuc_spins_1_.col(k) = rotateVector( rotation_vector_ , state.nuc_spins_1_.col(k) );
  }
  for (k = 0 ;  k < num_nuc_spins_2_ ; k++)
  {
    // calculates the rotation vector for nuc spin 2,k
    rotation_vector_ = ( nuc_couplings_2_(k) * time_step_ ) * state.getNormedSpinVector(2) ;
    state.nuc_spins_2_.col(k) = rotateVector( rotation_vector_ , state.nuc_spins_2_.col(k) );

  }
}

/** Default constructor for SpinStateEvolver
 *
 */
SpinStateEvolver::SpinStateEvolver()
{

}

/** Constructs an evolver object with parameters set from a parameters object
 *
 * @param [in] parameters_in a Parameters object that contains parameters for the spin dynamics
 */
SpinStateEvolver::SpinStateEvolver( Parameters const &parameters_in )
{
  if( parameters_in.containsData() )
  {
    // the nuclear couplings are set here
    nuc_couplings_1_ = parameters_in.getNucCouplings(0) ;
    nuc_couplings_2_ = parameters_in.getNucCouplings(1) ;

    // the number of nuclear spins on each radicla is extracted
    num_nuc_spins_1_ = (int) nuc_couplings_1_.size() ;
    num_nuc_spins_2_ = (int) nuc_couplings_2_.size() ;

    // the zeeman frequency vectors are set here
    omega_1_ = parameters_in.getOmega(0) ;
    omega_2_ = parameters_in.getOmega(1) ;

    // sets the exchange coupling
    exchange_coupling_ = parameters_in.getExchangeCoupling() ;

    // sets the rate constants and assigns the values of k_bar and delta_k
    singlet_rate_ = parameters_in.getSingletRate() ;
    triplet_rate_ = parameters_in.getTripletRate() ;
    k_bar = 0.25 * singlet_rate_ + 0.75 * triplet_rate_ ;
    delta_k = 0.25 * (singlet_rate_ - triplet_rate_ ) ;
    dephasing_rate_ = parameters_in.getDephasingRate() ;
    gamma_ = parameters_in.getGammaRelaxationRates() ;

    // sets anisotropic coupling parameters if set
    is_aniso_ = false ;
    if (parameters_in.anisoParamsSet())
    {
      aniso_nuc_couplings_1_ = parameters_in.getAnisoNucCouplings(0) ;
      aniso_nuc_couplings_2_ = parameters_in.getAnisoNucCouplings(1) ;
      aniso_e_coupling_tensor_ = (Matrix3d) parameters_in.getECouplingTensor() ;
//      cout << aniso_nuc_couplings_1_ << endl ;
//      cout << aniso_nuc_couplings_2_ << endl ;
      is_aniso_ = true ;
    }

  }
}

void SpinStateEvolver::printInfo()
{
  cout << "Omega 1:" << endl ;
  cout << omega_1_.transpose() << endl ;
  cout << "Omega 2:" << endl ;
  cout << omega_2_.transpose() << endl ;
  cout << "Singlet rate: " << singlet_rate_ << endl ;
  cout << "Triplet rate: " << triplet_rate_ << endl ;
  cout << "Exchange couplings: " << exchange_coupling_ << endl ;
  cout << "Num nuc 1: " << num_nuc_spins_1_ << endl ;
  cout << "Num nuc 2: " << num_nuc_spins_2_ << endl ;
  cout << "Nuc couplings 1:" << endl ;
  cout << nuc_couplings_1_ << endl ;
  cout << "Nuc couplings 2:" << endl ;
  cout << nuc_couplings_2_ << endl ;
}

const ArrayXd &SpinStateEvolver::getNucCouplings1() const
{
  return nuc_couplings_1_;
}

const ArrayXd &SpinStateEvolver::getNucCouplings2() const
{
  return nuc_couplings_2_;
}

void SpinStateEvolver::setNucCouplings1( const ArrayXd &nuc_couplings_1_ )
{
  SpinStateEvolver::nuc_couplings_1_ = nuc_couplings_1_;
}

void SpinStateEvolver::setNucCouplings2( const ArrayXd &nuc_couplings_2_ )
{
  SpinStateEvolver::nuc_couplings_2_ = nuc_couplings_2_;
}

/** Constructs a dephasing matrix, the action of a singlet-triplfet dephasing rate on the electron spin variables such
 *  that d/dt (v) = ... + D * v
 *
 * @return a matrix that acts on the electron spin variables to cause S-T dephasing
 */
MatrixXd SpinStateEvolver::constructDephasingMatrix()
{
  MatrixXd dephasing_matrix = MatrixXd::Zero(16,16) ;

  // d/dt S_1i = ... - (k_STD/2) ( S_1i - S_2i)
  // The diagonal part for the spin vectors
  dephasing_matrix.block(0,0,6,6).diagonal() = -(dephasing_rate_*0.5) * VectorXd::Ones(6) ;

  // The spin vector component parts off diagonal blocks
  dephasing_matrix.block(3,0,3,3).diagonal() = (dephasing_rate_*0.5) * VectorXd::Ones(3) ; // set off diagonal
  dephasing_matrix.block(0,3,3,3).diagonal() = (dephasing_rate_*0.5) * VectorXd::Ones(3) ; // set off diagonal

  int n, m, i, j ;

  // the off diagonal components of the tensor i=/=j
  // d/dt T_ij = - (k_STD/2)( T_ij - T_ji)
  for (i = 0 ; i < 3 ; i++)
  {
    for ( j = (i+1) ; j < 3 ; j++ )
    {
      n = (3*i) + j + 6 ; // T_ij index
      m = (3*j) + i + 6 ; // T_ji index
      dephasing_matrix(n,m) = 0.5 * dephasing_rate_ ;
      dephasing_matrix(m,n) = 0.5 * dephasing_rate_ ;
      dephasing_matrix(n,n) = -0.5 * dephasing_rate_ ;
      dephasing_matrix(m,m) = -0.5 * dephasing_rate_ ;
    }
  }
//  cout << recombination_matrix << endl;
  return dephasing_matrix;
}

/** Evolves the electron spin variables in the zeeman and instantaneous hyperfine fields for dt/2, with anisotropic couplings.
 *
 * @param state state for which the electron spin variables are evolved.
 */
void SpinStateEvolver::evolveElectronSpinsAniso( SpinState &state )
{
  // Update the omega_effs with the current hyperfine fields, for anisotropic couplings
  updateOmegaEffsAniso( state ) ;

  // create the rotation vector for spin 1
  rotation_vector_ =   ( 0.5 * time_step_ ) * omega_eff_1_  ;

  // rotate the S_1 vector
  new_vector_ = rotateVector( rotation_vector_ , state.electron_spin_variables_.segment(0,3) ) ;
  state.electron_spin_variables_.segment(0,3) = new_vector_ ;

  // Rotate the T components
  for (int i = 0 ; i <3 ; i++)
  {
    new_vector_ = rotateVector( rotation_vector_ , state.getTensorCol(i) ) ;
    state.setTensorCol( i , new_vector_ );
  }

  // create the rotation vector for spin 1
  rotation_vector_ =   ( 0.5 * time_step_ ) * omega_eff_2_  ;

  // rotate the S_2 vector
  new_vector_ = rotateVector( rotation_vector_ , state.electron_spin_variables_.segment(3,3) ) ;
  state.electron_spin_variables_.segment(3,3) = new_vector_ ;

  // Rotate the T components
  for (int i = 0 ; i <3 ; i++)
  {
    // performing a silly check by adding a - sign here
    new_vector_ = rotateVector( rotation_vector_ , state.getTensorRow(i) ) ;
    state.setTensorRow( i , new_vector_ );
  }
}

/** Evolves nuclear spins with anisotropic couplings for a time step of dt
 *
 * @param state spin state for which the nuclear spins will be evolved
 */
void SpinStateEvolver::evolveNuclearSpinsAniso( SpinState &state )
{
  int k , j ;
  Vector3d electron_spin_vector ;
  electron_spin_vector = state.getNormedSpinVector(1) ;
  for (k = 0 ;  k < num_nuc_spins_1_ ; k++)
  {
    // calculates the rotation vector for nuc spin 1,k, theta_j(k) = a_jl(k)S_l
    for (j = 0 ; j <3 ; j++)
    {
      rotation_vector_(j) = ( time_step_ ) * (aniso_nuc_couplings_1_(3*j,k) * electron_spin_vector(0)
        + aniso_nuc_couplings_1_(3*j+1,k) * electron_spin_vector(1)
        + aniso_nuc_couplings_1_(3*j+2,k) * electron_spin_vector(2)) ;
    }

    state.nuc_spins_1_.col(k) = rotateVector( rotation_vector_ , state.nuc_spins_1_.col(k) );
  }

  electron_spin_vector = state.getNormedSpinVector(2) ;
  for (k = 0 ;  k < num_nuc_spins_2_ ; k++)
  {
    // calculates the rotation vector for nuc spin 2,k
    for (j = 0 ; j <3 ; j++)
    {
      rotation_vector_(j) = ( time_step_ ) * (aniso_nuc_couplings_2_(3*j,k) * electron_spin_vector(0)
        + aniso_nuc_couplings_2_(3*j+1,k) * electron_spin_vector(1)
        + aniso_nuc_couplings_2_(3*j+2,k) * electron_spin_vector(2)) ;
    }
    state.nuc_spins_2_.col(k) = rotateVector( rotation_vector_ , state.nuc_spins_2_.col(k) );

  }
}

void SpinStateEvolver::updateOmegaEffsAniso( SpinState &state )
{
  // static zeeman terms
  omega_eff_1_ = omega_1_ ;
  omega_eff_2_ = omega_2_ ;
  int i , j, k;
  // adds the instantaneous hyperfine fields to both effective zeeman frequencies
  for (i = 0 ; i < num_nuc_spins_1_ ; i++)
  {
    for (j = 0 ; j < 3 ; j++)
    {
      for (k = 0 ; k < 3 ; k++)
      {
        // omega_j = omega_j(ext) + sum_i sum_k + a_jk(i) * I_k(i)
        omega_eff_1_(j) = omega_eff_1_(j) + aniso_nuc_couplings_1_(3*j+k,i) * state.nuc_spins_1_(k,i) ;
      }
    }
  }
  // repeats for second spin
  for (i = 0 ; i < num_nuc_spins_2_ ; i++)
  {
    for (j = 0 ; j < 3 ; j++)
    {
      for (k = 0 ; k < 3 ; k++)
      {
        omega_eff_2_(j) = omega_eff_2_(j) + aniso_nuc_couplings_2_(3*j+k,i) * state.nuc_spins_2_(k,i) ;
      }
    }
  }

}

void SpinStateEvolver::evolveStateAniso( SpinState &state )
{
  // Apply the recombination evolution operator for dt/2
//  state.electron_spin_variables_ = recomb_evolution_operator_ * state.electron_spin_variables_ ;
  // Apply the coupling evolution operator for dt/2
//  state.electron_spin_variables_ = coupling_evolution_operator_ * state.electron_spin_variables_ ;

  // Applies the combined recombination and coupling propgators for dt/2
  state.electron_spin_variables_ = evol_mat_step_1_ * state.electron_spin_variables_ ;
  // Evolve the electron spins in the hyperfine field of the nuclei at t for dt/2
  evolveElectronSpinsAniso( state );

  // Evolve the nuclear spins around the electron spins at t+dt/2 for dt
  evolveNuclearSpinsAniso( state );

  // Evolve the electron spins the the hyperfine field of the nuclear at t+dt for dt/2
  evolveElectronSpinsAniso( state );

  // Apply the coupling evolution operator for dt/2
//  state.electron_spin_variables_ = coupling_evolution_operator_ * state.electron_spin_variables_ ;
  // Apply the recombination evolution operator for dt/2
//  state.electron_spin_variables_ = recomb_evolution_operator_ * state.electron_spin_variables_ ;

  // Applies the combined coupling and recombination propgators for dt/2
  state.electron_spin_variables_ = evol_mat_step_2_ * state.electron_spin_variables_ ;
}

double SpinStateEvolver::getSingletRate() const
{
  return singlet_rate_;
}

double SpinStateEvolver::getTripletRate() const
{
  return triplet_rate_;
}

double SpinStateEvolver::getTimeStep() const
{
  return time_step_;
}

void SpinStateEvolver::evolveStateAnisoSw( SpinState &state )
{
  // Applies the combined recombination and coupling propgators for dt/2
  state.electron_spin_variables_ = evol_mat_step_1_ * state.electron_spin_variables_ ;
  // Evolve the electron spins in the hyperfine field of the nuclei at t for dt/2
  evolveElectronSpinsAnisoSw( state );

  // Evolve the electron spins the the hyperfine field of the nuclear at t+dt for dt/2
  evolveElectronSpinsAnisoSw( state );

  // Applies the combined coupling and recombination propgators for dt/2
  state.electron_spin_variables_ = evol_mat_step_2_ * state.electron_spin_variables_ ;
}

/** Evolves the electron spin variables in the zeeman and instantaneous hyperfine fields for dt/2, with anisotropic couplings.
 *
 * @param state state for which the electron spin variables are evolved.
 */
void SpinStateEvolver::evolveElectronSpinsAnisoSw( SpinState &state )
{

  // create the rotation vector for spin 1
  rotation_vector_ =   ( 0.5 * time_step_ ) * omega_eff_1_  ;

  // rotate the S_1 vector
  new_vector_ = rotateVector( rotation_vector_ , state.electron_spin_variables_.segment(0,3) ) ;
  state.electron_spin_variables_.segment(0,3) = new_vector_ ;

  // Rotate the T components
  for (int i = 0 ; i <3 ; i++)
  {
    new_vector_ = rotateVector( rotation_vector_ , state.getTensorCol(i) ) ;
    state.setTensorCol( i , new_vector_ );
  }

  // create the rotation vector for spin 1
  rotation_vector_ =   ( 0.5 * time_step_ ) * omega_eff_2_  ;

  // rotate the S_2 vector
  new_vector_ = rotateVector( rotation_vector_ , state.electron_spin_variables_.segment(3,3) ) ;
  state.electron_spin_variables_.segment(3,3) = new_vector_ ;

  // Rotate the T components
  for (int i = 0 ; i <3 ; i++)
  {
    new_vector_ = rotateVector( rotation_vector_ , state.getTensorRow(i) ) ;
    state.setTensorRow( i , new_vector_ );
  }
}

void SpinStateEvolver::evolveElectronSpinVectors( SpinState &state )
{
  // Update the omega_effs with the current hyperfine fields
  updateOmegaEffsAniso( state ) ;

  // create the rotation vector for spin 1 for dt/4
  rotation_vector_ =   ( 0.25 * time_step_ ) * (omega_eff_1_  + 2.0*aniso_e_coupling_tensor_ * state.electron_spin_variables_.segment(3,3)) ;

  // rotate the S_1 vector
  new_vector_ = rotateVector( rotation_vector_ , state.electron_spin_variables_.segment(0,3) ) ;
  state.electron_spin_variables_.segment(0,3) = new_vector_ ;

  // create the rotation vector for spin 2 for dt/2
  rotation_vector_ =   ( 0.5 * time_step_ ) *( omega_eff_2_ + 2.0*aniso_e_coupling_tensor_ * state.electron_spin_variables_.segment(0,3)) ;

  // rotate the S_2 vector
  new_vector_ = rotateVector( rotation_vector_ , state.electron_spin_variables_.segment(3,3) ) ;
  state.electron_spin_variables_.segment(3,3) = new_vector_ ;

  // create the rotation vector for spin 1 for dt/4
  rotation_vector_ =   ( 0.25 * time_step_ ) * (omega_eff_1_ + 2.0*aniso_e_coupling_tensor_ * state.electron_spin_variables_.segment(3,3) ) ;

  // rotate the S_1 vector
  new_vector_ = rotateVector( rotation_vector_ , state.electron_spin_variables_.segment(0,3) ) ;
  state.electron_spin_variables_.segment(0,3) = new_vector_ ;

}

void SpinStateEvolver::evolveStateAnisoVectors( SpinState &state )
{

  // Evolve the electron spins in the hyperfine field of the nuclei at t for dt/2
  evolveElectronSpinVectors( state );

  // Evolve the nuclear spins around the electron spins at t+dt/2 for dt
  evolveNuclearSpinsAniso( state );

  // Evolve the electron spins the the hyperfine field of the nuclear at t+dt for dt/2
  evolveElectronSpinVectors( state );

}

void SpinStateEvolver::evolveNuclearSpinsAnisoAltNorm( SpinState &state )
{
  int k , j ;
  Vector3d electron_spin_vector ;
  double norm_const = (state.getUnitOperator()) ;
  norm_const = 1.0/norm_const ;
  electron_spin_vector = state.getSpinVector(1) ;
  electron_spin_vector = norm_const * electron_spin_vector ;
//  cout << "vector 1:" << electron_spin_vector.transpose() << "unit op:" << state.getUnitOperator() << endl ;
  for (k = 0 ;  k < num_nuc_spins_1_ ; k++)
  {
    // calculates the rotation vector for nuc spin 1,k
    for (j = 0 ; j <3 ; j++)
    {
      rotation_vector_(j) = ( time_step_ ) * (aniso_nuc_couplings_1_(3*j,k) * electron_spin_vector(0)
        + aniso_nuc_couplings_1_(3*j+1,k) * electron_spin_vector(1)
        + aniso_nuc_couplings_1_(3*j+2,k) * electron_spin_vector(2)) ;
    }

    state.nuc_spins_1_.col(k) = rotateVector( rotation_vector_ , state.nuc_spins_1_.col(k) );
  }

  electron_spin_vector = state.getSpinVector(2) ;
  electron_spin_vector = norm_const * electron_spin_vector ;
//  cout << "vector 2:" << electron_spin_vector.transpose() << "unit op:" << state.getUnitOperator() << endl ;
  for (k = 0 ;  k < num_nuc_spins_2_ ; k++)
  {
    // calculates the rotation vector for nuc spin 2,k
    for (j = 0 ; j <3 ; j++)
    {
      rotation_vector_(j) = ( time_step_) * (aniso_nuc_couplings_2_(3*j,k) * electron_spin_vector(0)
        + aniso_nuc_couplings_2_(3*j+1,k) * electron_spin_vector(1)
        + aniso_nuc_couplings_2_(3*j+2,k) * electron_spin_vector(2)) ;
    }
    state.nuc_spins_2_.col(k) = rotateVector( rotation_vector_ , state.nuc_spins_2_.col(k) );

  }
}

void SpinStateEvolver::evolveStateAnisoAltNorm( SpinState &state )
{

  // Applies the combined recombination and coupling propgators for dt/2
  state.electron_spin_variables_ = evol_mat_step_1_ * state.electron_spin_variables_ ;
  // Evolve the electron spins in the hyperfine field of the nuclei at t for dt/2
  evolveElectronSpinsAniso( state );

  // Evolve the nuclear spins around the electron spins at t+dt/2 for dt
  evolveNuclearSpinsAnisoAltNorm(state );

  // Evolve the electron spins the the hyperfine field of the nuclear at t+dt for dt/2
  evolveElectronSpinsAniso( state );

  // Applies the combined coupling and recombination propgators for dt/2
  state.electron_spin_variables_ = evol_mat_step_2_ * state.electron_spin_variables_ ;
}

/* Creates the relaxation matrix for a Lindblad operator of the form:
 *  sum_i,alpha gamma_i,alpha (S_ialpha rho S_ialpha  - 1/4 rho)
 *
 */
MatrixXd SpinStateEvolver::constructRelaxationMatrix()
{
  int n, m, i, j ;
  MatrixXd relaxation_matrix = MatrixXd::Zero(16,16) ;

  // first create the decay rates for the ix operators
  ArrayXd decay_rates = ArrayXd::Zero(6) ;
  // decay rates for the S_1alpha operators
  decay_rates(0) = 0.5*(gamma_(1)+gamma_(2));
  decay_rates(1) = 0.5*(gamma_(2)+gamma_(0));
  decay_rates(2) = 0.5*(gamma_(0)+gamma_(1));
  // decay rates for the S_2alpha operators
  decay_rates(3) = 0.5*(gamma_(4)+gamma_(5));
  decay_rates(4) = 0.5*(gamma_(5)+gamma_(3));
  decay_rates(5) = 0.5*(gamma_(3)+gamma_(4));

  // 0 to 5 correspond to the S_ialpha operators, so set the diagonal of this block to -decay_rates
  for ( i = 0 ; i < 6 ; i++ )
  {
    relaxation_matrix(i,i) = -decay_rates(i) ;
  }

  // The T_ij element relaxes at a rate -k_1i - k_2j so set the diagonal elements of this matrix to that
  for ( i = 0 ; i < 3 ; i++ )
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      n = (3*i) + j + 6 ; // T_ij index
      relaxation_matrix(n,n) = -(decay_rates(i)+decay_rates(j+3)) ;
    }
  }

  return relaxation_matrix;
}




