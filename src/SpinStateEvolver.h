//
// Created by Tom on 04/12/2017.
//

#ifndef SEMICLASSICAL_SPINSTATEEVOLVER_H
#define SEMICLASSICAL_SPINSTATEEVOLVER_H

// tolerance for minimum angle for rotation, angle=0 breaks method otherwise
#define TOLERANCE 1.0e-12
#define USE_CAYLEY true

#include <iostream>
#ifdef HEADERONLY
  #include "Eigen/Dense"
  #include "unsupported/Eigen/MatrixFunctions"

#else
  #include <Eigen/Dense>
  #include <unsupported/Eigen/MatrixFunctions>
#endif

#include <cmath>

#include "SpinState.h"
#include "Parameters.h"



using namespace std ;
using namespace Eigen ;

/** Class that contains methods for propagating a semi-classical radical pair spin state.
 *
 *  A split operator/analytic rotation propagator is used to propagate the electron spin state.
 */
class SpinStateEvolver
{
public:
  // Constructor that sets all evolution variables
  SpinStateEvolver( const ArrayXd nuc_couplings_1_in,
                    const ArrayXd nuc_couplings_2_in,
                    const double exchange_coupling_in,
                    const double singlet_rate_in,
                    const double triplet_rate_in,
                    const Vector3d omega_1_in,
                    const Vector3d omega_2_in );

  // Constructor that sets all evolution variables with anisotropic coupling constants
  SpinStateEvolver( const ArrayXXd aniso_nuc_couplings_1_in,
                    const ArrayXXd aniso_nuc_couplings_2_in,
                    const Array33d coupling_tensor_in,
                    const double singlet_rate_in,
                    const double triplet_rate_in,
                    const Vector3d omega_1_in,
                    const Vector3d omega_2_in ) ;

  // default constructor
  SpinStateEvolver();

  // constructor reading in parameters from a Parameters object
  SpinStateEvolver( Parameters const &parameters_in );

  // sets the time step
  void setTimeStep( const double time_step_in ) ;

  // creates the electron spin evolution operators
  void createEvolutionOperators() ;

  // evolves an input state for time step dt
  void evolveState( SpinState &state );

  // evolves an input state for time step dt with anisotropic couplings
  void evolveStateAniso( SpinState &state );

  // evolves an input state for time step dt with anisotropic couplings
  void evolveStateAnisoSw( SpinState &state );

  // evolves an input state for time step dt with anisotropic couplings using the vector only method
  void evolveStateAnisoVectors( SpinState &state );

  // evolves an input state for time step dt with anisotropic couplings using the unit inverse alternative normalisation
  void evolveStateAnisoAltNorm( SpinState &state );

  // prints info
  void printInfo() ;


  // getters and setters for nuclear couplings
  const ArrayXd &getNucCouplings1() const;
  const ArrayXd &getNucCouplings2() const;
  void setNucCouplings1( const ArrayXd &nuc_couplings_1_ );
  void setNucCouplings2( const ArrayXd &nuc_couplings_2_ );

  // update effective Zeeman term for each spin
  void updateOmegaEffsAniso( SpinState &state ) ;

  double getSingletRate() const;
  double getTripletRate() const;
  double getTimeStep() const;
private:
  // the number of nuclear spins on each radical
  int num_nuc_spins_1_, num_nuc_spins_2_ ;

  // boolean for whether anisotropic couplings are used
  bool is_aniso_ ;
  // arrays containing the nuclear hyperfine coupling constants - assumed isotropic for the moment
  ArrayXd nuc_couplings_1_, nuc_couplings_2_ ;
  // arrays containing the anisotropic coupling tensors for the hyperfine constants
  ArrayXXd aniso_nuc_couplings_1_ , aniso_nuc_couplings_2_ ;
  // electron coupling tensor
  Matrix3d aniso_e_coupling_tensor_ ;
  // exchange coupling, rate constants and time step for propagation
  double exchange_coupling_ , singlet_rate_, triplet_rate_ , delta_k , k_bar, time_step_ , dephasing_rate_;
  // relaxation superoperator parameters (gamma_1x gamma_1y gamma_1z gamma_2x gamma_2y gamma_2z)
  ArrayXd gamma_ ;

  // zeeman frequency terms such that H_1 = omega_i . S_i + ...
  Vector3d omega_1_ , omega_2_ ;
  // instantaneous zeeman frequencies for each radical and rotation vector/rotated vector objects
  Vector3d omega_eff_1_, omega_eff_2_ , rotation_vector_, new_vector_ ;
  // matrix propagators for electron spin variable evolution.
  MatrixXd recomb_evolution_operator_ , coupling_evolution_operator_ , evol_mat_step_1_, evol_mat_step_2_;

  // rotates an input vector
  Vector3d rotateVector( Vector3d const &rotation_vector , Vector3d const &vector_in ) ;

  // methods to construct electron spin coupling and recombination propagators
  MatrixXd constructRecombinationMatrix() ;
  MatrixXd constructCouplingMatrix( Matrix3d const coupling_tensor ) ;
  MatrixXd constructRecombEvolutionOperator() ;
  MatrixXd constructCouplingEvolutionOperator() ;
  Matrix3d constructCouplingTensor() ;
  MatrixXd constructDephasingMatrix() ;
  MatrixXd constructRelaxationMatrix() ;

  // returns index for electron spin tensor component in the electron spin variables vector
  int tensorIndex(int const i, int const j) ;

  // methods for evolving nuclear and electron spins - including updating instantanous zeeman frequencies for electrons
  void evolveElectronSpins( SpinState &state ) ;
  void evolveNuclearSpins( SpinState &state ) ;
  void updateOmegaEffs( SpinState &state ) ;
  void evolveElectronSpinsAnisoSw( SpinState &state );

  void evolveElectronSpinVectors( SpinState &state  );

  // methods for evolving nuclear and electron spins with anisotropic couplings - including updating instantanous zeeman frequencies for electrons
  void evolveElectronSpinsAniso( SpinState &state ) ;
  void evolveNuclearSpinsAniso( SpinState &state ) ;
  void evolveNuclearSpinsAnisoAltNorm( SpinState &state ) ;


};


#endif //SEMICLASSICAL_SPINSTATEEVOLVER_H
