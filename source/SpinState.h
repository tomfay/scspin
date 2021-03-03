//
// Created by Tom on 04/12/2017.
//

#ifndef SEMICLASSICAL_SPINSTATE_H
#define SEMICLASSICAL_SPINSTATE_H

#define HALF_SQRT3 0.86602540378

#include <iostream>
#ifdef HEADERONLY
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#else
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#endif
#include <cmath>

using namespace std;
using namespace Eigen;

/** Class that contains information about the semi-classical spin state of a molecule.
 *
 * Contains methods for extracting state information and changing values of spin tensors and spin vectors.
 *
 * Note that the electron spin variables are stored as a vector:
 * (S_1/2, S_2/2, T_11, T_12, T_13, T_21, T_22, T_23, T_31, T_32, T_33, unit/4)
 *
 * and the nuclear spin vectors are stored in two 3 x num_nuc_spins_i_  matrices.
 */
class SpinState
{
public:
  // constructor for spin state taking in the number of nuclei on each radical
  SpinState( int const num_nuc_spins_1_in, int const num_nuc_spins_2_in ) ;

  // default constructor for the spin state object
  SpinState();

  // methods for extracting state survival probabilities
  double calculateSingletProbability() ;
  double calculateTripletProbability() ;
  double calculateSurvivalProbability() ;
  double calculateSingletProbabilityVectorOnly() ;
  double calculateTripletProbabilityVectorOnly() ;

  // Matrices that contain the electron spin variables
  VectorXd electron_spin_variables_ ;
  MatrixXd nuc_spins_1_, nuc_spins_2_ ;

  // Methods for obtains electron spin tensor rows and columns and settign their values
  Vector3d getTensorRow( int const i ) ;
  Vector3d getTensorCol( int const i ) ;
  void setTensorRow(int const i, Vector3d const &new_col) ;
  void setTensorCol(int const i, Vector3d const &new_col) ;

  // Methods for extracting the spin vectors and unit operator
  Vector3d getNormedSpinVector( int const i ) ;
  Vector3d getSpinVector( int const i ) ;
  double getUnitOperator() ;

  // Method for generating the index of T_ij in the electron spin variables object
  int getTensorIndex(int const i, int const j) ;

  // Arrays holding the spin vector lengths for the nuclear spin
  ArrayXd nuc_spin_lengths_1_ ,nuc_spin_lengths_2_ ;

  // number of nuclear spins on each radical
  int num_nuc_spins_1_, num_nuc_spins_2_ ;

  // prints info
  void printInfo();

  // set the nuclear spin vector lengths
  void setNucSpinLengths( ArrayXd const nuc_spin_lengths_1 , ArrayXd const nuc_spin_lengths_2 );

private:



  // spin state degeneracies for each radical
  ArrayXi nuc_spin_degeneracies_1_, nuc_spin_degeneracies_2_ ;

  // method to generate spin vector lengths from
  ArrayXd generateSpinLengths( ArrayXi const spin_degeneracies ) ;
};


#endif //SEMICLASSICAL_SPINSTATE_H
