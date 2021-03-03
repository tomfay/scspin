//
// Created by Tom on 04/12/2017.
//

#include "SpinState.h"

/** Constructor for the SpinState class where the number of nuclear spins on both radicals are specified.
 *
 * @param [in] num_nuc_spins_1_in number of nuclear spins on radical 1
 * @param [in] num_nuc_spins_2_in number of nuclear spins on radical 2
 */
SpinState::SpinState( int const num_nuc_spins_1_in, int const num_nuc_spins_2_in )
{
  // read in the number of nuc spins
  num_nuc_spins_1_ = num_nuc_spins_1_in ;
  num_nuc_spins_2_ = num_nuc_spins_2_in ;

  // assumes the degeneracy is 2 - i.e. all spin 1/2
  nuc_spin_degeneracies_1_ = 2 * ArrayXi::Ones(num_nuc_spins_1_) ;
  nuc_spin_degeneracies_2_ = 2 * ArrayXi::Ones(num_nuc_spins_2_) ;

  // set vector lenghts
  nuc_spin_lengths_1_ = generateSpinLengths( nuc_spin_degeneracies_1_ ) ;
  nuc_spin_lengths_2_ = generateSpinLengths( nuc_spin_degeneracies_2_ ) ;

  // create electronic state variables - 16 of them in total
  electron_spin_variables_ = VectorXd::Zero( 16 ) ;

  // create nuc spin variables
  nuc_spins_1_ = MatrixXd::Zero( 3, num_nuc_spins_1_ ) ;
  nuc_spins_2_ = MatrixXd::Zero( 3, num_nuc_spins_2_ ) ;
}

/** Generates an array of spin vector lengths from an array of spin degeneracies.
 *
 * @param spin_degeneracies the quantum degeneracy of each spin state (2I_k +1) .
 * @return an array of spin vector lengths sqrt( I_k (I_k+1) ).
 */
ArrayXd SpinState::generateSpinLengths( ArrayXi const spin_degeneracies )
{
  ArrayXd spin_lengths_out = ArrayXd::Zero( spin_degeneracies.size() ) ;
  double x ;
  for (int i = 0 ; i < spin_lengths_out.size() ; i++)
  {
    x = (double) spin_degeneracies(i) ;
    x = 0.5*(x - 1.0 ) ;
    x = sqrt(x*(x+1.0)) ;
    spin_lengths_out(i) = x ;
  }
  return spin_lengths_out;
}

/** Default constructor for SpinState class.
 *
 */
SpinState::SpinState()
{

}

/** Returns the singlet probability for the state P_S = unit/4 - (T_11 + T_22 + T_33)
 *
 * @return singlet probability P_S for the state
 */
double SpinState::calculateSingletProbability()
{
  return electron_spin_variables_(15) // unit/4
    - ( electron_spin_variables_(6)   // T_11
      + electron_spin_variables_(10)  // T_22
      + electron_spin_variables_(14)  // T_33
    ) ;
}

/** Returns the singlet probability for the state P_T =  3 * unit/4 + (T_11 + T_22 + T_33)
 *
 * @return triplet probability P_T for the state
 */
double SpinState::calculateTripletProbability()
{
  return (3.0 * electron_spin_variables_(15))
    + ( electron_spin_variables_(6)
      + electron_spin_variables_(10)
      + electron_spin_variables_(14) );
}

/** Returns the total survival probability i.e. the unit operator.
 *
 * @return the unit operator i.e. survival probabilty for the state.
 */
double SpinState::calculateSurvivalProbability()
{
  return 4.0 * electron_spin_variables_(15) ;
}

/** Returns the index in the electron spin variable vector corresponding to T_ij.
 *
 * @param [in] i electron spin 1 component index
 * @param [out] j electron spin 2 component index
 * @return index corresponding to T_ij in the electron spin variable vector
 */
int SpinState::getTensorIndex( int const i, int const j )
{
  return (3*i) + j + 6;
}

/** Returns the a vector of tensor components correspond to a row of the tensor (T_i1, T_i2, T_i3) = S_1i (S_2x, S_2y, S_2z)
 *
 * @param [in] i the row index to be returned.
 * @return row of the electron spin tensor (T_i1, T_i2, T_i3) = S_1i (S_2x, S_2y, S_2z)
 */
Vector3d SpinState::getTensorRow( int const i )
{
  // index of first element is at ( 3*i +6 ) returns this segment of the vector
  return electron_spin_variables_.segment(getTensorIndex(i,0),3) ;
}

// returns the (T_1i, T_2i, T_3i)
/** Returns a column of the electron spin tensor (T_1i, T_2i, T_3i) = (S_1x, S_1y, S_1z) S_2i.
 *
 * @param [in] i column index to be returned.
 * @return column of the electron spin tensor (T_1i, T_2i, T_3i) = (S_1x, S_1y, S_1z) S_2i
 */
Vector3d SpinState::getTensorCol( int const i )
{
  // returns the 3-vector.
  return Vector3d( electron_spin_variables_(getTensorIndex(0,i)),
                   electron_spin_variables_(getTensorIndex(1,i)),
                   electron_spin_variables_(getTensorIndex(2,i))) ;
}

/**  Sets a row of the electron spin tensor to a new value (T_i1, T_i2, T_i3)
 *
 * @param [in] i the row index that is to be re-assigned.
 * @param new_row the new value of the row (T_i1, T_i2, T_i3)
 */
void SpinState::setTensorRow( int i, const Vector3d &new_row )
{
  electron_spin_variables_.segment(getTensorIndex(i,0),3) = new_row ;
}

/** Sets a electron spin tensor column (T_1i, T_2i, T_3i) to a new value.
 *
 * @param [in] i the index of the column to be re-assigned.
 * @param new_col the new values for the column (T_1i, T_2i, T_3i)
 */
void SpinState::setTensorCol( int const i, Vector3d const &new_col )
{
  electron_spin_variables_(getTensorIndex(0,i)) = new_col(0) ;
  electron_spin_variables_(getTensorIndex(1,i)) = new_col(1) ;
  electron_spin_variables_(getTensorIndex(2,i)) = new_col(2) ;
}

/** Returns the electron spin vector i normalized to its physical vector length sqrt(3/4)
 *
 * @param [in] i the electron spin index 1 or 2
 * @return the electron spin vector normalized to its semi-classical length
 */
Vector3d SpinState::getNormedSpinVector( int const i )
{
  Vector3d normed_spin_vec ;
  normed_spin_vec = electron_spin_variables_.segment( 3*(i-1) , 3);
  normed_spin_vec = (HALF_SQRT3 / normed_spin_vec.norm()) * normed_spin_vec ;
  return normed_spin_vec ;
//  return HALF_SQRT3 * (electron_spin_variables_.segment(3*(i-1) , 3).normalized()) ;
}

/** Returns the value of the electron spin vector i.
 *
 * @param [in] i the index of the spin vector to be obtained.
 * @return the value of this spin vector (S_ix, S_iy, S_iz)
 */
Vector3d SpinState::getSpinVector( int const i )
{
  return 2.0 * electron_spin_variables_.segment(3*(i-1),3) ;
}

/** Returns the value of the unit operator.
 *
 * @return the value of the electron spin unit operator.
 */
double SpinState::getUnitOperator()
{
  return 4.0 * electron_spin_variables_(15) ;
}

void SpinState::printInfo()
{
  cout << "Num nuc 1: " << num_nuc_spins_1_ << endl ;
  cout << "Num nuc 2: " << num_nuc_spins_2_ << endl ;
  cout << "Lengths 1 :" << endl ;
  cout << nuc_spin_lengths_1_.transpose() << endl ;
  cout << "Lengths 2 :" << endl ;
  cout << nuc_spin_lengths_2_.transpose() << endl ;
  cout << "Nuc spins 1 dims: " << nuc_spins_1_.rows() << " x " << nuc_spins_1_.cols() << endl ;
  cout << "Nuc spins 2 dims: " << nuc_spins_2_.rows() << " x " << nuc_spins_2_.cols() << endl ;
}

void SpinState::setNucSpinLengths( ArrayXd const nuc_spin_lengths_1, ArrayXd const nuc_spin_lengths_2 )
{
  double x ;
  nuc_spin_lengths_1_ = nuc_spin_lengths_1 ;
  nuc_spin_lengths_2_ = nuc_spin_lengths_2 ;

  for (int k = 0 ; k< num_nuc_spins_1_ ; k++)
  {
    x = nuc_spin_lengths_1_(k) ;
    x = x*x ;
    x = 0.5 + sqrt(x + 0.25) ;
    x = 2.0*x+1.0 ;
    nuc_spin_degeneracies_1_(k) = (int) x;
  }

  for (int k = 0 ; k< num_nuc_spins_2_ ; k++)
  {
    x = nuc_spin_lengths_2_(k) ; // x = sqrt(I_k(I_k+1))
    x = x*x ; // x = I_k(I_k+1)
    x = -0.5 + sqrt(x + 0.25) ; // x = I_k
    x = 2.0*x+1.0 ; // x = 2*I_k+1
    nuc_spin_degeneracies_2_(k) = (int) x;
  }

  return ;
}

double SpinState::calculateSingletProbabilityVectorOnly()
{
  return 0.25 - 4.0 * ( electron_spin_variables_(0)*electron_spin_variables_(3)
    + electron_spin_variables_(1)*electron_spin_variables_(4)
    + electron_spin_variables_(2)*electron_spin_variables_(5) ) ;
}

double SpinState::calculateTripletProbabilityVectorOnly()
{
  return 0.75 + 4.0 * ( electron_spin_variables_(0)*electron_spin_variables_(3)
    + electron_spin_variables_(1)*electron_spin_variables_(4)
    + electron_spin_variables_(2)*electron_spin_variables_(5) ) ;
}
