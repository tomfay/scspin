//
// Created by Tom on 05/12/2017.
//

#include "KahanSum.h"

/** Constuctor with setting dimensions of the array to be summed.
 *
 * @param [in] dim_1_in dimension 1 size
 * @param [in] dim_2_in dimensions 2 size
 */
KahanSum::KahanSum( int const dim_1_in, int const dim_2_in )
{
  // Sets the junk to 0 initially
  c_array_ = ArrayXXd::Zero(dim_1_in, dim_2_in) ;
  t_array_ = c_array_ ;
  y_array_ = c_array_ ;
}

/** Adds the term x_in to the sum with the Kahan Sum algorithm.
 *
 * @param [in,out] sum_inout
 * @param [in] x_in
 */
void KahanSum::addToSum( ArrayXXd &sum_inout, ArrayXXd const &x_in )
{
  // This is the algorithm - see Wikipedia article on Kahan Sum for details.
  y_array_ = x_in - c_array_ ;
  t_array_ = sum_inout + y_array_ ;
  c_array_ = (t_array_ - sum_inout) - y_array_ ;
  sum_inout = t_array_ ;
}

/** Resets the state of the Kahan Sum object - setting all the junk to zero.
 *
 */
void KahanSum::reset()
{
  c_array_ = 0.0 ;
  t_array_ = c_array_ ;
  y_array_ = c_array_ ;
}

/** Default constructor. Does nothing.
 *
 */
KahanSum::KahanSum()
{

}
