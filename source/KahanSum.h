//
// Created by Tom on 05/12/2017.
//

#ifndef SEMICLASSICAL_KAHANSUM_H
#define SEMICLASSICAL_KAHANSUM_H

#ifdef HEADERONLY
  #include "Eigen/Dense"
#else
  #include <Eigen/Dense>
#endif
#include <iostream>

using namespace std ;
using namespace Eigen ;

/** A class for performing a Kahan sum on a 2-dimensional Eigen array.
 *
 * This is an algorithm for reducing loss of accuracy when adding a very small number to a very large number.
 */
class KahanSum
{
public:
  // Constructor with dimensions of the array to be summed.
  KahanSum(int const dim_1_in, int const dim_2_in);

  // Default constructor
  KahanSum() ;

  // Adds a value to the target sum
  void addToSum( ArrayXXd &sum_inout , ArrayXXd const &x_in );

  // Resets the algorithm, setting c = 0
  void reset();
private:
  // Arrays of junk used in the algorithm
  ArrayXXd c_array_, y_array_, t_array_ ;
};


#endif //SEMICLASSICAL_KAHANSUM_H
