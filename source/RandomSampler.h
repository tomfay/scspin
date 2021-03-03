
#ifndef CLASSICAL_RANDOMSAMPLER_H
#define CLASSICAL_RANDOMSAMPLER_H

#define SOBOL_MAX_DIM 4
#define MT_SEED 12345
#define SOBOL_SEED 12345


#include <iostream>
#include <random>
#include <cmath>
#ifdef HEADERONLY
  #include "Eigen/Dense"
#else
  #include <Eigen/Dense>
#endif
#include "sobol.h"


using namespace std;
using namespace Eigen;

/** Class wrapping some of the standard library random stuff.
 *
 * Uses the Mersenne Twister as a pseudo random number generator. Can generate gaussian and uniform uniform random
 * variables and random unit 3-vectors sampled from the surface of a sphere.
 *
 */
class RandomSampler
{
public:
  // Default constructor. Sets up the mersenne twister generators with a random seed.
  RandomSampler();

  // Generates a random number sampled from a normal distribution.
  double sampleGaussianVariable( double const standard_deviation, double const mean) ;

  // Generates a random number
  double sampleUniformVariable( double const lower_bound, double const upper_bound ) ;
  Vector3d sampleUnitSphere() ;

  // Generates a random sobol vector
  ArrayXd sampleUniformVariablesSobol( int const size );

  // Samples two spheres using a sobol sequence
  VectorXd sampleTwoSpheresSobol() ;
  void sampleTwoSpheresSobol( VectorXd &vector ) ;

  // Resets the seed
  void resetSeed();

  // Gets the seed
  int getSeed() ;

  // Sets the seed
  void setSeed(int seed_in) ;


private:
  // The mersenne twister pseudo random number generator from the std random library.
  mt19937 mersenne_twister_generator_;
  // Normal distrubition and uniform distribution objects from the std random library.
  normal_distribution<double> normal_dist_;
  uniform_real_distribution<double> uniform_dist_;
  // Sobol junk for quasi-random number generation
  long long int seed_sobol_ ;
  double sobol_vector_[SOBOL_MAX_DIM] ;


};

#endif //CLASSICAL_RANDOMSAMPLER_H
