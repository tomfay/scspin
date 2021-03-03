
#include "RandomSampler.h"

/** Constructor for the RandomSampler object.
 *
 * Uses a random seed for the mersenne twister generator.
 */
RandomSampler::RandomSampler()
{
//  random_device seed_generator ;
//  mersenne_twister_generator_.seed(seed_generator()) ;
  mersenne_twister_generator_.seed(MT_SEED) ;
  normal_dist_ = normal_distribution<double>(0.0,1.0) ;
  uniform_dist_ = uniform_real_distribution<double>(0.0,1.0);

  seed_sobol_ = SOBOL_SEED ;
}

/** Samples a random variable sampled from a normal distribution.
 *
 * @param [in] standard_deviation standard deviation of the distribution to be sampled from.
 * @param [in] mean mean of the distribution to be sampled from.
 * @return the sampled variable
 */
double RandomSampler::sampleGaussianVariable(double const standard_deviation, double const mean)
{
  return standard_deviation * normal_dist_(mersenne_twister_generator_) + mean ;
}

/** Samples a random variable from a uniform distribution.
 *
 * @param [in] lower_bound lower bound for the distribution
 * @param [in] upper_bound upper bound for the distribution
 * @return uniformly sampled random variable
 */
double RandomSampler::sampleUniformVariable(double const lower_bound, double const upper_bound)
{
  return (upper_bound - lower_bound) * uniform_dist_(mersenne_twister_generator_) + lower_bound ;
}

/** Samples a random unit 3-vector uniformly from the surface of a sphere.
 *
 * @return randomly sampled unit 3-vector from surface of sphere.
 */
Vector3d RandomSampler::sampleUnitSphere()
{
  double theta = acos( sampleUniformVariable( -1.0 , 1.0 ) ) ;
  double phi = sampleUniformVariable( 0.0, 2.0 * M_PI ) ;
  Vector3d unit_vector( sin(theta) * cos(phi) , sin(theta) * sin(phi)  , cos(theta) ) ;
  return unit_vector;
}

void RandomSampler::resetSeed()
{
  random_device seed_generator ;
  mersenne_twister_generator_.seed(seed_generator()) ;
}

/** Samples real numbers from [0,1) using the sobol quasi-random sequence
 *
 * @param [in] size size of the sobol vector to be generated
 * @return an array of quasi-random numbers
 */
ArrayXd RandomSampler::sampleUniformVariablesSobol( int const size )
{
  ArrayXd sobol_vector_out = ArrayXd::Zero( size );
  // generates a sobol vector
  i8_sobol( size , &seed_sobol_ , sobol_vector_) ;
  // copies into an output Eigen::ArrayXd
  for (int i = 0 ; i < size ; i++)
  {
    sobol_vector_out(i) = sobol_vector_[i] ;
  }
  return sobol_vector_out;
}

VectorXd RandomSampler::sampleTwoSpheresSobol()
{
  VectorXd vectors_out = VectorXd::Zero(6) ;
  // sample the angles
  ArrayXd sobol_vector = sampleUniformVariablesSobol(4) ;
  double theta,phi ;
  for (int i = 0 ;  i < 2 ; i++)
  {
    // theta_i
    theta = acos( 2.0 * sobol_vector(2*i) - 1.0 ) ;
    // phi_i
    phi = 2.0 * M_PI * sobol_vector(2*i+1) ;
    vectors_out( 3*i ) = sin(theta) * cos(phi) ;
    vectors_out( 3*i + 1 ) = sin(theta) * sin(phi) ;
    vectors_out( 3*i + 2 ) = cos(theta) ;
  }

  return vectors_out;
}

void RandomSampler::sampleTwoSpheresSobol( VectorXd &vectors )
{
  // sample the angles
  ArrayXd sobol_vector = sampleUniformVariablesSobol(4) ;
  double theta,phi ;
  for (int i = 0 ;  i < 2 ; i++)
  {
    // theta_i
    theta = acos( 2.0 * sobol_vector(2*i) - 1.0 ) ;
    // phi_i
    phi = 2.0 * M_PI * sobol_vector(2*i+1) ;
    vectors( 3*i ) = sin(theta) * cos(phi) ;
    vectors( 3*i +1 ) = sin(theta) * sin(phi) ;
    vectors( 3*i +2 ) = cos(theta) ;
  }
}

