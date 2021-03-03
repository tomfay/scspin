//
// Created by Tom on 07/12/2017.
//

#ifndef SEMICLASSICAL_PARAMETERS_H
#define SEMICLASSICAL_PARAMETERS_H

#include <iostream>
#ifdef HEADERONLY
  #include "Eigen/Dense"
#else
  #include <Eigen/Dense>
#endif
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

//#ifdef HEADERONLY
//#include "rapidjson/include/rapidjson/document.h"
//#else
//#include "rapidjson/document.h"
//#endif
#include "rapidjson/include/rapidjson/document.h"

using namespace std ;
using namespace rapidjson ;
using namespace Eigen ;

class Parameters
{
public:
  void readSpinParameters(string const json_filename) ;
  void readSpinParameters(Document &document) ;
  void readRadicals( Document &document ) ;
  void readElectronCouplings( Document &document ) ;
  void readRateConstants( Document &document ) ;
  Parameters() ;
  const ArrayXi &getNumNuc() const;
  const vector<ArrayXd> &getNucCouplings() const;
  const vector<ArrayXXd> &getAnisoNucCouplings() const;
  const vector<Vector3d> &getOmegas() const;
  double getExchangeCoupling() const;
  double getTripletRate() const;
  double getSingletRate() const;
  double getDephasingRate() const ;
  const ArrayXd getNucCouplings(int const i) const ;
  const ArrayXXd getAnisoNucCouplings( int const i ) const ;
  const Array33d getECouplingTensor() const ;
  const Vector3d getOmega(int const i) const ;
  const ArrayXd getNucVecLengths( int const i ) const ;
  const ArrayXd getGammaRelaxationRates() const ;
  bool containsData() const ;
  void readElectronCouplings( Value &electron_couplings );
  bool anisoParamsSet() const ;

private:
  ArrayXi num_nuc_ ;
  vector<ArrayXd> nuc_couplings_ ;
  vector<ArrayXXd> aniso_nuc_couplings_ ;
  vector<Vector3d> omegas_ ;
  vector<ArrayXd> nuc_vec_lengths_ ;
  double exchange_coupling_ ;
  double triplet_rate_ ;
  double singlet_rate_ ;
  double dephasing_rate_ ;
  ArrayXd gamma_relaxation_rates_ ;
  Array33d e_coupling_tensor_ ;
  bool parameters_read_ ;
  bool aniso_params_set_ ;
  bool checkTensorSymmetric(ArrayXd tensor_in);

};


#endif //SEMICLASSICAL_PARAMETERS_H
