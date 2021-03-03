//
// Created by Tom on 15/12/2017.
//

#ifndef SEMICLASSICAL_MARKOVSTATEPARAMETERS_H
#define SEMICLASSICAL_MARKOVSTATEPARAMETERS_H

#ifdef HEADERONLY
  #include "Eigen/Dense"
#else
  #include <Eigen/Dense>
#endif
#include <iostream>
#include <vector>

//#ifdef HEADERONLY
//  #include "rapidjson/include/rapidjson/document.h"
//#else
//  #include "rapidjson/document.h"
//#endif
#include "rapidjson/include/rapidjson/document.h"

using namespace std;
using namespace Eigen;
using namespace rapidjson ;


class MarkovStateParameters
{
public:
  unsigned int getNumGroups() const;
  const ArrayXi &getArrayNumStates() const;
  const ArrayXi &getArrayNumSites() const;
  const vector<ArrayXi> &getSiteIndices() const;
  const vector<ArrayXXd> &getSiteCouplings() const;
  const vector<MatrixXd> &getHoppingRates() const;
  int getNumHoppingSteps() const;
  void readMarkovParameters( Document &document ) ;
  MarkovStateParameters();
  void setNumHoppingSteps( int num_hopping_steps_ );
private:
  // the number of groups of independent sites undergoing hops e.g. methyl groups
  unsigned int num_groups_ ;

  // the number of states that each site can hop between e.g. 3 for a methyl group
  ArrayXi array_num_states_ ;

  // the number of sites (hfccs) hopped between in each group
  ArrayXi array_num_sites_ ;

  // the indices of the sites that are hopped between for each independent group of sites
  vector<ArrayXi> site_indices_ ;

  // site couplings for each
  vector<ArrayXXd> site_couplings_ ;

  // matrix of hopping rates between different sites for each group of independent sites
  vector<MatrixXd> hopping_rates_ ;

  // the number of time steps per spin dynamics step
  int num_hopping_steps_ ;

  // Reads a group's parameters
  void readGroupParameters( Document &document , int const i);
};


#endif //SEMICLASSICAL_MARKOVSTATEPARAMETERS_H
