//
// Created by Tom on 15/12/2017.
//

#include "MarkovStateParameters.h"

void MarkovStateParameters::readMarkovParameters( Document &document )
{
  if( document.HasMember("Groups"))
  {
    num_groups_ = document["Groups"].Size() ;
    site_indices_.resize(num_groups_) ;
    site_couplings_.resize(num_groups_) ;
    hopping_rates_.resize(num_groups_) ;
    array_num_sites_.resize(num_groups_) ;
    array_num_states_.resize(num_groups_) ;

    for (int i = 0 ; i < num_groups_ ; i++)
    {
      readGroupParameters( document , i );
//      cout << "Site indices " << i << endl;
//      cout << site_indices_[i]<< endl ;
//      cout << "Site Couplings " << i << endl;
//      cout << site_couplings_[i]<< endl ;
//      cout << "Hopping rates " << i << endl ;
//      cout << hopping_rates_[i]<< endl ;
    }
  }
}

void MarkovStateParameters::readGroupParameters( Document &document, int const i )
{
  Value &group = document["Groups"][i] ;

  // extract the number of sites
  array_num_sites_(i) = group["Site Indices"].Size() ;

  // extract the site indices
  site_indices_[i].resize(array_num_sites_(i)) ;
  for (int k = 0 ; k < array_num_sites_(i) ;  k++)
  {
    site_indices_[i](k) = group["Site Indices"][k].GetInt();
  }

  // extract the number of states
  array_num_states_(i) = group["Hopping Rates"].Size() ;

  // extract the hopping rates
  hopping_rates_[i].resize( array_num_states_(i) , array_num_states_(i) ) ;
  for (int j = 0 ; j < array_num_states_(i) ; j++)
  {
    for (int k = 0 ; k < array_num_states_(i) ; k++)
    {
      hopping_rates_[i](j,k) = group["Hopping Rates"][j][k].GetDouble();
    }
  }

  // extract the site couplings
  site_couplings_[i].resize( array_num_sites_(i) , array_num_states_(i) ) ;
  for (int j = 0; j < array_num_sites_(i) ; j++)
  {
    for (int k = 0 ; k < array_num_states_(i) ; k++)
    {
      site_couplings_[i](j,k) = group["Site Couplings"][j][k].GetDouble();
    }
  }
}

MarkovStateParameters::MarkovStateParameters()
{
  num_groups_ = 0 ;
  num_hopping_steps_ = 0 ;
}

void MarkovStateParameters::setNumHoppingSteps( int num_hopping_steps_ )
{
  MarkovStateParameters::num_hopping_steps_ = num_hopping_steps_;
}

unsigned int MarkovStateParameters::getNumGroups() const
{
  return num_groups_;
}

const ArrayXi &MarkovStateParameters::getArrayNumStates() const
{
  return array_num_states_;
}

const ArrayXi &MarkovStateParameters::getArrayNumSites() const
{
  return array_num_sites_;
}

const vector<ArrayXi> &MarkovStateParameters::getSiteIndices() const
{
  return site_indices_;
}

const vector<ArrayXXd> &MarkovStateParameters::getSiteCouplings() const
{
  return site_couplings_;
}

const vector<MatrixXd> &MarkovStateParameters::getHoppingRates() const
{
  return hopping_rates_;
}

int MarkovStateParameters::getNumHoppingSteps() const
{
  return num_hopping_steps_;
}
