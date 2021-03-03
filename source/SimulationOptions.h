//
// Created by Tom on 14/12/2017.
//

#ifndef SEMICLASSICAL_SIMULATIONOPTIONS_H
#define SEMICLASSICAL_SIMULATIONOPTIONS_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "Parameters.h"
#include "MarkovStateParameters.h"
//#ifdef HEADERONLY
//#include "rapidjson/include/rapidjson/document.h"
//#else
//#include "rapidjson/document.h"
//#endif
#include "rapidjson/include/rapidjson/document.h"
using namespace std ;

class SimulationOptions
{
public:
  void parseJson( string const json_filename ) ;
  SimulationOptions() ;
  bool isDiffSim() const;
  bool isWithHopping() const ;
  bool isAnisoSim() const ;
  bool isSwSim() const ;
  bool isVectorOnlySim() const ;
  bool isAltNormSim() const ;
  bool isNewLengthsSim() const ;


  int getNumSamples() const;
  int getNumSteps() const;
  double getTimeStep() const;
  Parameters getSpinParameters( int const system_index) const ;
  MarkovStateParameters getMarkovParameters( int const radical_index ) const ;

private:
  bool run_diff_sim_ ;
  bool run_hopping_ ;
  bool run_aniso_sim_ ;
  bool run_sw_sim_ ;
  bool run_vector_only_sim_ ;
  bool run_alt_norm_sim_ ;
  bool run_new_lengths_sim_ ;
  int num_samples_ ;
  int num_steps_ ;
  double time_step_ ;
  int num_systems_ ;
  vector<Parameters> spin_system_params_ ;
  vector<MarkovStateParameters> markov_state_params_ ;

  Document sim_document_ ;

  void readDiffSimOptions() ;
  void readHoppingOptions() ;
  void readAnisoSimOptions() ;
  void readSimOptions() ;




};


#endif //SEMICLASSICAL_SIMULATIONOPTIONS_H
