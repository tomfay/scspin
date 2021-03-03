//
// Created by Tom on 14/12/2017.
//

#include "SimulationOptions.h"

SimulationOptions::SimulationOptions()
{
  run_diff_sim_ = false ;
  run_hopping_ = false ;
  run_aniso_sim_ = false ;
  run_sw_sim_ = false ;
  run_vector_only_sim_ = false ;
  run_alt_norm_sim_ = false ;
  run_new_lengths_sim_ = false ;
}

void SimulationOptions::parseJson( string const json_filename )
{
  // read the file data in to a string
  ifstream json_file_stream(json_filename) ;
  stringstream buffer ;
  buffer << json_file_stream.rdbuf() ;

  // parse the string using rapidjson
  sim_document_.Parse( buffer.str().c_str() ) ;

  if (sim_document_.HasMember("Difference Simulation"))
  {
    readDiffSimOptions() ;
    run_diff_sim_ = true ;

    if (sim_document_["Difference Simulation"].HasMember("Hopping"))
    {
      readHoppingOptions() ;
      run_hopping_ = true ;
    }
  }

  if (sim_document_.HasMember("Anisotropic Simulation"))
  {
    readAnisoSimOptions() ;
  }

  if (sim_document_.HasMember("Simulation"))
  {
    readSimOptions() ;
    run_aniso_sim_ = false ;
  }

}


void SimulationOptions::readDiffSimOptions()
{
  Value &diff_sim_options = sim_document_["Difference Simulation"] ;

  // sets the number of spin system parameters to be read in
  num_systems_ = 2 ;
  spin_system_params_.resize(num_systems_) ;

  // sets the time step, number of samples and number of time steps
  time_step_ = diff_sim_options["Step Size"].GetDouble() ;
  num_samples_ = diff_sim_options["Samples"].GetInt() ;
  num_steps_ = diff_sim_options["Time Steps"].GetInt() ;

  // sets the spin system parameters for the two systems
  // creates a separate document to hold spin system info
  Document spin_system_doc ;
  // sets this as an object
  spin_system_doc.SetObject() ;
  // copies the system 1 radicals data into this new document
  spin_system_doc.CopyFrom(sim_document_["Difference Simulation"]["System 1"] , spin_system_doc.GetAllocator()) ;
  spin_system_params_[0].readSpinParameters( spin_system_doc ) ;
  spin_system_doc.RemoveAllMembers();
  spin_system_doc.CopyFrom(sim_document_["Difference Simulation"]["System 2"] , spin_system_doc.GetAllocator()) ;

  // sets the spin system 2 parameters from this document
  spin_system_params_[1].readSpinParameters( spin_system_doc ) ;

  // checks to see if hopping stuff is present and if so reads hopping parameters
  if (sim_document_["Difference Simulation"].HasMember("Hopping"))
  {
    run_hopping_ = true ;
  }

}

bool SimulationOptions::isDiffSim() const
{
  return run_diff_sim_;
}

int SimulationOptions::getNumSamples() const
{
  return num_samples_;
}

int SimulationOptions::getNumSteps() const
{
  return num_steps_;
}

double SimulationOptions::getTimeStep() const
{
  return time_step_;
}

Parameters SimulationOptions::getSpinParameters( int const system_index ) const
{
  return spin_system_params_[system_index];
}

void SimulationOptions::readHoppingOptions()
{

  Value &hopping_params = sim_document_["Difference Simulation"]["Hopping"] ;
  int radical_index ;
  int num_radicals = sim_document_["Difference Simulation"]["Hopping"].Size() ;
  num_radicals = min(2,num_radicals) ;
  markov_state_params_.resize(2);
  for (int i = 0 ; i < num_radicals ; i++)
  {
    radical_index = hopping_params[i]["Radical Index"].GetInt() ;
    markov_state_params_[radical_index].setNumHoppingSteps( hopping_params[i]["Time Steps"].GetInt() ) ;
    // creates a separate document to hold spin system info
    Document hopping_doc ;
    // sets this as an object
    hopping_doc.SetObject() ;
    // copies the system 1 radicals data into this new document
    hopping_doc.CopyFrom(sim_document_["Difference Simulation"]["Hopping"][i] , hopping_doc.GetAllocator()) ;
    // reads in parameters form new document
    markov_state_params_[radical_index].readMarkovParameters( hopping_doc ) ;
  }
}

bool SimulationOptions::isWithHopping() const
{
  return run_hopping_;
}

MarkovStateParameters SimulationOptions::getMarkovParameters( int const radical_index ) const
{
  return markov_state_params_[radical_index];
}

bool SimulationOptions::isAnisoSim() const
{
  return run_aniso_sim_;
}

bool SimulationOptions::isSwSim() const
{
  return run_sw_sim_;
}

void SimulationOptions::readAnisoSimOptions()
{
  Value &diff_sim_options = sim_document_["Anisotropic Simulation"] ;

  // sets the number of spin system parameters to be read in
  num_systems_ = 1 ;
  spin_system_params_.resize(num_systems_) ;

  // sets the time step, number of samples and number of time steps
  run_aniso_sim_ = true ;

  if (diff_sim_options.HasMember("Schulten Wolynes"))
  {
    run_sw_sim_ = diff_sim_options["Schulten Wolynes"].GetBool() ;
  }

  if (diff_sim_options.HasMember("Vectors Only"))
  {
    run_vector_only_sim_ = diff_sim_options["Vectors Only"].GetBool() ;
  }

  if (diff_sim_options.HasMember("Alt Norm"))
  {
    run_alt_norm_sim_ = diff_sim_options["Alt Norm"].GetBool() ;
  }

  if (diff_sim_options.HasMember("New Lengths"))
  {
    run_new_lengths_sim_ = diff_sim_options["New Lengths"].GetBool() ;
  }

  time_step_ = diff_sim_options["Step Size"].GetDouble() ;
  num_samples_ = diff_sim_options["Samples"].GetInt() ;
  num_steps_ = diff_sim_options["Time Steps"].GetInt() ;

  // sets the spin system parameters for the two systems
  // creates a separate document to hold spin system info
  Document spin_system_doc ;
  // sets this as an object
  spin_system_doc.SetObject() ;
  // copies the system 1 radicals data into this new document
  spin_system_doc.CopyFrom(sim_document_["Anisotropic Simulation"]["System 1"] , spin_system_doc.GetAllocator()) ;
  spin_system_params_[0].readSpinParameters( spin_system_doc ) ;
  spin_system_doc.RemoveAllMembers();

}

void SimulationOptions::readSimOptions()
{
  Value &diff_sim_options = sim_document_["Simulation"] ;

  // sets the number of spin system parameters to be read in
  num_systems_ = 1 ;
  spin_system_params_.resize(num_systems_) ;

  // sets the time step, number of samples and number of time steps
  time_step_ = diff_sim_options["Step Size"].GetDouble() ;
  num_samples_ = diff_sim_options["Samples"].GetInt() ;
  num_steps_ = diff_sim_options["Time Steps"].GetInt() ;

  // sets the spin system parameters for the two systems
  // creates a separate document to hold spin system info
  Document spin_system_doc ;
  // sets this as an object
  spin_system_doc.SetObject() ;
  // copies the system 1 radicals data into this new document
  spin_system_doc.CopyFrom(sim_document_["Simulation"]["System 1"] , spin_system_doc.GetAllocator()) ;
  spin_system_params_[0].readSpinParameters( spin_system_doc ) ;
  spin_system_doc.RemoveAllMembers();

}

bool SimulationOptions::isVectorOnlySim() const
{
  return run_vector_only_sim_;
}

bool SimulationOptions::isAltNormSim() const
{
  return run_alt_norm_sim_;
}

bool SimulationOptions::isNewLengthsSim() const
{
  return run_new_lengths_sim_;
}
