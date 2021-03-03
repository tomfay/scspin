#include <iostream>

#include "SpinStateEvolver.h"
#include "Simulation.h"
#include "SimulationOptions.h"
#include "MolecularState.h"
#include "Parameters.h"

using namespace std ;
using namespace Eigen ;

int main( int argc, char* argv[] )
{
  if(argc != 2)
  {
    cerr << "Incorrect input string." << endl;
    return 1;
  }


  // creates an object that reads in the simulation options and parameters
  SimulationOptions sim_options ;
  string json_filename(argv[1]);
  sim_options.parseJson(json_filename) ;

  // creates an object that performs the SC simulations
  Simulation sim ;

  // if the simulation is a difference simulation runs this with parameters from the sim_options read in
  if (sim_options.isDiffSim() && (!sim_options.isWithHopping()))
  {
    // runs the simulation with parameters defined in the simulation options object
    cout << "Running difference simulation." << endl ;
    sim.runDifferenceSimulation( sim_options.getNumSamples(), // number of Monte carlo samples
                                 sim_options.getNumSteps() ,  // number of time steps for evolution
                                 sim_options.getTimeStep() ,  // time step size for evolution
                                 sim_options.getSpinParameters(0),    // parameters for spin system 1
                                 sim_options.getSpinParameters(1) )  ; // parameters for spin system 2
  }
  else if (sim_options.isDiffSim() && sim_options.isWithHopping())
  {
    // runs the difference simulation with hopping as defined in the simulation with
    cout << "Running difference simulation with hopping." << endl ;
    sim.runDiffSimWithHopping( sim_options.getNumSamples(), // number of Monte carlo samples
                               sim_options.getNumSteps() ,  // number of time steps for evolution
                               sim_options.getTimeStep() ,  // time step size for evolution
                               sim_options.getSpinParameters(0),    // parameters for spin system 1
                               sim_options.getSpinParameters(1) ,   // parameters for spin system 2
                               sim_options.getMarkovParameters(0) ,
                               sim_options.getMarkovParameters(1) ) ;

    // runs the simulation with parameters defined in the simulation options object
//    cout << "Running difference simulation." << endl ;
//    sim.runDifferenceSimulation( sim_options.getNumSamples(), // number of Monte carlo samples
//                                 sim_options.getNumSteps() ,  // number of time steps for evolution
//                                 sim_options.getTimeStep() ,  // time step size for evolution
//                                 sim_options.getSpinParameters(0),    // parameters for spin system 1
//                                 sim_options.getSpinParameters(1) )  ; // parameters for spin system 2
  }
  else if ((!sim_options.isDiffSim())
    && (sim_options.isAnisoSim())
    && (!sim_options.isSwSim())
    && (!sim_options.isVectorOnlySim())
    && (!sim_options.isAltNormSim())
    && (!sim_options.isNewLengthsSim()))
  {
    // runs the difference simulation with hopping as defined in the simulation with
    cout << "Running anisotropic simulation" << endl ;
    sim.runAnisoSimulation( sim_options.getNumSamples(), // number of Monte carlo samples
                            sim_options.getNumSteps() ,  // number of time steps for evolution
                            sim_options.getTimeStep() ,  // time step size for evolution
                            sim_options.getSpinParameters(0)) ; // parameters for spin system
  }
  else if ((!sim_options.isDiffSim())
    && (sim_options.isAnisoSim())
    && (!sim_options.isSwSim())
    && (sim_options.isVectorOnlySim())
    && (!sim_options.isNewLengthsSim())
    && (!sim_options.isAltNormSim()))
  {
    // runs the difference simulation with hopping as defined in the simulation with
    cout << "Running anisotropic vector only simulation" << endl ;
    sim.runAnisoVectorOnlySimulation( sim_options.getNumSamples(), // number of Monte carlo samples
                            sim_options.getNumSteps() ,  // number of time steps for evolution
                            sim_options.getTimeStep() ,  // time step size for evolution
                            sim_options.getSpinParameters(0)) ; // parameters for spin system
  }
  else if ((!sim_options.isDiffSim())
    && (sim_options.isAnisoSim())
    && (sim_options.isSwSim())
    && (!sim_options.isVectorOnlySim())
    && (!sim_options.isNewLengthsSim())
    && (!sim_options.isAltNormSim()))
  {
    // runs the difference simulation with hopping as defined in the simulation with
    cout << "Running anisotropic Schulten Wolynes simulation" << endl ;
    sim.runAnisoSwSimulation( sim_options.getNumSamples(), // number of Monte carlo samples
                            sim_options.getNumSteps() ,  // number of time steps for evolution
                            sim_options.getTimeStep() ,  // time step size for evolution
                            sim_options.getSpinParameters(0)) ; // parameters for spin system
  }
  else if ((!sim_options.isDiffSim())
    && (sim_options.isAnisoSim())
    && (!sim_options.isSwSim())
    && (!sim_options.isVectorOnlySim())
    && (!sim_options.isNewLengthsSim())
    && (sim_options.isAltNormSim()))
  {
    // runs the difference simulation with hopping as defined in the simulation with
    cout << "Running anisotropic simulation with alternative normalisation" << endl ;
    sim.runAnisoAltNormSimulation( sim_options.getNumSamples(), // number of Monte carlo samples
                              sim_options.getNumSteps() ,  // number of time steps for evolution
                              sim_options.getTimeStep() ,  // time step size for evolution
                              sim_options.getSpinParameters(0)) ; // parameters for spin system
  }
  else if ((!sim_options.isDiffSim())
    && (sim_options.isAnisoSim())
    && (!sim_options.isSwSim())
    && (!sim_options.isVectorOnlySim())
    && (sim_options.isNewLengthsSim())
    && (!sim_options.isAltNormSim()))
  {
    // runs the difference simulation with hopping as defined in the simulation with
    cout << "Running anisotropic simulation with sqrt(5/12) vector lengths" << endl ;
    sim.runAnisoNewLengthsSimulation( sim_options.getNumSamples(), // number of Monte carlo samples
                                   sim_options.getNumSteps() ,  // number of time steps for evolution
                                   sim_options.getTimeStep() ,  // time step size for evolution
                                   sim_options.getSpinParameters(0)) ; // parameters for spin system
  }
  else if ((!sim_options.isDiffSim())&&(!sim_options.isAnisoSim()))
  {
    // runs the difference simulation with hopping as defined in the simulation with
    cout << "Running isotropic simulation" << endl ;
    sim.runSimulation( sim_options.getNumSamples(), // number of Monte carlo samples
                            sim_options.getNumSteps() ,  // number of time steps for evolution
                            sim_options.getTimeStep() ,  // time step size for evolution
                            sim_options.getSpinParameters(0)) ; // parameters for spin system
  }

  return 0;
}