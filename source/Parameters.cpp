//
// Created by Tom on 07/12/2017.
//

#include "Parameters.h"

void Parameters::readSpinParameters( string const json_filename )
{
  // read the file data in to a string
  ifstream json_file_stream(json_filename) ;
  stringstream buffer ;
  buffer << json_file_stream.rdbuf() ;

  // parse the string using rapidjson
  Document document ;
  document.Parse( buffer.str().c_str() ) ;


  // Read the Radicals section if present
  if (document.HasMember("Radicals"))
  {
    readRadicals( document );
    //cout << "Radical info read." << endl ;
  }

  // Read the Electron Couplings segment if present
  if (document.HasMember("Electron Couplings"))
  {
    readElectronCouplings( document ) ;
    //cout << "Coupling info read." << endl ;
  }

  // Read the Rate Constants section if present
  if (document.HasMember("Rate Constants"))
  {
    readRateConstants( document ) ;
    //cout << "Rate constant info read." << endl ;
  }

  // states that parameters have been read in
  parameters_read_ = true ;

}

void Parameters::readSpinParameters( Document &document )
{

  // Read the Radicals section if present
  if (document.HasMember("Radicals"))
  {
    readRadicals( document );
    //cout << "Radical info read." << endl ;
  }

  // Read the Electron Couplings segment if present
  if (document.HasMember("Electron Couplings"))
  {
    readElectronCouplings( document ) ;
    //cout << "Coupling info read." << endl ;
  }

  // Read the Rate Constants section if present
  if (document.HasMember("Rate Constants"))
  {
    readRateConstants( document ) ;
    //cout << "Rate constant info read." << endl ;
  }


  // states that parameters have been read in
  parameters_read_ = true ;

//  cout << "electron coupling tensor" << endl ;
//  cout <<  e_coupling_tensor_ << endl ;
//  cout << "aniso couplings 1" << endl ;
//  cout << aniso_nuc_couplings_[0] << endl ;
//  cout << "aniso couplings 2" << endl ;
//  cout << aniso_nuc_couplings_[1] << endl ;
//    cout << "vec lengths 1" << endl ;
//  cout << nuc_vec_lengths_[0] << endl ;
//  cout << "vec lengths 2" << endl ;
//  cout << nuc_vec_lengths_[1] << endl ;

}

/** Reads in parameters for the radicals - currently assumes only up to 2 radicals are specified
 *
 * @param document
 */
void Parameters::readRadicals( Document &document )
{
  int radical_index ;
  double x ;
  Value &radicals = document["Radicals"] ;
  cout << "Reading radical parameters." << endl ;
  //cout << "Number of radicals instances: " << radicals.Size() << endl ;
  for(int i = 0 ; i < radicals.Size() ; i++)
  {

    if (radicals[i].HasMember("Index") && (radicals[i]["Index"].GetInt()==0 || radicals[i]["Index"].GetInt()==1))
    {
      radical_index = radicals[i]["Index"].GetInt();
      cout << "Reading radical " << radical_index << " data." << endl ;
      //cout << "radical index: " << radical_index << endl ;
      num_nuc_(radical_index) = 0 ;
      if (radicals[i].HasMember("Omega"))
      {
        cout << "Reading omega." << endl ;
//        omegas_[radical_index] = Vector3d( radicals[i]["Omega"][0].GetDouble(), radicals[i]["Omega"][1].GetDouble(), radicals[i]["Omega"][2].GetDouble() );
        //cout << "omega: "<< omegas_[radical_index] << endl ;
        x = radicals[i]["Omega"][0].GetDouble() ;
        omegas_[radical_index](0) = x ;
        x = radicals[i]["Omega"][1].GetDouble() ;
        omegas_[radical_index](1) = x ;
        x = radicals[i]["Omega"][2].GetDouble() ;
        omegas_[radical_index](2) = x ;

      }
      if (radicals[i].HasMember("Hyperfines"))
      {
        cout << "Reading hyperfines." << endl ;
        num_nuc_(radical_index) = radicals[i]["Hyperfines"]["Number"].GetInt() ;
        //cout << "num nuc: " << num_nuc_(radical_index) << endl ;
        if (num_nuc_(radical_index)>0)
        {

          nuc_couplings_[radical_index].resize(num_nuc_(radical_index)) ;
          aniso_nuc_couplings_[radical_index].resize(9,num_nuc_(radical_index)) ;
          nuc_vec_lengths_[radical_index].resize(num_nuc_(radical_index)) ;

          if(radicals[i]["Hyperfines"].HasMember("Isotropic Constants"))
          {
            cout << "Reading isotropic hyperfines." << endl ;
            for (int j = 0 ; j< num_nuc_(radical_index) ;j++)
            {
              nuc_couplings_[radical_index](j) = radicals[i]["Hyperfines"]["Isotropic Constants"][j].GetDouble() ;

            }
          }

          if(radicals[i]["Hyperfines"].HasMember("Tensors"))
          {
            cout << "Reading anisotropic hyperfine tensors." << endl ;
            aniso_params_set_ = true ;
            for (int j = 0 ; j< num_nuc_(radical_index) ;j++)
            {
              for (int k = 0 ; k < 9 ; k++)
              {
                aniso_nuc_couplings_[radical_index](k,j) = radicals[i]["Hyperfines"]["Tensors"][j][k].GetDouble() ;
//                cout << aniso_nuc_couplings_[radical_index](k,j) << " ";
              }
//              cout << endl ;
              if(!checkTensorSymmetric(aniso_nuc_couplings_[radical_index].col(j)))
              {
                cout << "Hyperfine tensor not symmetric! Radical: " <<radical_index << "Nucleus: " << j <<endl;
              }
            }
          }

          if(radicals[i]["Hyperfines"].HasMember("Multiplicities"))
          {
            cout << "Reading nuclear spin multiplicities." << endl ;
            for (int j = 0 ; j< num_nuc_(radical_index) ;j++)
            {
              x = radicals[i]["Hyperfines"]["Multiplicities"][j].GetDouble() ;
              x = 0.5*(x-1.0) ;
              x = sqrt(x*(x+1.0)) ;
              nuc_vec_lengths_[radical_index](j) = x ;
            }
          }
          else
          {
            nuc_vec_lengths_[radical_index] = sqrt(0.5*1.5)*ArrayXd::Ones(num_nuc_(radical_index)) ;
          }

        }

      }
    }
  }
}

Parameters::Parameters()
{
  omegas_.resize(2) ;
  omegas_[0] = Vector3d::Zero() ;
  omegas_[1] = Vector3d::Zero() ;
  nuc_couplings_.resize(2);
  aniso_nuc_couplings_.resize(2);
  nuc_vec_lengths_.resize(2);
  num_nuc_ = ArrayXi::Zero(2);
  exchange_coupling_ = 0.0 ;
  triplet_rate_ = 0.0 ;
  singlet_rate_ = 0.0 ;
  dephasing_rate_ = 0.0 ;
  parameters_read_ = false ;
  aniso_params_set_ = false ;
  gamma_relaxation_rates_ = ArrayXd::Zero(6) ;
}

void Parameters::readElectronCouplings( Value &electron_couplings )
{
  if(electron_couplings.HasMember("Exchange"))
  {
    exchange_coupling_ = electron_couplings["Exchange"].GetDouble() ;
  }
  else
  {
    exchange_coupling_ = 0.0 ;
  }

}

void Parameters::readElectronCouplings( Document &document )
{
  cout << "Reading electron couplings." << endl ;
  if(document["Electron Couplings"].HasMember("Exchange"))
  {

    cout << "Reading isotropic, scalar coupling." << endl ;
    exchange_coupling_ = document["Electron Couplings"]["Exchange"].GetDouble() ;
  }
  if(document["Electron Couplings"].HasMember("Tensor"))
  {
    cout << "Reading electron coupling tensor." << endl ;
    aniso_params_set_ = true ;
    for (int j = 0  ; j < 3 ; j++)
    {
      for (int k = 0 ; k < 3 ; k++)
      {
        e_coupling_tensor_(k, j) = document["Electron Couplings"]["Tensor"][3*j+k].GetDouble() ;
      }
    }
    ArrayXd tensor_flat = ArrayXd::Zero(9) ;
    tensor_flat.segment(0,3) = e_coupling_tensor_.row(0) ;
    tensor_flat.segment(3,3) = e_coupling_tensor_.row(1) ;
    tensor_flat.segment(6,3) = e_coupling_tensor_.row(2) ;
    if(!checkTensorSymmetric(tensor_flat))
    {
      cout << "Electron Coupling Tensor not symmetric! "<<endl;
    }

  }
}

const ArrayXi &Parameters::getNumNuc() const
{
  return num_nuc_;
}

const vector<ArrayXd> &Parameters::getNucCouplings() const
{
  return nuc_couplings_;
}


const vector<Vector3d> &Parameters::getOmegas() const
{
  return omegas_;
}

double Parameters::getExchangeCoupling() const
{
  return exchange_coupling_;
}

double Parameters::getTripletRate() const
{
  return triplet_rate_;
}

double Parameters::getSingletRate() const
{
  return singlet_rate_;
}

double Parameters::getDephasingRate() const
{
  return dephasing_rate_;
}

void Parameters::readRateConstants( Document &document )
{
  cout << "Reading rate constants." << endl ;
  Value& rate_constants = document["Rate Constants"] ;
  if (rate_constants.HasMember("Singlet Rate"))
  {
    cout << "Reading singlet rate." << endl ;
    singlet_rate_ = rate_constants["Singlet Rate"].GetDouble() ;
  }
  if (rate_constants.HasMember("Triplet Rate"))
  {
    cout << "Reading triplet rate." << endl ;
    triplet_rate_ = rate_constants["Triplet Rate"].GetDouble() ;
  }
  if ( rate_constants.HasMember("Symmetric Rate")
    && !(rate_constants.HasMember("Singlet Rate") && rate_constants.HasMember("Triplet Rate")))
  {
    cout << "Reading symmetric rate." << endl ;
    singlet_rate_ = rate_constants["Symmetric Rate"].GetDouble() ;
    triplet_rate_ = singlet_rate_ ;
  }
  if (rate_constants.HasMember("Dephasing Rate"))
  {
    cout << "Reading singlet-triplet dephasing rate." << endl ;
    dephasing_rate_ = rate_constants["Dephasing Rate"].GetDouble() ;
  }
  if (rate_constants.HasMember("Spin Relaxation Rates"))
  {
    cout << "Reading spin relaxation rates." << endl ;
    for (int i = 0 ;  i < 6 ; i++)
    {
      gamma_relaxation_rates_(i) = rate_constants["Spin Relaxation Rates"][i].GetDouble() ;
    }
  }
}

const Vector3d Parameters::getOmega( int const i ) const
{
  return omegas_[i];
}

const ArrayXd Parameters::getNucCouplings( int const i ) const
{
  return nuc_couplings_[i];
}

bool Parameters::containsData() const
{
  return parameters_read_;
}

bool Parameters::anisoParamsSet() const
{
  return aniso_params_set_;
}

const ArrayXXd Parameters::getAnisoNucCouplings( int const i ) const
{
  return aniso_nuc_couplings_[i];
}

const Array33d Parameters::getECouplingTensor() const
{
  return e_coupling_tensor_ ;
}

const ArrayXd Parameters::getNucVecLengths( int const i ) const
{
  return nuc_vec_lengths_[i];
}


bool Parameters::checkTensorSymmetric(ArrayXd tensor_in)
{
  bool is_symmetric= true ;
  if( tensor_in( 1) != tensor_in( 3))
  {
    is_symmetric = false ;
  }
  if( tensor_in( 2) != tensor_in( 6))
  {
    is_symmetric = false ;
  }
  if( tensor_in( 5) != tensor_in( 7))
  {
    is_symmetric = false ;
  }
  return is_symmetric;
}

const ArrayXd Parameters::getGammaRelaxationRates() const
{
  return gamma_relaxation_rates_ ;
}


