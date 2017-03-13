#include <sde.h>
#include <string>
using namespace PJLib;

const double SDE::Pi=4.*atan(1.);


void LogNormalSDE::link_parameters(ParameterMaps *params,std::string filename)
{
  params->assignParameter(&sigma,filename+"-sigma","Volatility in "+filename+" :: ",0.2);
  params->assignParameter(&interest_rate,filename+"-interest_rate"," Interest rate "+filename+" :: ",0.05);
  params->assignParameter(&dividend,filename+"-dividend"," Dividend "+filename+" :: ",0.);
  process_name = filename;
}

void LogNormalSDE::erase_parameters()
{
  sde_params->eraseParameter(process_name+"-sigma");
  sde_params->eraseParameter(process_name+"-interest_rate");
  sde_params->eraseParameter(process_name+"-dividend");
}

void MRLogSinSDE::link_parameters(ParameterMaps *params,std::string filename)
{
  params->assignParameter(&theta,filename+"-theta","Long term average in "+filename+" :: ",10.);
  params->assignParameter(&mu,filename+"-mu","Mean reversion rate in "+filename+" :: ",0.1);
  params->assignParameter(&sigma,filename+"-sigma","Volatility in "+filename+" :: ",0.2);
  params->assignParameter(&psi,filename+"-psi","Phase shift in "+filename+" :: ",2.);
  params->assignParameter(&phi,filename+"-phi","Amplitude in "+filename+" :: ",0.375);
  params->assignParameter(&gamma,filename+"-gamma","Frequency in "+filename+" :: ",2.*Pi/24.);
  process_name = filename;
}

void MRLogSinSDE::erase_parameters()
{
  sde_params->eraseParameter(process_name+"-theta");
  sde_params->eraseParameter(process_name+"-mu");
  sde_params->eraseParameter(process_name+"-sigma");
  sde_params->eraseParameter(process_name+"-psi");
  sde_params->eraseParameter(process_name+"-phi");
  sde_params->eraseParameter(process_name+"-gamma");
}

