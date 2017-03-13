//  ## Stochastic Differential Equation

#ifndef _sde_h_included_
#define _sde_h_included_

#include <string>
#include <cmath>
#include <param_maps.h>

namespace PJLib{

class SDE
{

  protected:
    
    static const double Pi;
  
  public:
	
	// constructor
	SDE(){};
	// destructor
	~SDE(){};
	// default is the log normal distribution
	// dX = alpha(X,t) dt + beta(X,t) dW
	virtual double XplusdX(double X,double t,double dt,double dW)
	{return X + alpha(X,t)*dt + beta(X,t)*dW;};
	// function - multiplying the dt
	virtual double alpha(double X,double t){return 0.;};
	virtual double alpha(int i,double t){return 0.;};
	virtual double alpha(double X,int k,double theta){return 0.;};
	// function - multiplying the dW
	virtual double beta(double X,double t){return 1.;};
	virtual double beta(int i,double t){return 1.;};
	virtual double beta(double X,int k,double theta){return 1.;};
	
};

class LogNormalSDE: public SDE
{
  std::string process_name;
  // link parameters to the outside
  ParameterMaps *sde_params;
  void link_parameters(ParameterMaps *params,std::string filename);
  void erase_parameters();
  protected:
    // parameter for distribution
    double sigma,interest_rate,dividend;
    
  public:
    // constructors
    LogNormalSDE(ParameterMaps *params,std::string filename){sde_params=params;link_parameters(params,filename);};
    ~LogNormalSDE(){erase_parameters();};
    // drift of process
    virtual double alpha(double S,double tau){return (interest_rate-dividend)*S;};
    // variation of process
    virtual double beta(double S,double tau){return sigma*S;};
    // return process name
    std::string name(){return process_name;};
    double sde_sigma(){return sigma;};
    double sde_interest_rate(){return interest_rate;};
    double sde_dividend(){return dividend;};
};

class MRLogSinSDE: public SDE
{
  std::string process_name;
  // link parameters to the outside
  ParameterMaps *sde_params;
  void link_parameters(ParameterMaps *params,std::string filename);
  void erase_parameters();
  protected:
    // parameter for distribution
    double theta,mu,sigma,psi,phi,gamma;
  public:
    // constructors
    MRLogSinSDE(ParameterMaps *params,std::string filename){sde_params=params;link_parameters(params,filename);};
    ~MRLogSinSDE(){erase_parameters();};
    // mean reversion function
    virtual double f(double tau){return theta*(1.+phi*sin(gamma*(tau-psi)));};
    // differential of the mean reversion function
    virtual double fdash(double tau){return theta*phi*gamma*cos(gamma*(tau-psi));};
    // drift of process
    virtual double alpha(double X,double tau){return mu*(f(tau)-X) - fdash(tau);};
    // variation of process
    virtual double beta(double X,double tau){return sigma*X;};
    // return process name
    std::string name(){return process_name;};
};

class MRSinSDE: public MRLogSinSDE
{
  public:
    // constructors
    MRSinSDE(ParameterMaps *params,std::string filename):MRLogSinSDE(params,filename){};
    // variation of process
    double beta(double X,double tau){return sigma;};
};

class MRLSForwardSDE: public MRLogSinSDE
{
  public:
    // constructors
    MRLSForwardSDE(ParameterMaps *params,std::string filename):MRLogSinSDE(params,filename){};
    // variation of process
    double f(double t){return theta*(1.+phi*sin(Pi+gamma*(t-psi)));};
    double fdash(double tau){return theta*phi*gamma*cos(Pi+gamma*(tau-psi));};
    double a(double X,double t){return mu*(f(t)-X) + fdash(t);};
    double a_dash(double X,double t){return -mu;};
    double alpha(double X,double t){return 2.*sigma*sigma*X - a(X,t);};
    double source(double X,double t){return sigma*sigma + mu;};
};

}// namespace PJLib

#endif
