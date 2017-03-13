#ifndef _pj-option_h_included_
#define _pj-option_h_included_
// #define DEBUG
/* Set up a class for the option value */


#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "fixed_grid.h"
#include "param_vecs.h"
#include "param_maps.h"

namespace PJLib {

class Payoff
{
  
  protected:
	double strike;
	
  public:
	
	Payoff(){strike =0;};
	Payoff(double X){change_strike(X);};
	void change_strike(double X){strike=X;};
	virtual double operator()(double x){return std::max(strike - x,0.);};
	
};

class Interpolation
{
  public: double interpolate(double x,std::vector<double> y,std::vector<double> g);
};

class Option: public Interpolation
{
  
  protected:
	// bool U setup
	bool u_setup;
	// number of elements in the option
	int no_of_elements,accuracy_x;
	// elements stored in U
	double *U;
	double dx,x_min,x_max;
	// check compatability of options
	bool check(Option& a,Option& b){return (a.size_x()==b.size_x());};
	
  public:
	// default deconstructor
	~Option()
	{
	  delete [] U;
	};
	// default constructor, with one commit
	Option(){default_values();setup(no_of_elements,x_min,x_max);};
	Option(int n,double x_min_,double x_max_){default_values();setup(n,x_min_,x_max_);};
	Option(Fixed_Grid *x){default_values();setup(x);assign_value();};
	// supply initial conditions
	Option(Fixed_Grid *x,Payoff *vanilla){default_values();setup(x);assign_value(vanilla);};
	// default values
	void default_values(void){accuracy_x = 4;U = NULL;u_setup = false;no_of_elements =  4;x_min = 0.;x_max = 1.;};
	// link option values to the grid
	void setup(Fixed_Grid *x)
	{
	  setup(x->size(),x->min(),x->max());
	};
	void setup(int n,double x_min_,double x_max_){
	  x_min = x_min_;
	  x_max = x_max_;
	  dx = (x_max - x_min)/double(n-1);
	  update(n);
	};
	// default option value to zero
	void assign_value(){for(int i=0;i<no_of_elements;i++){U[i]=0.;};};
	// assign option value a payoff function;
	void assign_value(Payoff *vanilla){for(int i=0;i<no_of_elements;i++){U[i]=(*vanilla)(x_min+i*dx);};};
	// update 
	void update(int n_);
	// change accuracy
	void change_accuracy_x(int temp){accuracy_x=std::max(1,temp);};
	// interpolate value
	double operator()(double x);
	// return value
	double& operator[](int index){return U[index];};
	Option& operator=(double *rhs){for(int i=0;i<no_of_elements;i++)U[i]=rhs[i];return *this;};
	Option& operator=(double rhs){for(int i=0;i<no_of_elements;i++)U[i]=rhs;return *this;};
	Option& operator+=(double rhs){for(int i=0;i<no_of_elements;i++)U[i]+=rhs;return *this;};
	Option& operator=(Option& rhs){
	  #ifdef DEBUG
	  if(!check(*this,rhs)){std::cout << " incompatible option sizes...";return *this;}
	  #endif // DEBUG
	  for(int i=0;i<no_of_elements;i++)U[i]=rhs[i];return *this;
	};
	
	// return size of array
	int size_x(void){return no_of_elements;};
	double last(void){return U[no_of_elements-1];};
	// friend functions
	Option& max(Option& a,Option& b);
	inline Option& max(Option& a){for(int i=0;i<no_of_elements;i++)U[i]=std::max(a[i],U[i]);return *this;};
	
};

class Print
{
  
  protected:
	
	// parameters for printing
	int n,precision;
	double dx,x_min,x_max;
	// output file
	std::string filename;
	std::ofstream file;
	
  public:
	
	Print(){precision=8;setup_x(11,0.,1.);change_filename(std::string("default"));};
	Print(Fixed_Grid *x,std::string _filename){precision=8;setup_x(x); change_filename(_filename);	};
	void setup_x(int _n,double min,double max){x_points(_n);x_range(min,max);};
	void setup_x(Fixed_Grid *x){x_points(x->size());x_range(x->min(),x->max());};
	void change_precision(int temp){precision=temp;};
	void change_filename(std::string _filename){filename=_filename+std::string(".dat");};
	virtual void calc_grid(){dx = (x_max - x_min)/double(n-1);if(n<2)dx=0.;};
	void x_points(int _n){n=_n;};
	void x_range(double min,double max){x_min=min;x_max=max;};
	void print(Option &u);
	void link_printer_config(Param_vecs *print_config);
	void link_printer_config(ParameterMaps *print_config);
	
};


}// PJ namespace

#endif
