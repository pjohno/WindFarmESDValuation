#ifndef _option_h_included_
#define _option_h_included_
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

typedef PJLib::Param_vecs Param_vecs;
typedef PJLib::ParameterMaps ParameterMaps;

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

class Payoff2
{
  
  protected:
	double strike_x,strike_y;
	
  public:
	
	Payoff2(){strike_x=0;strike_y=0;};
	Payoff2(double X,double Y){change_strike(X,Y);};
	void change_strike(double X,double Y){strike_x=X;strike_y=Y;};
	virtual double operator()(double x,double y){return std::max(std::max(strike_y - y,strike_x - x),0.);};
	
};

class Payoff3
{
  
  protected:
	double strike_x,strike_y,strike_z;
	
  public:
	
	Payoff3(){strike_x=0;strike_y=0;strike_z=0;};
	Payoff3(double X,double Y,double Z){change_strike(X,Y,Z);};
	void change_strike(double X,double Y,double Z){strike_x=X;strike_y=Y;strike_z=Z;};
	virtual double operator()(double x,double y,double z){return std::max(std::max(strike_z - z,std::max(strike_y - y,strike_x - x)),0.);};
	
};

class Payoff4
{
  
  protected:
	double strike_x,strike_y,strike_z,strike_w;
	
  public:
	
	Payoff4(){strike_x=0;strike_y=0;strike_z=0;strike_w=0;};
	Payoff4(double X,double Y,double Z,double W){change_strike(X,Y,Z,W);};
	void change_strike(double X,double Y,double Z,double W){strike_x=X;strike_y=Y;strike_z=Z;strike_w=W;};
	virtual double operator()(double x,double y,double z,double w){return std::max(std::max(strike_w - w,std::max(strike_z - z,std::max(strike_y - y,strike_x - x))),0.);};
	
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
	Option(int n,double x_min_,double x_max_){default_values();setup(n,x_min_,x_max_);update(n);};
	Option(Fixed_Grid *x){default_values();setup(x);assign_value();};
	// supply initial conditions
	Option(Fixed_Grid *x,Payoff *vanilla){default_values();setup(x);assign_value(vanilla);};
	// default values
	void default_values(void){accuracy_x = 4;U = NULL;u_setup = false;no_of_elements =  4;x_min = 0.;x_max = 1.;};
	// link option values to the grid
	void setup(Fixed_Grid *x)
	{
	  setup(x->size(),x->min(),x->max());
	  update(x->size());
	};
	void setup(int n,double x_min_,double x_max_){
	  x_min = x_min_;
	  x_max = x_max_;
	  dx = (x_max - x_min)/double(n-1);
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


class Option2: public Interpolation
{
  
  protected:
	// bool U setup
	bool u_setup;
	// number of elements in the option
	int n,m,accuracy_x,accuracy_y;
	// elements stored in U
	double **U;
	double dx,x_min,x_max;
	double dy,y_min,y_max;
	// check compatibility of options
	bool check(Option2& a,Option2& b){return ((a.size_x()==b.size_x()) and (a.size_y()==b.size_y()));};
	
  public:
	// default deconstructor
	~Option2()
	{
	  if(u_setup){
		for(int i=0;i<n;i++)delete [] U[i];
		delete [] U;
	  }
	};
	// default constructor, with one commit
	Option2(){default_values();};
	Option2(Fixed_Grid *x,Fixed_Grid *y);
	Option2(Fixed_Grid *x,Fixed_Grid *y,Payoff2 *vanilla);
	void default_values(void);
	// link option values to the grid
	void setup_x(int n_,double x_min_,double x_max_){
	  x_min = x_min_;
	  x_max = x_max_;
	  dx = (x_max - x_min)/double(n_-1);
	};
	void setup_y(int n_,double x_min_,double x_max_){
	  y_min = x_min_;
	  y_max = x_max_;
	  dy = (y_max - y_min)/double(n_-1);
  	};
	
	// default option value to zero
	void assign_value(){for(int i=0;i<n;i++){for(int j=0;j<m;j++){U[i][j]=0.;};};};
	// assign option value a payoff function;
	void assign_value(Payoff2 *vanilla);
	// reform option value into new grid
	void update(Fixed_Grid *x,Fixed_Grid *y);
	// change accuracy
	void change_accuracy_x(int temp){accuracy_x=std::min(accuracy_x,temp);};
	void change_accuracy_y(int temp){accuracy_y=std::min(accuracy_y,temp);};
	// clear out option elements
	void clear()
	{
	  for(int i=0;i<n;i++)delete [] U[i];
	  delete [] U;
	};
	// interpolate value
	double operator()(double x,double y);
	double operator()(int i,double y);
	// return value
	inline double* operator[](int index){return U[index];};
	// copy value from pointer -- very risky!!!
	Option2& operator=(double **rhs){for(int i=0;i<n;i++)for(int j=0;j<m;j++)U[i][j]=rhs[i][j];return *this;};
	// set option equal to constant
	Option2& operator=(double rhs){for(int i=0;i<n;i++)for(int j=0;j<m;j++)U[i][j]=rhs;return *this;};
	Option2& operator+=(double rhs){for(int i=0;i<n;i++)for(int j=0;j<m;j++)U[i][j]+=rhs;return *this;};
	Option2& operator=(Option2& rhs){
	  #ifdef DEBUG
	  if(!check(*this,rhs)){std::cout << " incompatible option sizes...";return *this;}
	  #endif // DEBUG
	  for(int i=0;i<n;i++)for(int j=0;j<m;j++)U[i][j]=rhs[i][j];return *this;
	};
	// return size of array
	int size_x(void){return n;};
	int size_y(void){return m;};
	Option2& max(Option2& a,Option2& b);
	inline Option2& max(Option2& a){for(int i=0;i<n;i++)for(int j=0;j<m;j++)U[i][j]=std::max(a[i][j],U[i][j]);return *this;}
	
};

class Print2:public Print
{
  
  protected:
	
	// parameters for printing
	int m;
	double dy,y_min,y_max;
	
  public:
	
	Print2(){setup_y(11,0.,1.);};
	Print2(Fixed_Grid *x,Fixed_Grid *y,std::string _filename):Print(x,_filename){setup_y(y);}
	void setup_y(Fixed_Grid *y){y_points(y->size());y_range(y->min(),y->max());}
	void setup_y(int _m,double min,double max){y_points(_m);y_range(min,max);}
	virtual void calc_grid();
	void y_points(int _m){m=_m;}
	void y_range(double min,double max){y_min=min;y_max=max;}
	void print(Option2 &u);
	void print_y_x(Option2 &u);
	void print_U_and_V(Option2 &u,Option2 &v);
	void print_y_x_U_and_V(Option2 &u,Option2 &v);
	void link_printer_config2(Param_vecs *print_config);
	void link_printer_config2(ParameterMaps *print_config);
	
};

class Option3: public Interpolation
{

  private:
	
  protected:
	// bool U setup
	bool u_setup;
	// elements stored in U
	double ***U;
	// number of elements in the option
	int n,m,p,accuracy_x,accuracy_y,accuracy_z;
	double dx,x_min,x_max;
	double dy,y_min,y_max;
	double dz,z_min,z_max;
	// check compatibility of options
	bool check(Option3& a,Option3& b){
	  #ifdef DEBUG
	  std::cout << " sizes :: " << a.size_x() <<" " << b.size_x() << " " << a.size_y() << " " << b.size_y() << " " << a.size_z() << " " << b.size_z();
	  #endif
	  return ((a.size_x()==b.size_x()) and (a.size_y()==b.size_y()) and (a.size_z()==b.size_z()));
	};
		
  public:
	// default deconstructor
	~Option3()
	{
	  #ifdef DEBUG
	  std::cout << " clear option...";
	  #endif
	  if(u_setup){
		for(int i=0;i<n;i++){
		  for(int j=0;j<m;j++)delete [] U[i][j];
		  delete [] U[i];
		}
		delete [] U;
	  }
	};
	// default constructor, with one commit
	Option3(){default_values();};
	Option3(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z);
	Option3(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,Payoff3 *vanilla);
	void default_values(void);
	// link option values to the grid
	void setup_x(int n_,double x_min_,double x_max_);
	void setup_y(int n_,double x_min_,double x_max_);
	void setup_z(int n_,double x_min_,double x_max_);
	// default option value to zero
	void assign_value(){for(int i=0;i<n;i++){for(int j=0;j<m;j++){for(int k=0;k<p;k++){U[i][j][k]=0.;};};};};
	// assign option value a payoff function;
	void assign_value(Payoff3 *vanilla);
	// reform option value into new grid
	void update(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z);
	// change accuracy
	void change_accuracy_x(int temp){accuracy_x=std::min(accuracy_x,temp);};
	void change_accuracy_y(int temp){accuracy_y=std::min(accuracy_y,temp);};
	void change_accuracy_z(int temp){accuracy_z=std::min(accuracy_z,temp);};
	// clear out option elements
	void clear()
	{
	  for(int i=0;i<n;i++){
		for(int j=0;j<m;j++)delete [] U[i][j];
		delete [] U[i];
	  }
	  delete [] U;
	};
	// interpolate value
	double operator()(double x,double y,double z);
	double operator()(int i,double y,double z);
	double operator()(int i,int j,double z);
	// return value
	inline double** operator[](int index){return U[index];};
	inline double*** operator()(){return U;};
	Option3& intrplt(Option3& rhs){
	  double x,y,z;
/*	  rhs.change_accuracy_x(2);
	  rhs.change_accuracy_y(2);
	  rhs.change_accuracy_z(2);*/
	  for(int i=0;i<n;i++)
	  {
		x = x_min + i*dx;
		for(int j=0;j<m;j++)
		{
		  y = y_min + j*dy;
		  for(int k=0;k<p;k++)
		  {
			z = z_min + k*dz;
			U[i][j][k]=rhs(x,y,z);
		  }
		}
	  }
/*	  rhs.change_accuracy_x(4);
	  rhs.change_accuracy_y(4);
	  rhs.change_accuracy_z(4);*/
	  return *this;
	};
	  // set option value equal to another option value
	Option3& operator=(Option3& rhs){
	  #ifdef DEBUG
	  if(!check(*this,rhs)){std::cout << " incompatible option sizes...";return *this;}
	  std::cout << " option copy";
	  #endif // DEBUG
	  for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)U[i][j][k]=rhs[i][j][k];return *this;};
	Option3& operator+=(Option3& rhs){
	  #ifdef DEBUG
	  if(!check(*this,rhs)){std::cout << " incompatible option sizes...";return *this;}
	  std::cout << " option copy";
	  #endif // DEBUG
	  for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)U[i][j][k]=std::max(U[i][j][k],rhs[i][j][k]);return *this;};
	Option3& operator=(double*** rhs){for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)U[i][j][k]=rhs[i][j][k];return *this;};
	Option3& operator=(double constant){for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)U[i][j][k]=constant;return *this;};
	Option3& operator+=(double constant){for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)U[i][j][k]+=constant;return *this;};
	// return size of array
	int size_x(void){return n;};
	int size_y(void){return m;};
	int size_z(void){return p;};
	// maximum of two options
	Option3& max(Option3& a,Option3& b);
	inline Option3& max(Option3& a){for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)U[i][j][k]=std::max(U[i][j][k],a[i][j][k]);return *this;};
	
};

double rms(Option3& a,Option3& b);

class Print3:public Print2
{
  
  protected:
	
	// parameters for printing
	int p;
	double dz,z_min,z_max;
	
  public:
	
	Print3(){setup_z(11,0.,1.);}
	Print3(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,std::string _filename):Print2(x,y,_filename){setup_z(z);}
	void setup_z(int _p,double min,double max){z_points(_p);z_range(min,max);}
	void setup_z(Fixed_Grid *z){z_points(z->size());z_range(z->min(),z->max());}
	virtual void calc_grid();
	void z_points(int _p){p=_p;}
	void z_range(double min,double max){z_min=min;z_max=max;}
	void print(Option3 &u);
	void print_x_y(Option3 &u);
	void print_x_z(Option3 &u);
	void link_printer_config3(Param_vecs *print_config);
	void link_printer_config3(ParameterMaps *print_config);
	
};

class Option4: public Interpolation
{

  protected:
	// bool U setup
	bool u_setup;
	// elements stored in U
	double ****U;
	// number of elements in the option
	int n,m,p,q,accuracy_x,accuracy_y,accuracy_z,accuracy_w;
	double dx,x_min,x_max;
	double dy,y_min,y_max;
	double dz,z_min,z_max;
	double dw,w_min,w_max;
	// check compatibility of options
	bool check(Option4& a,Option4& b){return ((a.size_x()==b.size_x()) and (a.size_y()==b.size_y()) and (a.size_z()==b.size_z()) and (a.size_w()==b.size_w()));};
	
  public:
	// default deconstructor
	~Option4()
	{
	  #ifdef DEBUG
	  std::cout << " clear option...";
	  #endif
	  if(u_setup){
		for(int i=0;i<n;i++){
		  for(int j=0;j<m;j++){
			for(int k=0;k<p;k++)delete [] U[i][j][k];
			delete [] U[i][j];
		  }
		  delete [] U[i];
		}
		delete [] U;
	  }
	};
	// default constructor, with one commit
	Option4(){default_values();};
	Option4(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,Fixed_Grid *w);
	Option4(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,Fixed_Grid *w,Payoff4 *vanilla);
	void default_values(void);
	// link option values to the grid
	void setup_x(int size,double min,double max);
	void setup_y(int size,double min,double max);
	void setup_z(int size,double min,double max);
	void setup_w(int size,double min,double max);
	// default option value to zero
	void assign_value(){for(int i=0;i<n;i++){for(int j=0;j<m;j++){for(int k=0;k<p;k++){for(int l=0;l<q;l++){U[i][j][k][l]=0.;}}}}};
	// assign option value a payoff function;
	void assign_value(Payoff4 *vanilla);
	// reform option value into new grid
	void update(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,Fixed_Grid *w);
	// change accuracy
	void change_accuracy_x(int temp){accuracy_x=std::min(accuracy_x,temp);};
	void change_accuracy_y(int temp){accuracy_y=std::min(accuracy_y,temp);};
	void change_accuracy_z(int temp){accuracy_z=std::min(accuracy_z,temp);};
	void change_accuracy_w(int temp){accuracy_w=std::min(accuracy_w,temp);};
	// clear out option elements
	void clear()
	{
	  for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
		  for(int k=0;k<p;k++)delete [] U[i][j][k];
		  delete [] U[i][j];
		}
		delete [] U[i];
	  }
	  delete [] U;
	};
	// interpolate value
	double operator()(double x,double y,double z,double w);
	double operator()(int i,double y,double z,double w);
	double operator()(int i,int j,double z,double w);
	double operator()(int i,int j,int k,double w);
	// return value
	inline double*** operator[](int index){return U[index];};
	// copy option value
	Option4& operator=(Option4& rhs){
	  #ifdef DEBUG
	  if(!check(*this,rhs)){std::cout << " incompatible option sizes...";return *this;}
	  #endif // DEBUG
	  for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)for(int l=0;l<q;l++)U[i][j][k][l]=rhs[i][j][k][l];return *this;
	};
	// set option value equal to a constant
	Option4& operator=(double rhs){for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)for(int l=0;l<q;l++)U[i][j][k][l]=rhs;return *this;};
	Option4& operator+=(double rhs){for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)for(int l=0;l<q;l++)U[i][j][k][l]+=rhs;return *this;};
	// copy option value
	void copy(int i,Option3& rhs){
	  for(int j=0;j<m;j++)for(int k=0;k<p;k++)for(int l=0;l<q;l++)U[i][j][k][l]=rhs[j][k][l];
	};
	// interpolate values
	Option4& intrplt(Option4& rhs){
	  double x,y,z,w;
/*	  rhs.change_accuracy_x(2);
	  rhs.change_accuracy_y(2);
	  rhs.change_accuracy_z(2);
	  rhs.change_accuracy_w(2);*/
	  for(int i=0;i<n;i++)
	  {
		x = x_min + i*dx;
		for(int j=0;j<m;j++)
		{
		  y = y_min + j*dy;
		  for(int k=0;k<p;k++)
		  {
			z = z_min + k*dz;
			for(int l=0;l<q;l++)
			{
			  w = w_min + l*dw;
			  U[i][j][k][l]=rhs(x,y,z,w);
			}
		  }
		}
	  }	
/*	  rhs.change_accuracy_x(4);
	  rhs.change_accuracy_y(4);
	  rhs.change_accuracy_z(4);
	  rhs.change_accuracy_w(4);*/
	  return *this;
	};
	// return size of array
	int size_x(void){return n;};
	int size_y(void){return m;};
	int size_z(void){return p;};
	int size_w(void){return q;};
	
};

double rms(Option4& a,Option4& b);

class Print4:public Print3
{
  
  protected:
	
	// parameters for printing
	int q;
	double dw,w_min,w_max;
	
  public:
	
	Print4(){setup_w(11,0.,1.);};
	Print4(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,Fixed_Grid *w,std::string _filename):Print3(x,y,z,_filename){setup_w(w);}
	void setup_w(int size,double min,double max){w_points(size);w_range(min,max);}
	void setup_w(Fixed_Grid *w){w_points(w->size());w_range(w->min(),w->max());}
	virtual void calc_grid();
	void w_points(int size){q=size;}
	void w_range(double min,double max){w_min=min;w_max=max;}
	void print(Option4 &u);
	void print_y_z(Option4 &u);
	void print_x_y(Option4 &u);
	void link_printer_config4(Param_vecs *print_config);
	void link_printer_config4(ParameterMaps *print_config);
	
};


#endif
