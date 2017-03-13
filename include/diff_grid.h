#ifndef _diff_grid_h_included_
#define _diff_grid_h_included_

#include "fixed_grid.h"
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Stencil
{
  
  public:		
	int adjust;
	double aa,bb,cc;
		
  public:
		
	void assign_adjust(int adjust_);
	int adjuster(void);
	int operator()(){return adjust;};
	void setup_stencil(double aa_,double bb_,double cc_,int adjust_);
	inline void operator()(double aa_,double bb_,double cc_,int adjust_){setup_stencil(aa_,bb_,cc_,adjust_);};
	void return_stencil(double *aa_,double *bb_,double *cc_,int *adjust_);
	inline void operator()(double *aa_,double *bb_,double *cc_,int *adjust_){return_stencil(aa_,bb_,cc_,adjust_);};
	void print_stencil(void);
	
	Stencil();
	Stencil& operator=(double a){aa=a;bb=a;cc=a;adjust=0;return *this;};
	
};

class Node
{
  private:
	
	bool node_setup;
	
  protected:
	
	double g,Z;
	
  public:
	
	void set_Z(double Z_){Z=Z_;};
	void set_g(double g_){g=g_;};

	double ZZ(){return Z;};
	double operator()(){return g;};
	inline Node& operator=(double a){Z=a;g=a;return *this;};
	inline Node& operator+=(double a){Z=g+a;return *this;};
	inline double operator+(double a){return g+a;};
	inline double operator-(double a){return g-a;};
	
};


#endif
