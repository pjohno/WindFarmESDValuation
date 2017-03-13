/* Set up a class for the option value */
#ifndef _finite_diff_method_h_included_
#define _finite_diff_method_h_included_

#include "fixed_grid.h"
#include "nodes.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

class Methods
{
  private:
    
  protected:
	
	bool x_setup,y_setup;
	int n,m,p;
	double *a_X,*b_X,*c_X,*d_X,*w_X;
	double *a_Y,*b_Y,*c_Y,*d_Y,*w_Y;
	
	// use default setup
	void default_values();
	// setup  default storage
	void setup_storage(int n_,int m_,int p_);
	// resize storage
	void resize_X(int n_),resize_Y(int m_);
	// clear storage
	void clear_X(),clear_Y();
	// assign storage
	void assign_X(int n_),assign_Y(int m_);
	
  public:
    
	// deconstructor
	~Methods(){clear_X();clear_Y();}
	// constructor
	Methods();
	// constructor with initial setup for storage
	Methods(Fixed_Grid &X,Fixed_Grid &Y,Fixed_Grid &Z);
	// update vector sizes
	void update(Fixed_Grid &X,Fixed_Grid &Y,Fixed_Grid &Z);
	// generic tridiagonal solver for X grid
	void tridag_solver_X();
	// generic tridiagonal solver for Y grid
	void tridag_solver_Y();
	// generic tridiagonal solver with vectors	
	void tridag_solver(vector<double> *w,vector<double> aa,vector<double>bb,vector<double>cc,	vector<double> dd);
	// alternate difference implicit method, solving across X -> use X grid
	void ADI_imp_X(double dt,double ***w, Grid *GG);
	// alternate difference implicit method, solving across Y -> use Y grid
	void ADI_imp_Y(double dt,double ***w, Grid *GG);
	// generic interpolation scheme
	double interpolate(double x,vector<double> y,vector<double> g);
	// polynomial setup
	void polynomial(vector<double> y,vector<double> g,vector<double> *coefficients);
	// polynomial evaluation
	double poly_eval(double x,vector<double> coefficients);
	// find maximum value of the polynimial
	double poly_max(double x_min,double x_max,vector<double> coefficients);
	// perform richardson extrapolation
	double richardson_extrap(double y_1,double y_2,int n_1,int n_2,int p);
		
};

#endif
