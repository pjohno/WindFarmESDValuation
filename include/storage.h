#ifndef _storage_h_included_
#define _storage_h_included_

/* Set up a class for the option value */

#include "option.h"
#include "pj_data.h"
#include "param_vecs.h"
#include "param_maps.h"
#include "nodes.h"
#include "fixed_grid.h"
#include "finite_diff_methods.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

class Storage
{
  
protected:
  
  // store option name
  string option_name;
  // print progress to screen
  bool show_progress;
  // globals
  static const int Width=20,Iter_max=100000;
  static constexpr double Pi=3.141592653589793238462643383279502884197l;
  
  // size of arrays and important constants
  int n,m,kmax,n_zero_1,n_zero_2,time_steps;
  
  // change option types --> redundant here
  int option_type;
  // fixed frid space
  Fixed_Grid X,Y,Q;
  // option values
  Option3 option_new,option_old;
  // payoff
  Payoff3 payoff;
  // numerical methods
  Grid GG;
  Methods ADI_time_march;
  // printing
  Print3 file_print;
  
  // setup stencil constants
  virtual void setup_stencils(double tau);
  // scheme setup
  virtual void setup_scheme(void);
  // setup node values for half steps
  void adjust_scheme_X(void);
  void adjust_scheme_Y(void);
  // smooth initial condition with interpolation
  void smooth_initial_conditions(void);
  // default parameters for initial setup
  void generate_default_params(void);
  // link parameter class pointers to storage
  void assign_pointers(Param_vecs *my_parameters);
  // link parameter class pointers to storage
  void assign_pointers(ParameterMaps *my_parameters);
  // 
  void default_parameters(void);
  void default_parameters(Param_vecs *my_parameters);
  void default_parameters(ParameterMaps *my_parameters);
  // set global time
  void set_time(double tau);
  // main solver
  void main_solver(void);
  // reset initial conditions, either to zero or passing in some values
  void reset_initial_conditions(vector<vector<vector<double> > > *u_old,vector<vector<vector<double> > > *u_new);
  void reset_initial_conditions(double ***u_old,double ***u_new);
  void reset_initial_conditions(double ***U);
  void reset_initial_conditions(void);
  void reset_initial_conditions(Option3& U);
  // forcing functions
  double average_X_t_dash(double tau),average_Y_t_dash(double tau);
  
protected:
  
  // adjust the number of timesteps as a ratio to 
  int stability_factor;
  // parameters
  double dX,X_max,C_commit,dY,Y_max,dQ,Q_max,dt,interest_rate,mean_reversion_X,mean_reversion_Y
  ,average_X,average_Y,sigma_X,sigma_Y,kappa,lambda_C,lambda_D,X_C,X_D,contract_rate,
  tolerance,omega,theta,rho,maturity,penalty_charge,power_max,power_cut_out,lambda_cut,
  power_turn,contract_length,phase_X,phase_Y,amp_X,amp_Y,freq_X,freq_Y,period_X,period_Y,
  global_time,maximum_commit,minimum_commit,buy_sell_spread;
  string initial_results,final_results,plot_file;
  // setup the grid
  virtual void assign_grid(void);
  // first call to initial conditions 
  void assign_initial_conditions(void);
  // find the point at which f(x) - c switches sign
  double root_finder(double cc,double a,double b);
  // charging functions
  double L_C(double xx,double qq),L_D(double xx,double qq);
  // allow for different penalty functions
  virtual double penalties(double xx,double yy,double qq);
  // allow for different payment mechanisms
  virtual double payments(double yy);
  
public:
  
  virtual ~Storage(){};
  // constructors
  Storage();
  Storage(Param_vecs *my_parameters,bool progress_bar);
  Storage(ParameterMaps *my_parameters,bool progress_bar);
  
  // setup the grid and stuff
  void test(){assign_grid();assign_initial_conditions();};
  // test initial conditions
  void test_initial_conditions(void);
  // link parameter class to the storage parameters
  void generate_pointers(Param_vecs *my_parameters);
  void generate_pointers(ParameterMaps *my_parameters);
  // run the solver
  void runSolver(void);
  void runSolver(double ***U,double tau);
  void runSolver(double ***u_old,double ***u_new,double tau);
  void runSolver(vector<vector<vector<double> > > *u_old,vector<vector<vector<double> > > *u_new,double tau);
  void runSolver(Option3& U,double tau);
  
  void return_vector_u(int j,int k,vector<double> *u);
  void return_vector_uold(int j,int k,vector<double> *u);
  
  // change the number of timesteps via the stability factor
  void set_stability_factor(int i);
  // return the option values
  Option3& return_option_new(){return option_new;};
  Option3& return_option_old(){return option_old;};
  // return grids
  Fixed_Grid* return_grid_x(){return &X;};
  Fixed_Grid* return_grid_y(){return &Y;};
  Fixed_Grid* return_grid_q(){return &Q;};
  // return SDE functions
  double average_X_t(double tau),average_Y_t(double tau);
  // print power functions
  void print_initial_conditions(void);
  
  // return charging functions -- not used in the code
  double L_C(double xx,double qq,double cc),L_D(double xx,double qq,double cc),power_output(double xx,double qq,double cc),charge_rate(double xx,double qq,double cc);
  // different wind farm power functions
  virtual double Power_Curve(double xx);
  // allow for new power stack functions to be applied
  virtual double Power_Stack(double yy);
  
  void print_results(string filename);
  void print_results(string filename,int j);
  void print_results(string filename,int j,double param,ofstream *graph);
  void print_results_old(string filename);
	
};

#endif
