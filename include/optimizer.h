#ifndef _optimizer_h_included_
#define _optimizer_h_included_

/* Set up a class for the option value */
#include <param_vecs.h>
#include <storage.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

class Optimizer
{
	
	protected:

	  // no of commits possible
	  int no_of_commits;
	  // grid variables (for printing)
	  int grid_i,grid_j,grid_k;
	  double dx,dy,dz;
	  // check on whether vectors etc are setup
	  bool option_vector_isnt_setup,we_have_option_max,change_commits,run_starter;
	  // storage for individual commit options
	  std::vector<Storage*> total_option_value;
	  // parameter linking to the commit options
	  std::vector<Param_vecs*> all_parameters;
	  // option value, and commits
	  std::vector<std::vector<std::vector<double> > > *option_max_old,*option_max_new,*option_commit,*option_temp;
	  // change the number of commits
	  void resize_total_option(int no_of_options);
	  // interpolate option value
	  double option_value_1(double value_stock_1,int j,int k,int accuracy);
	  double option_value_2(double value_stock_1,double value_stock_2,int k,int accuracy);
	
	public:
		
	  // default deconstructor
	  ~Optimizer();
	  // default constructor, with one commit
	  Optimizer();
	  // constructor, specifying number of commits
	  Optimizer(int no_of_options,bool progress_bar);
	  // constructor, using no price variable (returning 1)
	  Optimizer(int no_of_options,bool progress_bar,int no_price);
	  // force charging during certain time periods
	  int force_commit(double tau);
	  // interpolate option value
	  double option_value_3(double value_stock_1,double value_stock_2,double value_stock_3,int accuracy);
	  
	  double rms_error(void),option_average(void);
	  
	  double estimate_perpetual_option(double time,double value);
	  
	  void test_sin_curves(void);
	  
	  void set_perpetual_option_conditions(double value);
	  
	  void adjust_perpetual_option_conditions(double value);
	  
	  virtual void perpetual_pricer(bool print_results);
	  
	  virtual void run_perpetual(void);
	  
	  void set_base_parameters(int i);
	  
	  void setup_grid(void);
	  
	  void test(void);
	  
	  void vary_parameter(int parameter_ref,double low,double high,int its,string filename);
	  
	  void grid_test(void);

	  void grid_test_dt(string filename,int grid_start,int grid_finish,int grid_incr,int param);
	  
	  void grid_test_X_max(void);
	  
	  void grid_test_X_max_per(void);
	  
	  void solve_one_option(void);
	  
	  void solve_for_multiple_option(int no_of_hours,string filename);
	  
	  virtual void solve_multiple_time(int no_of_hours,string filename,bool display,bool reset);
	  
	  void solve_for_all_C(double tau);
	  
	  void find_option_max(int current_commits);
	  // increase the number of timesteps via stability factor
	  void increase_timesteps(int i);
	  
	  void print_results(string filename,bool debug);
	  
	  void vary_two_parameters(int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,string filename);
	  
	  void vary_two_parameters_multi(int no_of_contracts,int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,string filename);
	  
	  void vary_tp_multi_new(int no_of_contracts,int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,string filename);
	  
	  void vary_perpetual(int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,string filename);
};

#endif
