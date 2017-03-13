#ifndef _delay_ops_test_h_included_
#define _delay_ops_test_h_included_

/* Set up a class for the option value */

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "storage.h"

class Param_vecs2;
class Option3;
class Fixed_Grid;

class Delay_Op2
{
  // do all linking
  void all_linking();
  // link up print configuration
  void link_print_config();
  // initial conditions set to zero
  bool first_run;  
  // default file names
  static const std::string str_default_parameters,str_default_store,str_default_wind,str_default_price;
  static const std::string str_default_trading,str_default_grid,data_dir;
  static const std::string str_default_option3,str_default_option4;
  
  protected:
	
	// default directories
	static const std::string config_dir;
	// initial setup function
	void input_values(std::string base,std::string store,std::string wind,std::string price,std::string trading,std::string grid_size);
	// initial setup function
	void input_values(std::string base);
	// no of commits possible
	int no_of_commits;
	// markets have a delay
	int market_delay;
	// debug messages
	int debug;
	// grid variables
	Fixed_Grid X,Y,Q,C;
	// storage for individual commit option
	Storage *option_pricer;
	// parameter linking to the commit option
	ParameterMaps all_parameters,print_config,print_config4,store_config,wind_config,trading_config,grid_config,price_config;
	// single option max
	Option3 option_max,option_max_old;
	// option value, and commits
	Option4 option_commits,option_commits_new,option_commits_old;
	// printing 
	Print3 file_print;
	Print4 file_print4;
	// update storage of grids and option values
	void update_storage();
	// base parameters
	void default_values();
	// base print configs
	void default_print_values();
	// estimate perpetual coupon option
	double estimate_perpetual_option(double time,double value);
	// find rms error between two option values
	double rms_error();
	// solve for multiple times
	void solve(int no_of_hours,std::string filename,bool display,bool reset);
	
  public:
	// default deconstructor
	~Delay_Op2();
	// default constructor, with one commit
	Delay_Op2();
	// use configuration files
	Delay_Op2(std::string base,std::string store,std::string wind,std::string price,std::string trading,std::string grid_size);
	// base parameters from file
	void default_values(std::string filename);
	// print configs from file
	void default_print_values(std::string filename,std::string filename4);
	// solve for multiple times
	void solve(int no_of_hours,std::string filename);
	// perpetual pricer
	// Try and extrapolate to find the perpetual option value quickly
	double perpetual_pricer(bool print_results,bool run_starter);
	// print all parameters to screen
	void print_params();
	// return annual returns
	double annual_returns();
	
};


#endif
