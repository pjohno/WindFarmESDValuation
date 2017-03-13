#pragma once

/* Set up a class for the option value */

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "storage.h"
#include "param_vecs.h"
using namespace PJLib;

class Option3;
class Fixed_Grid;

class Delay_Op
{
  // do all linking
  void all_linking(){link_print_config();link_store_config();link_wind_config();link_trading_config();link_price_config();link_grid_config();};
  // link up trading configuration
  void link_trading_config();
  // link up print configuration
  void link_print_config();
  // link up store configuration
  void link_store_config();
  // link up wind configuration
  void link_wind_config();
  // link up price configuration
  void link_price_config();
  // link up grid configuration
  void link_grid_config();
  // initial conditions set to zero
  bool first_run;  
  
  protected:
	
	// initial setup function
	void input_values(std::string base,std::string store,std::string wind,std::string price,std::string trading,std::string grid_size);
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
	Param_vecs2 all_parameters,print_config,print_config4,store_config,wind_config,trading_config,grid_config,price_config;
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
	// base trading configs
	void default_trading_values();
	// base store configs
	void default_store_values();
	// base wind configs
	void default_wind_values();
	// base grid configs
	void default_grid_values();
	// base price configs
	void default_price_values();
	// estimate perpetual coupon option
	double estimate_perpetual_option(double time,double value);
	// find rms error between two option values
	double rms_error();
	// solve for multiple times
	void solve(int no_of_hours,std::string filename,bool display,bool reset);
	
  public:
	// default deconstructor
	~Delay_Op();
	// default constructor, with one commit
	Delay_Op();
	// use configuration files
	Delay_Op(std::string base,std::string store,std::string wind,std::string price,std::string trading,std::string grid_size);
	// base parameters from file
	void default_values(std::string filename);
	// print configs from file
	void default_print_values(std::string filename,std::string filename4);
	// store configs from file
	void default_store_values(std::string filename);
	// wind configs from file
	void default_wind_values(std::string filename);
	// grid configs from file
	void default_grid_values(std::string filename);
	// price configs from file
	void default_price_values(std::string filename);
	// trade configs from file
	void default_trading_values(std::string filename);
	// solve for multiple times
	void solve(int no_of_hours,std::string filename);
	// perpetual pricer
	// Try and extrapolate to find the perpetual option value quickly
	double perpetual_pricer(bool print_results,bool run_starter);
	// print all parameters to screen
	void print_params(std::ostream &output);
	// return annual returns
	double annual_returns();
	
};

