#pragma once

/* Set up a class for the option value */
#include <generic_SDE.h>
#include <fixed_grid.h>
#include <option.h>
#include <optimizer.h>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

class Simulator: public Optimizer
{
  
  private:
	bool run_starter;
	// allows checking on whether vectors and pointers are setup
	bool commit_strategy_setup,generate_distributions;
	// phase space
	double X,Y,Q,C,T;
	// number of points in simulation
	int no_of_points,no_of_runs;
	// timestep size
	double dt;
	// commit strategy in (X;Y;Q;hour)
	std::vector<std::vector<std::vector<std::vector<double>>>> commit_strategy;
	// monte distributions for a given day
	Generic_SDE *monte_X,*monte_Y;
	
  public:
	
	// deconstructor
	~Simulator();
	// constructor
	Simulator();
	Simulator(int n,int m,int p,int q,int no_of_options,bool progress_bar,int no_price);
	// overwrite solver function to put strategy into commit_strategy
	void solve_multiple_time(int no_of_hours,std::string filename,bool display,bool reset);
	// override perpetual pricer to not call base parameters
	void run_perpetual(void);
	// override parameter iterator so that base parameters maintained
	void vary_perpetual(int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,std::string filename);
	// write strategy to vector
	void write_strategy_to_vector(int hour);
	// generate a days battery usage
	void generate_simulated_day(double store_level,std::string filename);
	// output from the wind farm + battery
	double output(void);
	// Average expected value of X
	double average_x(double tau);
	// Average expected value of Y
	double average_y(double tau);
	// set up a smoothing store
	void smoothing_store(void);
	// set up an arbitraging store
	void arbitrage_store(void);
	// set up with no store
	void no_store(void);
	// set test case
	void test_case(void);
	// work out implied hours of power per year
	double hours_per_year(void);
	// allow access to changing parameters outside the class
	void changeParameter(double x,int parameter_ref);
	// reset initial conditions -- do this unless parameter is non-invasive
	void reset_initial_conditions(void){run_starter=true;};
	
};

