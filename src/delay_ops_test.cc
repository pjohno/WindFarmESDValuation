#include "delay_ops_test.h"
#include "storage.h"
#include "option.h"
#include "no_price.h"
#include "finite_diff_methods.h"
#include "countdown.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#define DEBUG

const std::string Delay_Op2::config_dir="Config/";
const std::string Delay_Op2::data_dir="Data/";
const std::string Delay_Op2::str_default_parameters="default_parameters.in";
const std::string Delay_Op2::str_default_option3="print_option3.in";
const std::string Delay_Op2::str_default_option4="print_option4.in";
const std::string Delay_Op2::str_default_store="default_store.in";
const std::string Delay_Op2::str_default_wind="default_wind.in";
const std::string Delay_Op2::str_default_trading="default_trading.in";
const std::string Delay_Op2::str_default_price="default_price.in";
const std::string Delay_Op2::str_default_grid="default_grid.in";
  
Delay_Op2::~Delay_Op2()
{
  #ifdef DEBUG
  std::cout << " Delay_Op2 clear parameters...\n";std::cout.flush();
  #endif
  all_parameters.clear();
  #ifdef DEBUG
  std::cout << " delete storage option...\n";
  #endif
  delete option_pricer;
  #ifdef DEBUG
  std::cout << " clean up complete...\n";
  #endif
}  

Delay_Op2::Delay_Op2()
{
  option_pricer = new Storage(&all_parameters,false);
  all_linking();first_run=true;
  input_values(str_default_parameters);
}

Delay_Op2::Delay_Op2(std::string base,std::string store,std::string wind,std::string price,std::string trading,std::string grid_size)
{
  option_pricer = new Storage(&all_parameters,false);
  all_linking();first_run=true;
  // use 
  input_values(base,store,wind,price,trading,grid_size);
}

void Delay_Op2::all_linking()
{
  link_print_config();
  all_parameters.assignParameter(&no_of_commits,"no_of_commits"," No of commits :: ",5,1,1000);
  all_parameters.assignParameter(&market_delay,"market_delay"," Market Delay :: ",0);
  all_parameters.assignParameter(&debug,"debug"," Show debugging messages :: ",0);
}

void Delay_Op2::input_values(std::string base)
{
  #ifdef DEBUG
  std::cout << " Use the following configuration files...\n";
  std::cout << base <<" "<<str_default_option3<< ""<<str_default_option4 << "\n";
  #endif
  default_values(base);
  default_print_values(str_default_option3,str_default_option4);
}

void Delay_Op2::input_values(std::string base,std::string store,std::string wind,std::string price,std::string trading,std::string grid_size)
{
  #ifdef DEBUG
  std::cout << " Use the following configuration files...\n";
  std::cout << base << " " <<store<<" "<<wind<<" "<<price<<" "<<trading<<" "<<grid_size<<" "<<str_default_option3<< ""<<str_default_option4 << "\n";
  #endif
  default_values(base);
  default_values(store);
  default_values(wind);
  default_values(price);
  default_values(trading);
  default_values(grid_size);
  default_print_values(str_default_option3,str_default_option4);
}

void Delay_Op2::update_storage()
{
  #ifdef DEBUG
  std::cout << " Update grids...\n";
  #endif
  // update grids
  X.setup_x(all_parameters.returnInteger("n"),0.,all_parameters.returnDouble("X_max"));
  Y.setup_x(all_parameters.returnInteger("m"),0.,all_parameters.returnDouble("Y_max"));
  Q.setup_x(all_parameters.returnInteger("kmax"),0.,all_parameters.returnDouble("Q_max")*all_parameters.returnDouble("power_max"));
  C.setup_x(all_parameters.returnInteger("no_of_commits"),all_parameters.returnDouble("minimum_commit"),all_parameters.returnDouble("maximum_commit"));
  #ifdef DEBUG
  std::cout << " Update options...\n";
  #endif
  if(!first_run)
  {
	#ifdef DEBUG
	std::cout << " Interpolate options...\n";
	#endif
	// interpolate from old values
	option_max_old=option_max;
	option_max.update(&X,&Y,&Q);
	option_max.intrplt(option_max_old);
	option_max_old.update(&X,&Y,&Q);
	
	option_commits=option_commits_new;
	option_commits_new.update(&C,&X,&Y,&Q);
	option_commits_new.intrplt(option_commits);
	option_commits.update(&C,&X,&Y,&Q);
	option_commits_old.update(&C,&X,&Y,&Q);
  }
  else
  {
	// update option storage
	option_max.update(&X,&Y,&Q);
	option_max_old.update(&X,&Y,&Q);
	option_commits_new.update(&C,&X,&Y,&Q);
	option_commits.update(&C,&X,&Y,&Q);
	option_commits_old.update(&C,&X,&Y,&Q);
  }
  
  #ifdef DEBUG
  std::cout << " Storage update successful...\n";
  #endif
  
}

void Delay_Op2::default_print_values()
{
  // update printer
  file_print.setup_z(&Q);
  file_print.setup_x(&X);
  file_print.setup_y(&Y);
  // 
  file_print4.setup_x(&C);
  file_print4.setup_w(&Q);
  file_print4.setup_y(&X);
  file_print4.setup_z(&Y);
  
  std::ofstream default_file((config_dir+str_default_option3).c_str());
  print_config.print_all_params(&default_file);
  default_file.close();
  
  default_file.open((config_dir+str_default_option4).c_str());
  print_config4.print_all_params(&default_file);
  default_file.close();
  
  #ifdef DEBUG
  std::cout << " Print update successful...\n";
  #endif

}

void Delay_Op2::link_print_config()
{
  // printing option max
  file_print.link_printer_config(&print_config);
  file_print.link_printer_config2(&print_config);
  file_print.link_printer_config3(&print_config);
  // printing option commits
  file_print4.link_printer_config(&print_config4);
  file_print4.link_printer_config2(&print_config4);
  file_print4.link_printer_config3(&print_config4);
  file_print4.link_printer_config4(&print_config4);
  #ifdef DEBUG
  std::cout << " link print config complete...\n";
  #endif
}

void Delay_Op2::default_print_values(std::string filename,std::string filename4)
{
  
  if(print_config.read_all_params(config_dir+filename) and print_config4.read_all_params(config_dir+filename4))return;
  else 
  {
	std::cout << " Use default print values...\n";
	default_print_values();
  }
}

void Delay_Op2::default_values(std::string filename)
{
  if(all_parameters.read_all_params(config_dir+filename))return;
  else
  {
	std::cout << " File read failed using default parameter values...\n";
	default_values();
  }
}

void Delay_Op2::default_values(void)
{
  // default values
  all_parameters.reset();
  // write default values to file
  std::ofstream default_file((config_dir+str_default_parameters).c_str());
  if(!default_file)
  {
	std::cout << " no default file created..\n";
	return;
  }
  all_parameters.print_all_params(&default_file);
  default_file.close();
}

void Delay_Op2::solve(int no_of_hours,string filename)
{	
  update_storage();
  solve(no_of_hours,filename,true,true);
}

void Delay_Op2::solve(int no_of_hours,string filename,bool display,bool reset)
{	
  if(reset or first_run)
  {
	option_max=0.;
	option_max_old=0.;
	option_commits=0.;
	option_commits_new=0.;
	option_commits_old=0.;
  }
  first_run=false;
  
  if(display)
	{
	  std::cout << " Solving for a multi-time option.\n No of commits :: " 
	  << C.size();
	  std::cout << "\n No of contracts :: " << no_of_hours << std::endl;
	}
	
	Countdown screen_output;
	
	double tau;
	stringstream str_stream;
	string temp,temp4;
	
	all_parameters.changeParameter(all_parameters.returnDouble("contract_length"),"maturity");
	
	for(int hours=0;hours<no_of_hours;hours++)
	{
	  #ifdef DEBUG
	  std::cout << " setup initial conditions for all C \n";
	  #endif
	  // initial conditions
	  if(market_delay)option_commits = option_commits_new;
	  else for(int i=0;i<C.size();i++)option_commits.copy(i,option_max);
	  // Write progress to screen...
	  if(display)screen_output.update(hours,no_of_hours);
	  // progress written
	  str_stream.str("");
	  str_stream << no_of_hours - hours - 1 ;
	  temp = data_dir + filename + "." + str_stream.str();
	  temp4 = data_dir + filename + ".delay." + str_stream.str();
	  tau = -hours*all_parameters.returnDouble("contract_length"); // change contract length
	  // solve for all commits
	  // set current commit level
	  all_parameters.changeParameter(C[0],"C_commit");
	  #ifdef DEBUG
	  std::cout << " i :: " << 0 << " j :: " << 0 << " tau :: " << tau << "\n";
	  #endif
	  option_pricer->runSolver(option_commits[0],tau);
	  // reset option max
	  option_max = option_pricer->return_option_new();
	  if(market_delay){
		// if we have the delayed option, need to run C.size() times for all possible initial conditions 'option_commits[j]'
		for(int j=1;j<C.size();j++)
		{
		  #ifdef DEBUG
		  std::cout << " i :: " << 0 << " j :: " << j << " tau :: " << tau << "\n";
		  #endif
		  option_pricer->runSolver(option_commits[j],tau);
		  // += is overloaded to maximum of two options
		  option_max += option_pricer->return_option_new();
		}
		// initial conditions for next timestep, maximised over j which is choice of C at previous timestep
		option_commits_new.copy(0,option_max);
	  }
	  for(int i=1;i<C.size();i++)
	  {
		//set current commit level
		all_parameters.changeParameter(C[i],"C_commit");
		if(market_delay)
		{
		  #ifdef DEBUG
		  std::cout << " i :: " << i << " j :: " << 0 << " tau :: " << tau << "\n";
		  #endif
		  // option having chosen C_0 at previous timestep
		  option_pricer->runSolver(option_commits[0],tau);
		  //reset option max
		  option_max = option_pricer->return_option_new();
		}
		#ifdef DEBUG
		std::cout << " i :: " << i << " j :: " << i << " tau :: " << tau << "\n";
		#endif
		option_pricer->runSolver(option_commits[i],tau);
		// += is overloaded to maximum of two options
		option_max += option_pricer->return_option_new();
		if(market_delay)
		{
		  for(int j=1;j<C.size();j++)
		  {
			// already calculated with i==j
			if(i==j)continue;
			#ifdef DEBUG
			std::cout << " i :: " << i << " j :: " << j << " tau :: " << tau << "\n";
			#endif
			option_pricer->runSolver(option_commits[j],tau);
			option_max += option_pricer->return_option_new();
		  }
		  // initial conditions for next timestep, maximised over j which is choice of C at previous timestep for commit C_i
		  option_commits_new.copy(i,option_max);
		}
	  }
	  if(display)
	  {
		if(market_delay)
		{
		  file_print4.change_filename(temp4);
		  file_print4.print_y_z(option_commits_new);
		}
		else
		{
		  file_print.change_filename(temp);
		  file_print.print_x_y(option_max);
		}
	  }
	  
	}

	if(display)
	{
		screen_output.update(no_of_hours,no_of_hours);
		std::cout << " \n Option solved.\n";
	}

}

// Try and extrapolate to find the perpetual option value quickly
double Delay_Op2::perpetual_pricer(bool print_results,bool run_starter)
{
  update_storage();
  
  double total_time,average,per_val,average_old,tolerance;
  double x,y,q,c;
  int contracts_in_a_day;
  // tolerance for convergence
  tolerance = 1.e-6;
  // use middle grid point
  x = option_pricer->average_X_t(0.);
  y = option_pricer->average_Y_t(0.);
  q = (Q.max()-Q.min())/2.;
  c = (C.max()-C.min())/2.;
  
  // contracts = period/contract_length
  contracts_in_a_day = int( all_parameters.returnDouble("period_Y")
  /all_parameters.returnDouble("contract_length") +0.5);
  // if the option is not periodic method won't work
  if(24%int(all_parameters.returnDouble("contract_length")+0.5)!=0)return 0.;
  // If running for the first time, run a single day option to estimate value
  if(run_starter){
	// reset initial conditions
	if(market_delay)option_commits_new = 0.;
	else option_max = 0.;
	// solve once for initial estimate
	solve(contracts_in_a_day,string("test"),false,true);
	// time = contracts * contract_length
	total_time = contracts_in_a_day*all_parameters.returnDouble("contract_length");
	// calculate average
	if(market_delay)average = option_commits(c,x,y,q);
	else average = option_max(x,y,q);
	// perpetual = average/interest_rate
	per_val = estimate_perpetual_option(total_time,average);
	// debug messages
	if(debug){
	  cout << " Average value : " << average << " after " << total_time << " hours\n";
	  cout << "Resulting perpetual option :: " << per_val << endl;
	}
	// set up new initial condition
	//      w = per_val
	if(market_delay)option_commits_new = per_val;
	else option_max = per_val;
  }
  // run through a few iteration at first to let the solution settle into periodic solution
  average = per_val;
  for(int loop=0;loop<6;loop++){
	// average old allows us to calculate the value added in a day
	average_old = average;
	// solve for the day
	solve(contracts_in_a_day,string("perpetual."),false,false);
	// calculate average
	if(market_delay)average = option_commits_new(c,x,y,q);
	else average = option_max(x,y,q);
	// debug messages
	if(debug){
	  cout << setprecision(8);
	  cout << "it : " << loop << " :: average :: " << 
	  average << " :: ratio :: " << 1.-average/average_old << endl;
	}
  }
  
  for(int big_loop=0;big_loop<200;big_loop++){
	// calculate value of average coupon payment each day
	total_time = contracts_in_a_day*all_parameters.returnDouble("contract_length");
	// value added in one day is average - average_old
	average = average - average_old;
	// per_val is the value of the perpetual option adding that value every day
	per_val = estimate_perpetual_option(total_time,average);
	// debug messages
	if(debug){
	  cout << " Average value : " << average << " after " << total_time << " hours\n";
	  cout << "Resulting perpetual option :: " << per_val << endl;
	}
	// break if the adjustment is sufficiently small
	// calculate average
	if(market_delay)average = option_commits_new(c,x,y,q);
	else average = option_max(x,y,q);
	if(rms_error()<tolerance)break;
	
	// adjust initial conditions so that
	//      w = w + per_val
	if(market_delay)option_commits_new += per_val;
	else option_max += per_val;
	// loop again to settle solution
	for(int loop=0;loop<4;loop++){
	  average_old = average;
	  //
	  if(market_delay)option_commits_old = option_commits_new;
	  else option_max_old = option_max;
	  solve(contracts_in_a_day,string("perpetual."),false,false);
	  if(market_delay)average = option_commits_new(c,x,y,q);
	  else average = option_max(x,y,q);
	  if(debug){
		cout << setprecision(8);
		cout << "it : " << loop << " :: average :: " << 
		average << " :: ratio :: " << 1.-average/average_old <<
		" RMS error :: " << rms_error() << endl;
	  }
	}
	
  }
  
  if(print_results)
  {
	solve(contracts_in_a_day,string("perpetual"),true,false);
	std::cout << "\n Annual return on the option is :: " << annual_returns() << "\n";
  }
  
  return annual_returns();
  
}

double Delay_Op2::rms_error()
{
  if(market_delay)return rms(option_commits_new,option_commits_old);
  else return rms(option_max,option_max_old);
}

double Delay_Op2::estimate_perpetual_option(double time,double value)
{
	
	double interest_rate_per_hour;
	// interest rate per hour
	interest_rate_per_hour = all_parameters.returnDouble("interest_rate")/24./365.;
	
	return value/time/interest_rate_per_hour;
	
}

void Delay_Op2::print_params()
{
  all_parameters.print_all_params(&std::cout);
}

double Delay_Op2::annual_returns()
{
  double x,y,c,q;
  x = option_pricer->average_X_t(0.);
  y = option_pricer->average_Y_t(0.);
  q = (Q.max()-Q.min())/2.;
  c = (C.max()-C.min())/2.;
  if(market_delay) return
	all_parameters.returnDouble("interest_rate")*option_commits(c,x,y,q);
  else return 
	all_parameters.returnDouble("interest_rate")*option_max(x,y,q);
}

