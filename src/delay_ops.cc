#include "param_vecs.h"
#include "delay_ops.h"
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

// #define DEBUG

const std::string root_dir=std::string(DATA_LOCATION);
const std::string config_dir=root_dir+"Config/";
const std::string data_dir=root_dir+"Data/";
const std::string str_default_parameters="default_parameters.in";
const std::string str_default_option3="print_option3.in";
const std::string str_default_option4="print_option4.in";
const std::string str_default_store="default_store.in";
const std::string str_default_wind="default_wind.in";
const std::string str_default_trading="default_trading.in";
const std::string str_default_price="default_price.in";
const std::string str_default_grid="default_grid.in";

Delay_Op::~Delay_Op()
{
  #ifdef DEBUG
  std::cout << " Delay_Op clear parameters...\n";std::cout.flush();
  #endif
  all_parameters.clear_vecs();
  #ifdef DEBUG
  std::cout << " delete storage option...\n";
  #endif
  delete option_pricer;
  #ifdef DEBUG
  std::cout << " clean up complete...\n";
  #endif
}  

Delay_Op::Delay_Op()
{
  option_pricer = new Storage(&all_parameters,false);
  all_linking();first_run=true;
  input_values(str_default_parameters,str_default_store,str_default_wind,str_default_price,str_default_trading,str_default_grid);
}

Delay_Op::Delay_Op(std::string base,std::string store,std::string wind,std::string price,std::string trading,std::string grid_size)
{
  option_pricer = new Storage(&all_parameters,false);
  all_linking();first_run=true;
  // use 
  input_values(base,store,wind,price,trading,grid_size);
}

void Delay_Op::input_values(std::string base,std::string store,std::string wind,std::string price,std::string trading,std::string grid_size)
{
  #ifdef DEBUG
  std::cout << " Use the following configuration files...\n";
  std::cout << base << " " <<store<<" "<<wind<<" "<<price<<" "<<trading<<" "<<grid_size<<" "<<str_default_option3<< ""<<str_default_option4 << "\n";
  #endif
  default_values(base);
  default_store_values(store);
  default_wind_values(wind);
  default_price_values(price);
  default_trading_values(trading);
  default_grid_values(grid_size);
  default_print_values(str_default_option3,str_default_option4);
}

void Delay_Op::update_storage()
{
  #ifdef DEBUG
  std::cout << " Update grids...\n";
  #endif
  // update grids
  X.setup_x(all_parameters.returnInteger(1),0.,all_parameters.returnDouble(1));
  Y.setup_x(all_parameters.returnInteger(2),0.,all_parameters.returnDouble(2));
  Q.setup_x(all_parameters.returnInteger(3),0.,all_parameters.returnDouble(3)*all_parameters.returnDouble(20));
  C.setup_x(all_parameters.returnInteger(5),all_parameters.returnDouble(34),all_parameters.returnDouble(33));
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

void Delay_Op::default_print_values()
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

void Delay_Op::link_print_config()
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

void Delay_Op::default_print_values(std::string filename,std::string filename4)
{
  if(print_config.read_all_params(config_dir+filename) and print_config4.read_all_params(config_dir+filename4))return;
  else 
  {
	std::cout << " Use default print values...\n";
	default_print_values();
  }
}

void Delay_Op::link_store_config()
{
  // link up store config to 
  store_config.assignParameter(all_parameters.the_reals[3].return_pointer()," Maximum size of store ");
  store_config.assignParameter(all_parameters.the_reals[12].return_pointer()," Energy Efficiency of store ");
  store_config.assignParameter(all_parameters.the_reals[13].return_pointer()," Rate of charge ");
  store_config.assignParameter(all_parameters.the_reals[14].return_pointer()," Rate of discharge ");
  store_config.assignParameter(all_parameters.the_reals[15].return_pointer()," Smoothing on charge ");
  store_config.assignParameter(all_parameters.the_reals[16].return_pointer()," Smoothing on discharge ");
  #ifdef DEBUG
  std::cout << " link store parameters complete...\n";
  #endif
}

void Delay_Op::default_store_values(std::string filename)
{
  if(store_config.read_all_params(config_dir+filename))return;
  else 
  {
	std::cout << " Use default store values...\n";
	default_store_values();
  }
}


void Delay_Op::default_store_values()
{
  std::ofstream default_file((config_dir+str_default_store).c_str());
  store_config.print_all_params(&default_file);
  default_file.close();
}


void Delay_Op::link_wind_config()
{
  // link up store config to 
  wind_config.assignParameter(all_parameters.the_reals[1].return_pointer()," Maximum wind speed X_max ");
  wind_config.assignParameter(all_parameters.the_reals[6].return_pointer()," Mean Reversion in X ");
  wind_config.assignParameter(all_parameters.the_reals[7].return_pointer()," Long-term average in X ");
  wind_config.assignParameter(all_parameters.the_reals[8].return_pointer()," Variance in X ");
  wind_config.assignParameter(all_parameters.the_reals[27].return_pointer()," Amplitude in sin X ");
  wind_config.assignParameter(all_parameters.the_reals[28].return_pointer()," Period in sin X ");
  wind_config.assignParameter(all_parameters.the_reals[29].return_pointer()," Phase in X ");
  
  #ifdef DEBUG
  std::cout << " link wind parameters complete...\n";
  #endif
}

void Delay_Op::default_wind_values(std::string filename)
{
  if(wind_config.read_all_params(config_dir+filename))return;
  else 
  {
	std::cout << " Use default wind values...\n";
	default_wind_values();
  }
}


void Delay_Op::default_wind_values()
{
  std::ofstream default_file((config_dir+str_default_wind).c_str());
  wind_config.print_all_params(&default_file);
  default_file.close();
}

void Delay_Op::link_trading_config()
{
  market_delay=0;
  no_of_commits=5;
  all_parameters.assignParameter(&no_of_commits,"No of commits :: ");
  // link up store config to 
  trading_config.assignParameter(&market_delay," Market Delay ");
  trading_config.assignParameter(all_parameters.the_reals[5].return_pointer()," Interest Rate ");
  trading_config.assignParameter(all_parameters.the_reals[22].return_pointer()," Contract Length ");
  trading_config.assignParameter(all_parameters.the_reals[20].return_pointer()," Maximum Power of the Farm ");
  trading_config.assignParameter(all_parameters.the_reals[34].return_pointer()," Minimum Commit ");
  trading_config.assignParameter(all_parameters.the_reals[33].return_pointer()," Maximum Commit ");
  trading_config.assignParameter(all_parameters.the_reals[35].return_pointer()," Buy/Sell Spread ");
  
  #ifdef DEBUG
  std::cout << " link trading parameters complete...\n";
  #endif
}

void Delay_Op::default_trading_values(std::string filename)
{
  if(trading_config.read_all_params(config_dir+filename))return;
  else 
  {
	std::cout << " Use default trading values...\n";
	default_trading_values();
  }
}


void Delay_Op::default_trading_values()
{
  std::ofstream default_file((config_dir+str_default_trading).c_str());
  trading_config.print_all_params(&default_file);
  default_file.close();
}

void Delay_Op::link_price_config()
{
  // link up store config to 
  price_config.assignParameter(all_parameters.the_reals[2].return_pointer()," Maximum wind price Y_max ");
  price_config.assignParameter(all_parameters.the_reals[9].return_pointer()," Mean Reversion in Y ");
  price_config.assignParameter(all_parameters.the_reals[10].return_pointer()," Average Y ");
  price_config.assignParameter(all_parameters.the_reals[11].return_pointer()," Variance in Y ");
  price_config.assignParameter(all_parameters.the_reals[30].return_pointer()," Amplitude in sin Y ");
  price_config.assignParameter(all_parameters.the_reals[31].return_pointer()," Period in sin Y ");
  price_config.assignParameter(all_parameters.the_reals[32].return_pointer()," Phase in sin Y ");
  
  #ifdef DEBUG
  std::cout << " link price parameters complete...\n";
  #endif
}

void Delay_Op::default_price_values(std::string filename)
{
  if(price_config.read_all_params(config_dir+filename))return;
  else 
  {
	std::cout << " Use default price values...\n";
	default_price_values();
  }
}


void Delay_Op::default_price_values()
{
  std::ofstream default_file((config_dir+str_default_price).c_str());
  price_config.print_all_params(&default_file);
  default_file.close();
}

void Delay_Op::link_grid_config()
{
  debug=0;
  // link up store config to 
  grid_config.assignParameter(all_parameters.the_ints[1].return_pointer()," Points in X ");
  grid_config.assignParameter(all_parameters.the_ints[2].return_pointer()," Points in Y ");
  grid_config.assignParameter(all_parameters.the_ints[3].return_pointer()," Points in Q ");
  grid_config.assignParameter(&no_of_commits," Points in C ");
  grid_config.assignParameter(&debug," Show debugging messages");
  
  #ifdef DEBUG
  std::cout << " link grid parameters complete...\n";
  #endif
}

void Delay_Op::default_grid_values(std::string filename)
{
  if(grid_config.read_all_params(config_dir+filename))return;
  else 
  {
	std::cout << " Use default grid values...\n";
	default_grid_values();
  }
}


void Delay_Op::default_grid_values()
{
  std::ofstream default_file((config_dir+str_default_grid).c_str());
  grid_config.print_all_params(&default_file);
  default_file.close();
}

void Delay_Op::default_values(std::string filename)
{
  if(all_parameters.read_all_params(config_dir+filename))return;
  else 
  {
	std::cout << " Use default parameter values...\n";
	default_values();
  }
}

void Delay_Op::default_values(void)
{
  
  // Set up base parameters...
  all_parameters.changeParameter(101,1); // change nodes in X	
  all_parameters.changeParameter(11,2); // change nodes in Y
  all_parameters.changeParameter(11,3); // change nodes in Q
  all_parameters.changeParameter(9,5); // change no of commits
  all_parameters.changeParameter(0.5,0); // change commit level - percentage of power output
  all_parameters.changeParameter(100.,1); // change X max
  all_parameters.changeParameter(100.,2); // change Y max
  all_parameters.changeParameter(1.,3); // change Q max - percentage of power output
  all_parameters.changeParameter(1.,4); // change maturity
  all_parameters.changeParameter(0.05,5); // change interest rate
  all_parameters.changeParameter(0.1,6); // change mean rev in X
  all_parameters.changeParameter(8.,7); // change average in X
  all_parameters.changeParameter(0.2,8); // change variance in X
  all_parameters.changeParameter(0.04,9); // change mean rev in Y
  all_parameters.changeParameter(40.0,10); // change average in Y
  all_parameters.changeParameter(0.075,11); // change variance in Y
  all_parameters.changeParameter(0.8366,12); // change energy efficiency
  all_parameters.changeParameter(0.25,13); // change charge - percentage of power output
  all_parameters.changeParameter(0.25,14); // change discharge - percentage of power output
  all_parameters.changeParameter(0.2,15); // change smoothing on charge
  all_parameters.changeParameter(0.2,16); // change smoothing on discharge
  all_parameters.changeParameter(0.9,17); // change contract sell price 
  all_parameters.changeParameter(1.1,18); // change contract bid price
  all_parameters.changeParameter(25.,19); // change power cut out
  all_parameters.changeParameter(1.,20); // maximum power output 
  all_parameters.changeParameter(10.,21); // change power turn
  all_parameters.changeParameter(1.,22); // change contract length
  all_parameters.changeParameter(0.375,27); // amp of X
  all_parameters.changeParameter(-2.,29); // phase of X
  all_parameters.changeParameter(0.375,30); // amp of Y
  all_parameters.changeParameter(-14.,32); // phase of Y
  all_parameters.changeParameter(1.,33); // maximum commit
  all_parameters.changeParameter(0.01,34); // minimum commit
  all_parameters.changeParameter(1.,35); // buy sell spread
  
  std::ofstream default_file((config_dir+str_default_parameters).c_str());
  all_parameters.print_all_params(&default_file);
  default_file.close();
}

void Delay_Op::solve(int no_of_hours,string filename)
{
  
  update_storage();
  solve(no_of_hours,filename,true,true);
  
}

void Delay_Op::solve(int no_of_hours,string filename,bool display,bool reset)
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
	  tau = -hours*all_parameters.returnDouble(22); // change contract length
	  // solve for all commits
	  // set current commit level
	  all_parameters.changeParameter(C[0],0);
	  #ifdef DEBUG
	  std::cout << " Commit :: " << C[0] << " tau :: " << tau << "\n";
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
		all_parameters.changeParameter(C[i],0);
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
double Delay_Op::perpetual_pricer(bool print_results,bool run_starter)
{
  update_storage();
  
  double total_time,average,per_val=1.,average_old=1.,tolerance;
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
  contracts_in_a_day = int( all_parameters.returnDouble(31)
  /all_parameters.returnDouble(22) +0.5);
  // if the option is not periodic method won't work
  if(24%int(all_parameters.returnDouble(22)+0.5)!=0)return 0.;
  // If running for the first time, run a single day option to estimate value
  if(run_starter){
	// reset initial conditions
	if(market_delay)option_commits_new = 0.;
	else option_max = 0.;
	// solve once for initial estimate
	solve(contracts_in_a_day,string("test"),false,true);
	// time = contracts * contract_length
	total_time = contracts_in_a_day*all_parameters.returnDouble(22);
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
	total_time = contracts_in_a_day*all_parameters.returnDouble(22);
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

double Delay_Op::rms_error()
{
  if(market_delay)return rms(option_commits_new,option_commits_old);
  else return rms(option_max,option_max_old);
}

double Delay_Op::estimate_perpetual_option(double time,double value)
{
	
	double interest_rate_per_hour;
	// interest rate per hour
	interest_rate_per_hour = all_parameters.returnDouble(5)/24./365.;
	
	return value/time/interest_rate_per_hour;
	
}

void Delay_Op::print_params(std::ostream &output)
{
  all_parameters.print_all_params(&output);
}

double Delay_Op::annual_returns()
{
  double x,y,c,q;
  x = option_pricer->average_X_t(0.);
  y = option_pricer->average_Y_t(0.);
  q = (Q.max()-Q.min())/2.;
  c = (C.max()-C.min())/2.;
  if(market_delay) return
	all_parameters.returnDouble(5)*option_commits(c,x,y,q);
  else return 
	all_parameters.returnDouble(5)*option_max(x,y,q);
}

