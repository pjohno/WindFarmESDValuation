#include "optimizer.h"
#include "storage.h"
#include "no_price.h"
#include "finite_diff_methods.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

Test_Optimizer::Test_Optimizer(int no_of_options,bool progress_bar):Optimizer(no_of_options,progress_bar)
{
  
};

Test_Optimizer::Test_Optimizer(int no_of_options,bool progress_bar,bool no_price):Optimizer(no_of_options,progress_bar,no_price)
{

};

void Test_Optimizer::perpetual_pricer(bool print_results)
{
  
  bool debug=true;
  double total_time,average,per_val,average_old,tolerance;
  int contracts_in_a_day;
  
  tolerance = 1.e-6;
  
  // contracts = period/contract_length
  contracts_in_a_day = int( all_parameters[0]->returnDouble(31)
  /all_parameters[0]->returnDouble(22) +0.5);
  if(24%int(all_parameters[0]->returnDouble(22)+0.5)!=0)return;
  // solve once for initial estimate
  solve_multiple_time(contracts_in_a_day,string("test."),true,false);
  // time = contracts * contract_length
  total_time = contracts_in_a_day*all_parameters[0]->returnDouble(22);
  // 	cout << we_have_option_max<<endl;
  // calculate average
  average = option_average();
  // perpetual = average/interest_rate
  per_val = estimate_perpetual_option(total_time,average);
  
  if(debug){
	cout << " Average value : " << average << " after " << total_time << " hours\n";
	cout << "Resulting perpetual option :: " << per_val << endl;
  }
  // set up new initial condition
  //      w = per_val
  set_perpetual_option_conditions(per_val);
  // run through a few iteration at first to let the solution settle into periodic solution
  average = per_val;
  for(int loop=0;loop<10;loop++){
	average_old = average;
	solve_multiple_time(contracts_in_a_day,string("perpetual."),false,true);
	average = option_average();
	if(debug){
	  cout << setprecision(8);
	  cout << "it : " << loop << " :: average :: " << 
	  average << " :: ratio :: " << 1.-average/average_old << endl;
	}
  }
  
  for(int big_loop=0;big_loop<200;big_loop++){
	// calculate value of average coupon payment each day
	total_time = contracts_in_a_day*all_parameters[0]->returnDouble(22);
	average = average - average_old;
	per_val = estimate_perpetual_option(total_time,average);
	if(debug){
	  cout << " Average value : " << average << " after " << total_time << " hours\n";
	  cout << "Resulting perpetual option :: " << per_val << endl;
	}
	// break if the adjustment is sufficiently small
	average = option_average();
	if(fabs(per_val/average)<tolerance)break;
	
	// adjust initial conditions so that
	//      w = w + per_val
	adjust_perpetual_option_conditions(per_val);
	// loop again to settle solution
	for(int loop=0;loop<5;loop++){
	  average_old = average;
	  *option_temp = *option_max_new;
	  solve_multiple_time(contracts_in_a_day,string("perpetual."),false,true);
	  average = option_average();
	  if(debug){
		cout << setprecision(8);
		cout << "it : " << loop << " :: average :: " << 
		average << " :: ratio :: " << 1.-average/average_old <<
		" RMS error :: " << rms_error() << endl;
	  }
	  
	}
	
  }
  
  if(print_results)solve_multiple_time(contracts_in_a_day,string("perpetual."),true,true);
	
}

