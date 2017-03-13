#ifndef _test_optimizer_h_included_
#define _test_optimizer_h_included_

/* Set up a class for the option value */

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

class Optimizer

class Test_Optimizer:public Optimizer
{

  public:
	
	// constructor, specifying number of commits
	Test_Optimizer(int no_of_options,bool progress_bar);
	// constructor, using no price variable (returning 1)
	Test_Optimizer(int no_of_options,bool progress_bar,bool no_price);
	
	void perpetual_pricer(bool print_results);
	

};
