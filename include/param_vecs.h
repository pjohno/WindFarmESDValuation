#ifndef _param_vecs_h_included_
#define _param_vecs_h_included_

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "pj_data.h"

// #define DEBUG
namespace PJLib {

class Param_vecs
{

  public:
	
	// create a container for all pointers to the integer parameters
	std::vector<PJ_data<int> > the_ints;
	// create a container for all pointers to the real parameters
	std::vector<PJ_data<double> > the_reals;
	// create a container for all pointers to the string parameters
	std::vector<PJ_data<std::string> > the_strings;

	// constructor
	Param_vecs(){};
	// empty all parameters
	void clear_vecs(void);
	// set up new parameter - int
	void assignParameter(int *xx,std::string name);
	// - double
	void assignParameter(double *xx,std::string name);
	// - string
	void assignParameter(std::string *xx,std::string name);
	// change parameter value - int
	void changeParameter(int xx,int parameter_ref);
	// - double
	void changeParameter(double xx,int parameter_ref);
	// - string
	void changeParameter(std::string xx,int parameter_ref);
	// return value - int
	int returnInteger(int parameter_ref);
	// - double
	double returnDouble(int parameter_ref);
	// - string
	std::string returnString(int parameter_ref);
	// print all parameters to stream
	virtual void print_all_params(std::ostream *output);
	
};

class Param_vecs2:public Param_vecs
{
  public:
	// print all parameters to stream
	void print_all_params(std::ostream *output);
	// read all parameters to stream
	bool read_all_params(std::string filename);
  
};

}//namespace

#endif
