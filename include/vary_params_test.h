#ifndef _vary_params_h_included_
#define _vary_params_h_included_

/* Set up a class for the option value */

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "delay_ops_test.h"

class Vary_Params2: public Delay_Op2
{
  
  // linking
  void all_linking(){link_filenames();};
  void link_filenames();
  static const std::string str_default_filenames,data_dir;
  
  protected:
	
	// iteration parameters
	int grid_consistent,grid_or_param,n1,n2;
	std::string p1,p2;
	Fixed_Grid X1,X2;
	double min1,min2,max1,max2;
	// store config filenames
	ParameterMaps filenames;
	// config filenames
	std::string base,store,wind,price,trading,grid,farm;
	// results filename
	std::string results_filename;
	// use defaults
	void default_filenames();
	// update configuration
	void update();
	// read parameters from file
	void default_filenames(std::string filename);
	
  public:
	
	// deconstructor
	~Vary_Params2(){};
	// constructor
	Vary_Params2();
	// iterate over parameters
	virtual void iterator();
	// iterate over grid
	virtual void iterator_grid();
	// input parameters
	void input(std::string filename);
	// run iteration from file
	void runIterator(){if(grid_or_param)iterator_grid();else iterator();};
	// update local grids
	void update_local_grids(Fixed_Grid &X,Fixed_Grid &Y,Fixed_Grid &Q,Fixed_Grid &C);
	// check for consistency in grid
	void check_grid_consistency(Fixed_Grid &X,Fixed_Grid &Y,Fixed_Grid &Q,Fixed_Grid &C);

};

#endif
