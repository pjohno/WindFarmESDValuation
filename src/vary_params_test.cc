
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "vary_params_test.h"
#include "fixed_grid.h"

#define DEBUG

const std::string Vary_Params2::str_default_filenames="default_filenames.in";
const std::string Vary_Params2::data_dir="Parameter_Sensitivity/";

Vary_Params2::Vary_Params2()
{
  all_linking();
  input(str_default_filenames);
}

void Vary_Params2::link_filenames()
{
  // link parameters for iteration
  filenames.assignParameter(&grid_or_param,"grid_or_param"," Iterate over the grid or the parameters ",0);
  filenames.assignParameter(&grid_consistent,"grid_consistent"," Check for grid consistency at each iteration ",0);
  filenames.assignParameter(&p1,"p1"," Choose parameter for first iteration ","X_max");
  filenames.assignParameter(&p2,"p2"," Choose parameter for second iteration ","Y_max");
  filenames.assignParameter(&n1,"n1"," Number of its (1) ",2);
  filenames.assignParameter(&n2,"n2"," Number of its (2) ",2);
  filenames.assignParameter(&min1,"min1"," Minimum Value (1) ",50.);
  filenames.assignParameter(&min2,"min2"," Minimum Value (2) ",80.);
  filenames.assignParameter(&max1,"max1"," Maximum value (1) ",100.);
  filenames.assignParameter(&max2,"max2"," Maximum value (2) ",160.);
  //  link filename parameters to use
  filenames.assignParameter(&base,"base"," Base Parameters ","default.in");
  filenames.assignParameter(&store,"store"," Store Parameters ","default.in");
  filenames.assignParameter(&wind,"wind"," Wind Parameters ","default.in");
  filenames.assignParameter(&price,"price"," Price Parameters ","default.in");
  filenames.assignParameter(&trading,"trading"," Trading Parameters ","default.in");
  filenames.assignParameter(&grid,"grid"," Grid Parameters ","default.in");
  filenames.assignParameter(&farm,"farm"," Wind Farm Parameters ","default.in");
  filenames.assignParameter(&results_filename,"results_filename"," Results Filename ","default");
  #ifdef DEBUG
  std::cout << " link filename parameters complete...\n";
  #endif
}

void Vary_Params2::default_filenames(std::string filename)
{
  if(filenames.read_all_params(config_dir+filename))return;
  else 
  {
	std::cout << " Use default filename values...\n";
	default_filenames();
  }
}

void Vary_Params2::default_filenames()
{
  filenames.reset();
  std::ofstream default_file((config_dir+str_default_filenames).c_str());
  filenames.print_all_params(&default_file);
  default_file.close();
}

void Vary_Params2::update()
{
  default_values(base);
  default_values(store);
  default_values(wind);
  default_values(price);
  default_values(trading);
  default_values(grid);
  default_values(farm);
  X1.setup_x(n1,min1,max1);
  X2.setup_x(n2,min2,max2);
}

void Vary_Params2::iterator()
{
  Fixed_Grid X_hat,Y_hat,Q_hat,C_hat;
  update_local_grids(X_hat,Y_hat,Q_hat,C_hat);
  
  bool debug_here=true;
  // open streams for writing data to file
  std::ofstream vary_param,log_file;
  std::string results;
  results = data_dir+"vary_param."+results_filename+".dat";
  vary_param.open(results.c_str());
  vary_param << setprecision(8);
  // use a log file to check progress
  results = results_filename+".log";
  log_file.open(results.c_str());
  filenames.print_all_params(&log_file);
  
  // iterate over parameter 1
  for(int i = 0;i<X1.size();i++)
  {
	log_file << "\n Iteration... \n";
	for(int j = 0;j<X2.size();j++)
	{
	  // change parameters
	  all_parameters.changeParameter(X1[i],p1);
	  all_parameters.changeParameter(X2[j],p2);
	  // check parameters
	  log_file << "Parameters :: ";
	  all_parameters.print(p1,&log_file);log_file << " :: ";
	  all_parameters.print(p2,&log_file);log_file <<"\n";
	  log_file.flush();
	  if(debug_here){
		std::cout << "Parameters :: ";
		all_parameters.print(p1,&std::cout);std::cout << " :: ";
		all_parameters.print(p2,&std::cout);std::cout <<"\n";
		std::cout.flush();
	  }
	  if(grid_consistent)check_grid_consistency(X_hat,Y_hat,Q_hat,C_hat);
	  // write annual returns to file
	  vary_param << X1[i] << " " << X2[j] << " ";
	  vary_param << perpetual_pricer(false,false) << "\n";
	  vary_param.flush();
	}
	vary_param << " \n" ;
  }
  vary_param.close();
  log_file.close();
}

void Vary_Params2::iterator_grid()
{
  bool debug_here=true;
  // open streams for writing data to file
  std::ofstream vary_param,log_file;
  std::string results;
  results = data_dir+"vary_grid."+results_filename+".dat";
  vary_param.open(results.c_str());
  vary_param << setprecision(8);
  // use a log file to check progress
  results = results_filename+".log";
  log_file.open(results.c_str());
  filenames.print_all_params(&log_file);
  
  // iterate over parameter 1
  for(int i = 0;i<X1.size();i++)
  {
	log_file << "\n Iteration... \n";
	// change parameters
	all_parameters.changeParameter(int(X1[i]+0.5),p1);
	// check parameters
	log_file << "Parameters :: ";
	all_parameters.print(p1,&log_file);log_file <<"\n";
	log_file.flush();
	if(debug_here){
	  std::cout << "Parameters :: ";
	  all_parameters.print(p1,&std::cout);std::cout <<"\n";
	  std::cout.flush();
	}
	// write annual returns to file
	vary_param << X1[i] << " ";
	vary_param << perpetual_pricer(false,false) << "\n";
  }
  vary_param.close();
  log_file.close();
  
}

void Vary_Params2::input(std::string filename)
{
  default_filenames(filename);
  update();
}


void Vary_Params2::update_local_grids(Fixed_Grid &X_,Fixed_Grid &Y_,Fixed_Grid &Q_,Fixed_Grid &C_)
{
  X_.setup_x(all_parameters.returnInteger("n"),0.,all_parameters.returnDouble("X_max"));
  Y_.setup_x(all_parameters.returnInteger("m"),0.,all_parameters.returnDouble("Y_max"));
  Q_.setup_x(all_parameters.returnInteger("kmax"),0.,all_parameters.returnDouble("Q_max")*all_parameters.returnDouble("power_max"));
  C_.setup_x(all_parameters.returnInteger("no_of_commits"),all_parameters.returnDouble("minimum_commit"),all_parameters.returnDouble("maximum_commit"));
}

void Vary_Params2::check_grid_consistency(Fixed_Grid &X_,Fixed_Grid &Y_,Fixed_Grid &Q_,Fixed_Grid &C_)
{
  int n,m,p,q;
  double dx,dy,dq,dc;
  double m_prec=1.e-14;

  dx = all_parameters.returnDouble("X_max")/(all_parameters.returnInteger("n")-1);
  dy = all_parameters.returnDouble("Y_max")/(all_parameters.returnInteger("m")-1);
  dq = all_parameters.returnDouble("Q_max")*all_parameters.returnDouble("power_max")
  /(all_parameters.returnInteger("kmax")-1);
  dc = ( all_parameters.returnDouble("maximum_commit") - all_parameters.returnDouble("minimum_commit") )
  /(all_parameters.returnInteger("no_of_commits")-1);
  
  #ifdef DEBUG
  std::cout << " checking for grid inconsistencies...\n";
  std::cout << dx <<" X " << X_.dx() << " " << dy <<" Y " << Y_.dx() << " "  << dq <<" Q " << Q_.dx() << " "  << dc <<" C " << C_.dx() << "\n";
  #endif
  if(std::abs(dx-X_.dx())>m_prec)
  {
	n = 1 + int(all_parameters.returnDouble("X_max")/X_.dx() +0.5);
	n = max(2,n);
	all_parameters.changeParameter(n,"n");
  }
  if(std::abs(dy-Y_.dx())>m_prec)
  {
	m = 1 + int(all_parameters.returnDouble("Y_max")/Y_.dx() +0.5);
	m = max(2,m);
	all_parameters.changeParameter(m,"m");
  }
  if(std::abs(dq-Q_.dx())>m_prec)
  {
	p = 1 + int(all_parameters.returnDouble("Q_max")*all_parameters.returnDouble("power_max")/Q_.dx() +0.5);
	p = max(2,p);
	all_parameters.changeParameter(p,"kmax");
  }
  if(std::abs(dc-C_.dx())>m_prec)
  {
	q = 1 + int((all_parameters.returnDouble("maximum_commit") - all_parameters.returnDouble("minimum_commit") )/C_.dx() +0.5);
	q = max(2,q);
	all_parameters.changeParameter(q,"no_of_commits");
  }
  #ifdef DEBUG
  dx = all_parameters.returnDouble("X_max")/(all_parameters.returnInteger("n")-1);
  dy = all_parameters.returnDouble("Y_max")/(all_parameters.returnInteger("m")-1);
  dq = all_parameters.returnDouble("Q_max")*all_parameters.returnDouble("power_max")
  /(all_parameters.returnInteger("kmax")-1);
  dc = ( all_parameters.returnDouble("maximum_commit") - all_parameters.returnDouble("minimum_commit") )
  /(all_parameters.returnInteger("no_of_commits")-1);
  std::cout << " corrected grid inconsistencies...\n";
  std::cout << dx <<" X " << X_.dx() << " "  << dy <<" Y " << Y_.dx() << " "  << dq <<" Q " << Q_.dx() << " "  << dc <<" C " << C_.dx() << "\n";
  #endif
  
}

