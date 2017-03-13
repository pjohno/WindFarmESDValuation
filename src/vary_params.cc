
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "vary_params.h"
#include "fixed_grid.h"

#define DEBUG

const std::string root_dir=std::string(std::getenv("HOME"))+"/Codes/Delayed_Wind_Farm_2016/data/";
const std::string config_dir=root_dir+"Config/";
const std::string data_dir=root_dir+"Parameter_Sensitivity/";
const std::string str_default_filenames="default_filenames.in";

Vary_Params::Vary_Params(const std::string &filename)
{
  all_linking();
  input(filename);
}

void Vary_Params::link_filenames()
{
  // link parameters for iteration
  filenames.assignParameter(&grid_or_param," Iterate over the grid or the parameters ");
  filenames.assignParameter(&p1," Choose parameter for first iteration ");
  filenames.assignParameter(&p2," Choose parameter for second iteration ");
  filenames.assignParameter(&n1," Number of its (1) ");
  filenames.assignParameter(&n2," Number of its (2) ");
  filenames.assignParameter(&min1," Minimum Value (1) ");
  filenames.assignParameter(&min2," Minimum Value (2) ");
  filenames.assignParameter(&max1," Maximum value (1) ");
  filenames.assignParameter(&max2," Maximum value (2) ");
  //  link filename parameters to use
  filenames.assignParameter(&base," Base Parameters ");
  filenames.assignParameter(&store," Store Parameters ");
  filenames.assignParameter(&wind," Wind Parameters ");
  filenames.assignParameter(&price," Price Parameters ");
  filenames.assignParameter(&trading," Trading Parameters ");
  filenames.assignParameter(&grid," Grid Parameters ");
  filenames.assignParameter(&results_filename," Results Filename ");
  #ifdef DEBUG
  std::cout << " link filename parameters complete...\n";
  #endif
}

void Vary_Params::default_filenames(std::string filename)
{
  if(filenames.read_all_params(config_dir+filename))return;
  else 
  {
	std::cout << " Use default filename values...\n";
	default_filenames();
  }
}

void Vary_Params::default_filenames()
{
  grid_or_param = 0;
  p1 = 1;
  p2 = 2;
  n1 = 2;
  n2 = 2;
  min1 = 50.;
  min2 = 40.;
  max1 = 100.;
  max2 = 80.;
  base = "default_parameters.in";
  store="default_store.in";
  wind="default_wind.in";
  trading="default_trading.in";
  price="default_price.in";
  grid="default_grid.in";
  results_filename="default";
  
  std::ofstream default_file((config_dir+str_default_filenames).c_str());
  filenames.print_all_params(&default_file);
  default_file.close();
}

void Vary_Params::update()
{
  input_values(base,store,wind,price,trading,grid);
  X1.setup_x(n1,min1,max1);
  X2.setup_x(n2,min2,max2);
}

void Vary_Params::iterator()
{
  Fixed_Grid X_hat,Y_hat,Q_hat,C_hat;
  grid_consistent = 1;
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
  print_params(log_file);
  
  // iterate over parameter 1
  for(int i = 0;i<X1.size();i++)
  {
	log_file << "\n Iteration... \n";
	for(int j = 0;j<X2.size();j++)
	{
	  // change parameters
	  all_parameters.changeParameter(X1[i],p1);
	  all_parameters.changeParameter(X2[j],p2);
	  
	  if(p2==3){
	    int new_Q_Size = Q_hat.size() * ( X2[j]/Q_hat.max() );
	    all_parameters.changeParameter(new_Q_Size,3);
	  }
	  
	  // check parameters
	  log_file << "Parameters :: ";
	  all_parameters.the_reals[p1].print(&log_file);log_file << " :: ";
	  all_parameters.the_reals[p2].print(&log_file);log_file <<"\n";
	  log_file.flush();
	  if(debug_here){
		std::cout << "Parameters :: ";
		all_parameters.the_reals[p1].print(&std::cout);std::cout << " :: ";
		all_parameters.the_reals[p2].print(&std::cout);std::cout <<"\n";
		std::cout.flush();
	  }
	  if(grid_consistent)check_grid_consistency(X_hat,Y_hat,Q_hat,C_hat);
	  // write annual returns to file
	  vary_param << X1[i] << " " << X2[j] << " ";
	  vary_param << perpetual_pricer(false,false) << "\n";
	  vary_param.flush();
	}
	vary_param << "\n\n" ;
  }
  vary_param.close();
  log_file.close();
}

void Vary_Params::iterator_grid()
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
	all_parameters.the_ints[p1].print(&log_file);log_file <<"\n";
	log_file.flush();
	if(debug_here){
	  std::cout << "Parameters :: ";
	  all_parameters.the_ints[p1].print(&std::cout);std::cout <<"\n";
	  std::cout.flush();
	}
	// write annual returns to file
	vary_param << X1[i] << " ";
	vary_param << perpetual_pricer(false,false) << "\n";
  }
  vary_param.close();
  log_file.close();
  
}

void Vary_Params::input(std::string filename)
{
  default_filenames(filename);
  update();
}


void Vary_Params::update_local_grids(Fixed_Grid &X_,Fixed_Grid &Y_,Fixed_Grid &Q_,Fixed_Grid &C_)
{
  X_.setup_x(all_parameters.returnInteger(1),0.,all_parameters.returnDouble(1));
  Y_.setup_x(all_parameters.returnInteger(2),0.,all_parameters.returnDouble(2));
  Q_.setup_x(all_parameters.returnInteger(3),0.,all_parameters.returnDouble(3)*all_parameters.returnDouble(20));
  C_.setup_x(all_parameters.returnInteger(5),all_parameters.returnDouble(34),all_parameters.returnDouble(33));
}

void Vary_Params::check_grid_consistency(Fixed_Grid &X_,Fixed_Grid &Y_,Fixed_Grid &Q_,Fixed_Grid &C_)
{
  int n,m,p,q;
  double dx,dy,dq,dc;
  double m_prec=1.e-14;

  dx = all_parameters.returnDouble(1)/(all_parameters.returnInteger(1)-1);
  dy = all_parameters.returnDouble(2)/(all_parameters.returnInteger(2)-1);
  dq = all_parameters.returnDouble(3)*all_parameters.returnDouble(20)
  /(all_parameters.returnInteger(3)-1);
  dc = ( all_parameters.returnDouble(33) - all_parameters.returnDouble(34) )
  /(all_parameters.returnInteger(5)-1);
  
  #ifdef DEBUG
  std::cout << " checking for grid inconsistencies...\n";
  std::cout << dx <<" X " << X_.dx() << " " << dy <<" Y " << Y_.dx() << " "  << dq <<" Q " << Q_.dx() << " "  << dc <<" C " << C_.dx() << "\n";
  #endif
  if(std::abs(dx-X_.dx())>m_prec)
  {
	n = 1 + int(all_parameters.returnDouble(1)/X_.dx() +0.5);
	n = max(2,n);
	all_parameters.changeParameter(n,1);
  }
  if(std::abs(dy-Y_.dx())>m_prec)
  {
	m = 1 + int(all_parameters.returnDouble(2)/Y_.dx() +0.5);
	m = max(2,m);
	all_parameters.changeParameter(m,2);
  }
  if(std::abs(dq-Q_.dx())>m_prec)
  {
	p = 1 + int(all_parameters.returnDouble(3)*all_parameters.returnDouble(20)/Q_.dx() +0.5);
	p = max(2,p);
	all_parameters.changeParameter(p,3);
  }
  if(std::abs(dc-C_.dx())>m_prec)
  {
	q = 1 + int((all_parameters.returnDouble(33) - all_parameters.returnDouble(34) )/C_.dx() +0.5);
	q = max(2,q);
	all_parameters.changeParameter(q,5);
  }
  #ifdef DEBUG
  dx = all_parameters.returnDouble(1)/(all_parameters.returnInteger(1)-1);
  dy = all_parameters.returnDouble(2)/(all_parameters.returnInteger(2)-1);
  dq = all_parameters.returnDouble(3)*all_parameters.returnDouble(20)
  /(all_parameters.returnInteger(3)-1);
  dc = ( all_parameters.returnDouble(33) - all_parameters.returnDouble(34) )
  /(all_parameters.returnInteger(5)-1);
  std::cout << " corrected grid inconsistencies...\n";
  std::cout << all_parameters.returnInteger(1) << " " << dx <<" X " << " "  << all_parameters.returnInteger(2)<< " " << dy <<" Y "  << all_parameters.returnInteger(3) << " " << dq <<" Q \n";
  #endif
  
}

