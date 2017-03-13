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
using namespace std;

Optimizer::Optimizer()
{
  no_of_commits = 1;
	all_parameters.push_back(new Param_vecs());
	total_option_value.push_back(new Storage(all_parameters[0],true));
	option_vector_isnt_setup = true;
	run_starter = true;
	we_have_option_max = false;
	change_commits = false;
	option_max_old = NULL;
	option_max_new = NULL;
	option_commit = NULL;
	option_temp = NULL;
}

Optimizer::Optimizer(int no_of_options,bool progress_bar)
{
  no_of_commits = no_of_options;
  for(int i=0;i<no_of_options;i++)
  {
	all_parameters.push_back(new Param_vecs());
	total_option_value.push_back(new Storage(all_parameters[i],progress_bar));
  }
  option_vector_isnt_setup = true;
  run_starter = true;
  option_max_old = NULL;
  option_max_new = NULL;
  option_commit = NULL;
  option_temp = NULL;
  we_have_option_max = false;
  change_commits = false;
}

Optimizer::Optimizer(int no_of_options,bool progress_bar,int no_price)
{
  no_of_commits = no_of_options;
  for(int i=0;i<no_of_options;i++)
  {
	all_parameters.push_back(new Param_vecs());
	if(no_price==1)total_option_value.push_back(new No_Price(all_parameters[i],progress_bar));
	else if(no_price==2)total_option_value.push_back(new Constant_Contract(all_parameters[i],progress_bar));
	else total_option_value.push_back(new Storage(all_parameters[i],progress_bar));
  }
  option_vector_isnt_setup = true;
  run_starter = true;
  option_max_old = NULL;
  option_max_new = NULL;
  option_commit = NULL;
  option_temp = NULL;
  we_have_option_max = false;
  change_commits = false;
}

Optimizer::~Optimizer()
{
	bool debug=false;
	delete option_max_old;
	delete option_max_new;
	delete option_commit;
	delete option_temp;
	for(int i=0;i<total_option_value.size();i++)
	{
			delete all_parameters[i];
			delete total_option_value[i];
			if(debug)std::cout << "\n Deconstructing... <" << i << ">\n";
	}
}

void Optimizer::resize_total_option(int no_of_options)
{

  bool debug=false;
  if( no_of_options == no_of_commits)return;
  if(debug)std::cout << "Resizing total option...\n First delete " << no_of_commits << " options";
  for(int i=0;i<no_of_commits;i++)
  {
	delete all_parameters[i];
	delete total_option_value[i];
	if(debug)std::cout << "\n Deconstructing... <" << i << ">\n";
  }
  all_parameters.clear();
  total_option_value.clear();
  if(debug)std::cout << " Set new number of commits...\n";
  no_of_commits = no_of_options;
  if(debug)std::cout << " Push back new options...\n";
  for(int i=0;i<no_of_options;i++)
  {
	all_parameters.push_back(new Param_vecs());
	total_option_value.push_back(new Storage(all_parameters[i],false));
	set_base_parameters(i);
  }
  
}

void Optimizer::set_base_parameters(int i)
{
// Set up base parameters...
	
	all_parameters[i]->changeParameter(101,1); // change nodes in X	
	all_parameters[i]->changeParameter(11,2); // change nodes in Y
	all_parameters[i]->changeParameter(11,3); // change nodes in Q
	all_parameters[i]->changeParameter(0.5,0); // change commit level - percentage of power output
	all_parameters[i]->changeParameter(100.,1); // change X max
	all_parameters[i]->changeParameter(100.,2); // change Y max
	all_parameters[i]->changeParameter(1.,3); // change Q max - percentage of power output
	all_parameters[i]->changeParameter(1.,4); // change maturity
	all_parameters[i]->changeParameter(0.05,5); // change interest rate
	all_parameters[i]->changeParameter(0.1,6); // change mean rev in X
	all_parameters[i]->changeParameter(8.,7); // change average in X
	all_parameters[i]->changeParameter(0.2,8); // change variance in X
	all_parameters[i]->changeParameter(0.04,9); // change mean rev in Y
	all_parameters[i]->changeParameter(40.0,10); // change average in Y
	all_parameters[i]->changeParameter(0.075,11); // change variance in Y
	all_parameters[i]->changeParameter(0.8366,12); // change energy efficiency
	all_parameters[i]->changeParameter(0.25,13); // change charge - percentage of power output
	all_parameters[i]->changeParameter(0.25,14); // change discharge - percentage of power output
	all_parameters[i]->changeParameter(0.2,15); // change smoothing on charge
	all_parameters[i]->changeParameter(0.2,16); // change smoothing on discharge
	all_parameters[i]->changeParameter(0.9,17); // change contract sell price 
	all_parameters[i]->changeParameter(1.1,18); // change contract bid price
	all_parameters[i]->changeParameter(25.,19); // change power cut out
	all_parameters[i]->changeParameter(1.,20); // maximum power output 
	all_parameters[i]->changeParameter(10.,21); // change power turn
	all_parameters[i]->changeParameter(1.,22); // change contract length
	all_parameters[i]->changeParameter(0.375,27); // amp of X
	all_parameters[i]->changeParameter(-2.,29); // phase of X
	all_parameters[i]->changeParameter(0.375,30); // amp of Y
	all_parameters[i]->changeParameter(-14.,32); // phase of Y
	all_parameters[i]->changeParameter(1.,33); // maximum commit
	all_parameters[i]->changeParameter(0.01,34); // minimum commit
	all_parameters[i]->changeParameter(1.,35); // buy sell spread
	
	
	
	test();
	// -----------------------------------------------------------------------
}

void Optimizer::vary_parameter(int parameter_ref,double low,double high,int its,string filename)
{

	int j_value,q_max;
	double param,dx,dy,d_its;
	vector<double> temp_empty,temp_full;
	ofstream vary_param,file_param;

	file_param.open("parameters_base.dat");
	set_base_parameters(0);	

	dx = all_parameters[0]->returnDouble(1) / (double(all_parameters[0]->returnInteger(1)) - 1.);
	dy = all_parameters[0]->returnDouble(2) / (double(all_parameters[0]->returnInteger(2)) - 1.);
	j_value = int( 0.5*all_parameters[0]->returnDouble(2) / dy + 0.5);
//	j_value = int( 40. / dy + 0.5);
	q_max = all_parameters[0]->returnInteger(3) - 1;

//	cout << j_value << " " << q_max;

	vary_param.open(filename.c_str());
	vary_param << setprecision(8);
	
	if(its==1){d_its = (high - low)/2.;}
	else{(high - low)/double(its-1);}
	
	for(int i = 0;i<its;i++)
	{	
		cout << "\n Iteration... \n";
		param = double(i)*d_its + low;
		cout << "Parameter :: "<<parameter_ref << " Value :: " << param << endl;
		all_parameters[0]->changeParameter(param,parameter_ref); // change parameter...
		all_parameters[0]->print_all_params(&file_param);
		total_option_value[0]->runSolver();
		total_option_value[0]->return_vector_u(j_value,0,&temp_empty);
		total_option_value[0]->return_vector_u(j_value,q_max,&temp_full);
		for(int ii=0;ii<temp_empty.size();ii++)
		{
			vary_param.setf(ios::left,ios::adjustfield);
			vary_param << ii*dx << " " << param << " " << temp_empty[ii] ;
			vary_param << " " << temp_full[ii] << endl;
		}
		vary_param << "\n\n";
	}
	vary_param.close();
	file_param.close();
	
}

void Optimizer::vary_two_parameters(int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,string filename)
{

	bool file_setup=true;

	int j_value,q_max;
	double param_1,param_2,dx,dy;
	vector<double> temp_empty,temp_full;
	ofstream vary_param,log_file;
	// Set base parameters
	set_base_parameters(0);	
	// set up variables for printing... 
	dx = all_parameters[0]->returnDouble(1) / (double(all_parameters[0]->returnInteger(1)) - 1.);
	dy = all_parameters[0]->returnDouble(2) / (double(all_parameters[0]->returnInteger(2)) - 1.);
	j_value = int( 0.5*all_parameters[0]->returnDouble(2) / dy + 0.5);
	q_max = all_parameters[0]->returnInteger(3) - 1;
	// Open file for printing
	vary_param.open(filename.c_str());
	log_file.open("log_option_pricer.log");
	vary_param << setprecision(8);

	vector<int> position;
					
	for(int i = 0;i<its_1;i++)
	{

		// iterate over parameter 1
		param_1 = (high_1 - low_1)*double(i)/double(its_1-1) + low_1;
		all_parameters[0]->changeParameter(param_1,parameter_ref_1); // change parameter...

		for(int j = 0;j<its_2;j++)
		{

			log_file << "\n Iteration... \n";
			// iterate over parameter 2
			param_2 = (high_2 - low_2)*double(j)/double(its_2-1) + low_2;
			all_parameters[0]->changeParameter(param_2,parameter_ref_2); // change parameter...
	
			// check parameters 
			log_file << "Parameter :: "<<parameter_ref_1 << " Value :: " << param_1;
			log_file << " Parameter :: "<<parameter_ref_2 << " Value :: " << param_2 << endl;
			log_file.flush();			
				
			// run solver with zero initial conditions
			total_option_value[0]->runSolver();
			// return option value for specific j value with both empty and full tank
			total_option_value[0]->return_vector_u(j_value,0,&temp_empty);
			total_option_value[0]->return_vector_u(j_value,q_max,&temp_full);
			log_file << " Solved ";			
			log_file.flush();			
			
			// setup write to file...
			if(file_setup){
			for(int ii=0;ii<temp_empty.size();ii++){
				vary_param.setf(ios::left,ios::adjustfield);
				vary_param.width(20);
				vary_param << ii*dx;
				vary_param.width(20);
				vary_param << param_1 ;  
				vary_param.width(20);
				vary_param << param_2 ;
				vary_param.width(20);
				vary_param << temp_empty[ii] ;
				vary_param.width(20);
				vary_param << temp_full[ii] << endl;
				position.push_back(vary_param.tellp());
				for(int iii = 0;iii<its_1;iii++){	for(int jjj = 0;jjj<its_2;jjj++){
					vary_param.width(100);
					vary_param << " " << endl;
					}
				vary_param.width(100);
				vary_param << " " << endl;
				}
				file_setup=false;
			}}else{
			for(int ii=temp_empty.size()-1;ii>=0;ii--)
			{
				vary_param.seekp(position[ii]);
				vary_param.setf(ios::left,ios::adjustfield);
				vary_param.width(20);
				vary_param << ii*dx ;
				vary_param.width(20);
				vary_param << param_1 ;  
				vary_param.width(20);
				vary_param << param_2 ;
				vary_param.width(20);
				vary_param << temp_empty[ii] ;
				vary_param.width(20);
				vary_param << temp_full[ii] << endl;
				position[ii] = vary_param.tellp();
			}}	
		}
		for(int ii =temp_empty.size()-1;ii>=0;ii--)
		{
			vary_param.seekp(position[ii]);
			vary_param.width(100);
			vary_param << " " << endl;
			position[ii] = vary_param.tellp();
		}
	}
	vary_param.close();
	
}


void Optimizer::vary_two_parameters_multi(int no_of_contracts,int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,string filename)
{
	
	bool file_setup=true;
	
	int j_value,q_max;
	double param_1,param_2,dx,dy;
	vector<double> temp_empty,temp_full;
	ofstream vary_param,log_file;
	// Set base parameters
	for(int ii=0;ii<total_option_value.size();ii++){
	set_base_parameters(ii);	}
	// set up variables for printing... 
	dx = all_parameters[0]->returnDouble(1) / (double(all_parameters[0]->returnInteger(1)) - 1.);
	dy = all_parameters[0]->returnDouble(2) / (double(all_parameters[0]->returnInteger(2)) - 1.);
	j_value = int( all_parameters[0]->returnDouble(10) / dy + 0.5);
	q_max = all_parameters[0]->returnInteger(3) - 1;
	// Open file for printing
	vary_param.open(filename.c_str());
	log_file.open("log_option_pricer.log");
	vary_param << setprecision(8);
	
	vector<int> position;
	
	for(int i = 0;i<its_1;i++)
	{
		
		// iterate over parameter 1
		param_1 = (high_1 - low_1)*double(i)/double(its_1-1) + low_1;
		for(int j = 0;j<its_2;j++)
		{
			
			log_file << "\n Iteration... \n";
			// iterate over parameter 2
			param_2 = (high_2 - low_2)*double(j)/double(its_2-1) + low_2;
			for(int ii=0;ii<total_option_value.size();ii++){
				all_parameters[ii]->changeParameter(param_1,parameter_ref_1); // change parameter...
				all_parameters[ii]->changeParameter(param_2,parameter_ref_2); // change parameter...
			}
			// check parameters 
			log_file << "Parameter :: "<<parameter_ref_1 << " Value :: " << param_1;
			log_file << " Parameter :: "<<parameter_ref_2 << " Value :: " << param_2 << endl;
			log_file.flush();			
			
			// run solver with zero initial conditions
			solve_multiple_time(no_of_contracts,string("test.dat"),false,false);
			// return option value for specific j value with both empty and full tank
			total_option_value[0]->return_vector_u(j_value,0,&temp_empty);
			total_option_value[0]->return_vector_u(j_value,q_max,&temp_full);
			log_file << " Solved ";			
			log_file.flush();			
			
			// setup write to file...
			if(file_setup){
				for(int ii=0;ii<temp_empty.size();ii++){
					vary_param.setf(ios::left,ios::adjustfield);
					vary_param.width(20);
					vary_param << ii*dx;
					vary_param.width(20);
					vary_param << param_1 ;  
					vary_param.width(20);
					vary_param << param_2 ;
					vary_param.width(20);
					vary_param << temp_empty[ii] ;
					vary_param.width(20);
					vary_param << temp_full[ii] << endl;
					position.push_back(vary_param.tellp());
					for(int iii = 0;iii<its_1;iii++){	for(int jjj = 0;jjj<its_2;jjj++){
						vary_param.width(100);
						vary_param << " " << endl;
					}
					vary_param.width(100);
					vary_param << " " << endl;
					}
					file_setup=false;
				}}else{
					for(int ii=temp_empty.size()-1;ii>=0;ii--)
					{
						vary_param.seekp(position[ii]);
						vary_param.setf(ios::left,ios::adjustfield);
						vary_param.width(20);
						vary_param << ii*dx ;
						vary_param.width(20);
						vary_param << param_1 ;  
						vary_param.width(20);
						vary_param << param_2 ;
						vary_param.width(20);
						vary_param << temp_empty[ii] ;
						vary_param.width(20);
						vary_param << temp_full[ii] << endl;
						position[ii] = vary_param.tellp();
					}}	
		}
		for(int ii =temp_empty.size()-1;ii>=0;ii--)
		{
			vary_param.seekp(position[ii]);
			vary_param.width(100);
			vary_param << " " << endl;
			position[ii] = vary_param.tellp();
		}
	}
	vary_param.close();
	
}

void Optimizer::vary_tp_multi_new(int no_of_contracts,int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,string filename)
{
	
	bool file_setup=true,debug=false;
	int temp_int;
	double param_1,param_2,x_value,x_value_star,y_value,q_empty,q_full;

	ofstream vary_param,log_file;
	// Set base parameters
	for(int ii=0;ii<total_option_value.size();ii++)
	{
	set_base_parameters(ii);	
	}
	x_value = all_parameters[0]->returnDouble(7);
	y_value = all_parameters[0]->returnDouble(10);
	q_empty = 0.;
	q_full = all_parameters[0]->returnDouble(3);
	// Open file for printing
	vary_param.open(filename.c_str());
	log_file.open("log_option_pricer.log");
	vary_param << setprecision(8);
	
	
	for(int i = 0;i<its_1;i++)
	{
		
		// iterate over parameter 1
		param_1 = (high_1 - low_1)*double(i)/double(its_1-1) + low_1;
		for(int j = 0;j<its_2;j++)
		{
			
			log_file << "\n Iteration... \n";
			// iterate over parameter 2
			param_2 = (high_2 - low_2)*double(j)/double(its_2-1) + low_2;
			for(int ii=0;ii<total_option_value.size();ii++){
				all_parameters[ii]->changeParameter(param_1,parameter_ref_1); // change parameter...
				all_parameters[ii]->changeParameter(param_2,parameter_ref_2); // change parameter...
			}
			// check parameters 
			log_file << "Parameter :: "<<parameter_ref_1 << " Value :: " << param_1;
			log_file << " Parameter :: "<<parameter_ref_2 << " Value :: " << param_2 << endl;
			log_file.flush();			
			
			// run solver with zero initial conditions
			solve_multiple_time(no_of_contracts,string("test.dat"),false,false);
			log_file << " Solved ";			
			log_file.flush();			
			
			if(debug){print_results(string("test.dat"),debug);cin >> temp_int;}
			
			// setup write to file...
			x_value_star = total_option_value[0]->average_X_t(-no_of_contracts*all_parameters[0]->returnDouble(22));
			vary_param.setf(ios::left,ios::adjustfield);
			vary_param << param_1 << " " << param_2 << " ";
			vary_param << option_value_3(x_value,y_value,q_empty,4) << " " ;
			vary_param << option_value_3(x_value,y_value,q_full,4) << " " ;
			vary_param << option_value_3(x_value_star,y_value,q_empty,4) << " " ;
			vary_param << option_value_3(x_value_star,y_value,q_full,4) << " \n" ;
		}
		vary_param << " \n" ;
	}
	vary_param.close();
	
}

void Optimizer::vary_perpetual(int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,string filename)
{
	
	bool file_setup=true,debug=false;
	int temp_int;
	double param_1,param_2,x_value,x_value_star,y_value,y_value_star,q_empty,q_full,dX,d_its_1;
	vector<double> grid_size(3);
	ofstream vary_param,log_file;
	string results_file;
	// Set base parameters
	for(int ii=0;ii<total_option_value.size();ii++)
	{
		set_base_parameters(ii);	
	}
	x_value = all_parameters[0]->returnDouble(7);
	y_value = all_parameters[0]->returnDouble(10);
	q_empty = 0.;
	q_full = all_parameters[0]->returnDouble(3);
	if(parameter_ref_1<4 and parameter_ref_1>0)
	{
		grid_size[parameter_ref_1-1] = all_parameters[0]->returnDouble(parameter_ref_1);
		dX = grid_size[parameter_ref_1-1]/
		double(all_parameters[0]->returnInteger(parameter_ref_1)-1);
	}
	if(parameter_ref_2<4 and parameter_ref_2>0)
	{
		grid_size[parameter_ref_2-1] = all_parameters[0]->returnDouble(parameter_ref_2);
		dX = grid_size[parameter_ref_2-1]/
		double(all_parameters[0]->returnInteger(parameter_ref_2)-1);
	}
	// Open file for printing
	results_file = "vary_param."+filename+".dat";
	vary_param.open(results_file.c_str());
	results_file = filename+".log";
	log_file.open(results_file.c_str());
	vary_param << setprecision(8);
	
	if(its_1==1){d_its_1 = (high_1-low_1)/2.;}
	else{d_its_1 = (high_1-low_1)/double(its_1-1);}
	double max_commit,dC;
	if(parameter_ref_2 == 33){
	  max_commit = all_parameters[0]->returnDouble(33);
	  dC = max_commit / (no_of_commits-1);
	  if(debug){cout << "No of commits.. " << no_of_commits << " with dC = " << dC <<"\n";}
	}
	  
	for(int i = 0;i<its_1;i++)
	{
		
		// iterate over parameter 1
		param_1 = d_its_1*double(i) + low_1;
		for(int j = 0;j<its_2;j++)
		{
			
			log_file << "\n Iteration... \n";
			// iterate over parameter 2
			param_2 = (high_2 - low_2)*double(j)/double(its_2-1) + low_2;
			// If we are changing the maximum commit try to alter number of options to compensate
			if(parameter_ref_2==33){
			  int temp;
			  temp = int(param_2/dC + 1.5);
			  if(debug){cout << "No of commits.. " << temp << " with dC = " << dC <<"\n";}
			  resize_total_option(temp);
			}
			
			for(int ii=0;ii<total_option_value.size();ii++){
				all_parameters[ii]->changeParameter(param_1,parameter_ref_1); // change parameter...
				all_parameters[ii]->changeParameter(param_2,parameter_ref_2); // change parameter...
			}
			// check parameters 
			log_file << "Parameters :: ";
			all_parameters[0]->the_reals[parameter_ref_1].print(&log_file);
			log_file << " :: ";
			all_parameters[0]->the_reals[parameter_ref_2].print(&log_file);
			log_file << endl;
			log_file.flush();			
			if(debug){cout <<"Parameters :: ";
			all_parameters[0]->the_reals[parameter_ref_1].print(&cout);
			cout << " :: ";
			all_parameters[0]->the_reals[parameter_ref_2].print(&cout);
			cout << endl;
			}
			if(parameter_ref_1<4 and parameter_ref_1>0)
			{
				grid_size[parameter_ref_1-1] = all_parameters[0]->returnDouble(parameter_ref_1);
				temp_int = int(grid_size[parameter_ref_1-1]/dX + 1.5);
				for(int ii=0;ii<total_option_value.size();ii++){
					all_parameters[ii]->changeParameter(temp_int,parameter_ref_1); // change parameter...
				}
			}
			if(parameter_ref_2<4 and parameter_ref_2>0)
			{
				grid_size[parameter_ref_2-1] = all_parameters[0]->returnDouble(parameter_ref_2);
				temp_int = int(grid_size[parameter_ref_2-1]/dX + 1.5);
				for(int ii=0;ii<total_option_value.size();ii++){
					all_parameters[ii]->changeParameter(temp_int,parameter_ref_2); // change parameter...
				}
			}
			// run solver with zero initial conditions
			perpetual_pricer(false);
			log_file << " Solved \n";			
			log_file.flush();			
			
			if(debug){print_results(string("test.dat"),debug);cin >> temp_int;}
			
			// setup write to file...
			x_value_star = total_option_value[0]->average_X_t(0.);
			y_value_star = total_option_value[0]->average_Y_t(0.);
			q_full = all_parameters[0]->returnDouble(3);
			vary_param.setf(ios::left,ios::adjustfield);
			vary_param << grid_i << " " << grid_j << " " << grid_k << " "<< no_of_commits << " ";
			vary_param << all_parameters[0]->returnInteger(4) << " " << param_1 << " " << param_2 << " ";
			vary_param << option_value_3(x_value,y_value,q_empty,4) << " " ;
			vary_param << option_value_3(x_value,y_value,q_full,4) << " " ;
			vary_param << option_value_3(x_value_star,y_value_star,q_empty,4) << " " ;
			vary_param << option_value_3(x_value_star,y_value_star,q_full,4) << " \n" ;
		}
		vary_param << " \n" ;
		vary_param.flush();
	}
	vary_param.close();
	
}

void Optimizer::grid_test_X_max(void)
{
	
	int ii;
	double i_max;
	
	stringstream str_stream;
	string temp;
	ofstream x_max;
	
	x_max.open("x_max_data.dat");
	all_parameters[0]->changeParameter(21,2);
	all_parameters[0]->changeParameter(21,3);
	for(int grid = 0;grid < 201;grid+=20)
	{
		str_stream.str("");
		str_stream << grid;
		temp = "test_grid." + str_stream.str() + ".dat";
		ii = 200+ grid + 1;
		i_max = 50. + grid/4;
		all_parameters[0]->changeParameter(ii,1); // change nodes in X
		all_parameters[0]->changeParameter(i_max,1); // change X_max
		total_option_value[0]->runSolver();
// 		total_option_value[0]->print_results(temp);
// 		x_max << "Value :: " << total_option_value[0]->option_value_3(25.,1.,1.,4);
// 		x_max << total_option_value[0]->option_value_3(45.,1.,1.,4) << endl;
	}

	x_max.close();

}

void Optimizer::grid_test(void)
{
	
	int ii,jj,kk;
	
	ofstream grid_data,grid_plot;
	
	// set base parameters
	set_base_parameters(0);
	all_parameters[0]->changeParameter(24.,4); // change maturity
	
	
	// get values in the middle of the grid  
	double x_value,y_value,q_value;
	
	grid_data.open("grid_data.dat");
	grid_plot.open("grid_plot.dat");
	
	for(int grid = 0;grid < 501;grid+=50)
	{
		cout << grid <<endl;
		// change grid size
		ii = 100+ grid + 1;
		jj = 10+grid/10+ 1;
		kk = 10+grid/10 + 1;
		all_parameters[0]->changeParameter(ii,1);
		all_parameters[0]->changeParameter(jj,2);
		all_parameters[0]->changeParameter(kk,3);
		// generate values in the center of the grid
		x_value = 0.5*all_parameters[0]->returnDouble(1);
		y_value = 0.5*all_parameters[0]->returnDouble(2);
		q_value = all_parameters[0]->returnDouble(3);
		
		// run solver with new grid
		total_option_value[0]->runSolver();

		// Output grid data to file
		grid_data << " i :: " << ii << " j :: " << jj << " k :: " << kk << " Value :: " ;
		grid_data << total_option_value[0]->return_option_new()(x_value,y_value,0.) << " " << total_option_value[0]->return_option_new()(x_value,y_value,q_value)<< endl;
		grid_plot.setf(ios::scientific,ios::floatfield);
		grid_plot << setprecision(8);
		grid_plot << ii << " " << total_option_value[0]->return_option_new()(x_value,y_value,0.) << " " << total_option_value[0]->return_option_new()(x_value,y_value,q_value)<< endl;
	}
	
	grid_data.close();
}

void Optimizer::grid_test_dt(string filename,int grid_start,int grid_finish,int grid_incr,int param)
{
	
	int ii,jj,kk,ll;
	int no_of_its;
	double temp,dtemp,temp1,temp2,value1,value2;
	string data,plot;

	ofstream grid_data,grid_plot;
	
	data = "grid_data." + filename + ".dat";
	plot = "grid_plot." + filename + ".dat";
	grid_data.open(data.c_str());
	grid_plot.open(plot.c_str());
	
		ii = 101;
		jj = 101;
		kk = 26;
		ll = 100;
		all_parameters[0]->changeParameter(ii,1);
		all_parameters[0]->changeParameter(jj,2);
		all_parameters[0]->changeParameter(kk,3);
		all_parameters[0]->changeParameter(ll,4);

	no_of_its = (grid_finish - grid_start) / grid_incr;
	
	temp = 1./grid_finish;
	dtemp = (1./grid_start - 1./grid_finish)/no_of_its;
	cout << dtemp << " " << no_of_its << endl;
	for(int grid = grid_start;grid <= grid_finish;grid+= grid_incr)
	{
		
		ii = int(1./temp + 0.5);
		cout << ii << " " << temp << endl;
		temp = temp + dtemp;
		
//		all_parameters[0]->changeParameter(ii,param);
		all_parameters[0]->changeParameter(10*ii+1,1);
		all_parameters[0]->changeParameter(2*ii+1,2);
		all_parameters[0]->changeParameter(ii+1,3);
		all_parameters[0]->changeParameter(ii,4);
		total_option_value[0]->runSolver();

		value1 = total_option_value[0]->return_option_new()(5.,40.,2.5);
		value2 = total_option_value[0]->return_option_new()(15.,40.,2.5);

		grid_data << " i :: " <<  all_parameters[0]->returnInteger(1);
		grid_data << " j :: " << all_parameters[0]->returnInteger(2);
		grid_data << " k :: " << all_parameters[0]->returnInteger(3);
    grid_data << " l :: " << all_parameters[0]->returnInteger(4);
		grid_data << " Value :: " ;
		grid_data << value1 << " " << value2<< endl;
		grid_plot.setf(ios::scientific,ios::floatfield);
		if(grid==grid_start){temp1=value1;temp2=value2;}
		grid_plot << setprecision(8);
		grid_plot << ii << " " << value1 << " " << value2 << 
		" " << temp1 << " " << temp2 << endl;
	}
	
	grid_data.close();
	grid_plot.close();
}

void Optimizer::test(void)
{

	ofstream params;
	params.open("parameters.dat");
	total_option_value[0]->test_initial_conditions();
	all_parameters[0]->print_all_params(&params);
	params.close();
}

void Optimizer::test_sin_curves(void)
{
  double pi=3.141592653589793238462643383279502884197,tau;
	
  ofstream sin_graph;
  sin_graph.open("sin_graph.dat");
  for(int i=0;i<48;i++){
	tau = -i/2.;
	sin_graph << tau << " " << total_option_value[0]->average_X_t(tau)
	<< " " << total_option_value[0]->average_Y_t(tau) << endl;
  }

}

void Optimizer::solve_one_option(void)
{
	set_base_parameters(0);
	//all_parameters[0]->changeParameter(5,2);
	all_parameters[0]->changeParameter(24.,4);
	total_option_value[0]->runSolver();
	total_option_value[0]->print_results(string("results.test.x.dat"));
	total_option_value[0]->print_results_old(string("results.test.y.dat"));
}

void Optimizer::solve_for_multiple_option(int no_of_hours,string filename)
{
	for(int ii=0;ii<total_option_value.size();ii++)
	{
		set_base_parameters(ii);	
	}
	
	solve_multiple_time(no_of_hours,filename,true,false);
	
}

void Optimizer::solve_multiple_time(int no_of_hours,string filename,bool display,bool reset)
{	
	if(display)
	{
	std::cout << " Solving for a multi-time option.\n No of commits :: " 
	<< total_option_value.size();
	std::cout << "\n No of contracts :: " << no_of_hours << std::endl;
	}
	
	// reset initial conditions
	we_have_option_max = reset;
	
	double tau;
	int current_commits;
	stringstream str_stream;
	string temp;
	char buffer[43] = { 0 };
	char percent[8] = "0.0%%";
	int ii,jj,i_step;
	buffer[0] = '[';i_step = 0;
	
	for(int i=0;i<no_of_hours;i++)
	{
		// Write progress to screen...
		ii = int(double(i)/double(no_of_hours) * 40.);jj = i % 4;	
		if(jj == 0)buffer[ii + 1] = '\\';
		else if(jj == 1)buffer[ii + 1] = '|';
		else if(jj == 2)buffer[ii + 1] = '/';
		else buffer[ii + 1] = '-';
		
		if(display)
		{
			if(i_step<ii){for(int k=1;k<ii+1;k++)buffer[k] = '-';i_step = ii;}
			for(int k = ii + 2; k < 41; k++)buffer[k] = ' ';buffer[41] = ']';
			sprintf(percent, "%3.2f%%", double(i) / double(no_of_hours) * 100.0);
			printf("%s%s\r", buffer, percent);
			cout.flush();
		}
		
		// progress written
		str_stream.str("");
		str_stream << no_of_hours - i - 1 ;
		temp = filename + str_stream.str() + ".dat";
		tau = -i*all_parameters[0]->returnDouble(22); // change contract length
		// solve for all commits
		solve_for_all_C(tau);
		// calculate number of commits
		change_commits=false;
		if(change_commits)
		{
		  current_commits = force_commit(tau);
		}
		else
		{
			current_commits = total_option_value.size();
		}
		// fin option max over commits
		find_option_max(current_commits);
		if(display)
		{
			// print option value
			//if(i%4==3)
			print_results(temp,false);
		}
		
	}

	//str_stream.str("");
	//str_stream << total_option_value.size();
	//temp = "test." + str_stream.str() + ".dat";
// 	if(display)
// 	{
// 		temp = filename + ".dat";
// 		print_results(temp,false);
// 	}
// 	
	if(display)
	{
			for(int k = 1; k < 41; k++)buffer[k] = '-';
		buffer[41] = ']';
		sprintf(percent, "%3.2f%%", 100.0);
		printf("%s%s", buffer, percent);
	
		std::cout << " \n Option solved.\n";
	}

}

void Optimizer::solve_for_all_C(double tau)
{
	bool debug=false;
	double commit_rate;
	
	if(we_have_option_max)
	{
		for(int i=0;i<total_option_value.size();i++)
		{
			total_option_value[i]->runSolver(option_max_old,option_max_new,tau);
		}
	}
	else
	{
	  double commit_max,commit_min;
	  commit_max = all_parameters[0]->returnDouble(33);
	  commit_min = all_parameters[0]->returnDouble(34);
	  for(int i=0;i<total_option_value.size();i++)
		{	
		  if(total_option_value.size()==1)
		  {
			commit_rate = commit_min + 0.5*(commit_max - commit_min);
		  }
		  else
		  {
			commit_rate = max(commit_min,i*commit_max/(total_option_value.size()-1));
		  }
		  if(debug)std::cout << " Commit rate :: " << commit_rate << "\n";
		  all_parameters[i]->changeParameter(commit_rate,0);
		  all_parameters[i]->changeParameter(all_parameters[0]->returnDouble(22),4);
		  // 			all_parameters[i]->print_all_params(&cout);
		  total_option_value[i]->runSolver();
		}
		setup_grid();
	  we_have_option_max = true;
	}
	
}

void Optimizer::find_option_max(int current_commits)
{
	
	vector<double> temp;
	if(option_vector_isnt_setup)
	{
		setup_grid();
	}
	
	for(int ii=0;ii<total_option_value.size();ii++)
	{
		if(grid_i != all_parameters[ii]->returnInteger(1)){cout << "\nConflicting parameters... Optimizer::find_option_max\n";return;}
		if(grid_j != all_parameters[ii]->returnInteger(2)){cout << "\nConflicting parameters... Optimizer::find_option_max\n";return;}
		if(grid_k != all_parameters[ii]->returnInteger(3)){cout << "\nConflicting parameters... Optimizer::find_option_max\n";return;}
	}
	
	for(int j = 0;j<grid_j;j++)
	{
		for(int k = 0;k<grid_k;k++)
		{
			// return  Option new values
			total_option_value[0]->return_vector_u(j,k,&temp);
			for(int i=0;i<grid_i;i++)
			{
				(*option_max_new)[i][j][k] = temp[i];
				(*option_commit)[i][j][k] = all_parameters[0]->returnDouble(0);
			}
			for(int ii=1;ii<current_commits;ii++)
			{
				total_option_value[ii]->return_vector_u(j,k,&temp);
				for(int i=0;i<grid_i;i++)
				{
					if((*option_max_new)[i][j][k] < temp[i])
					{
						(*option_max_new)[i][j][k] = temp[i];
						(*option_commit)[i][j][k] = all_parameters[ii]->returnDouble(0);
					}
				}
			}
			// return Option old values
			total_option_value[0]->return_vector_uold(j,k,&temp);
			for(int i=0;i<grid_i;i++)
			{
				(*option_max_old)[i][j][k] = temp[i];
			}
			for(int ii=1;ii<current_commits;ii++)
			{
				total_option_value[ii]->return_vector_uold(j,k,&temp);
				for(int i=0;i<grid_i;i++)
				{
					(*option_max_old)[i][j][k] = max(temp[i],(*option_max_old)[i][j][k]);
				}
			}
		}
	}
	
}

void Optimizer::setup_grid(void)
{
	vector<double> temp;
	
	grid_i =  all_parameters[0]->returnInteger(1);
	grid_j = all_parameters[0]->returnInteger(2);
	grid_k = all_parameters[0]->returnInteger(3);
	dx = all_parameters[0]->returnDouble(1)/double(grid_i-1);
	dy = all_parameters[0]->returnDouble(2)/double(grid_j-1);
	dz = all_parameters[0]->returnDouble(3)*all_parameters[0]->returnDouble(20)/double(grid_k-1);
	
	delete option_max_old;
	delete option_max_new;
	delete option_commit;
	delete option_temp;
	
	temp.resize(grid_k);
	vector<vector<double> > temp_j(grid_j,temp);
	option_max_old =  new vector<vector<vector<double> > > (grid_i,temp_j);
	option_max_new =  new vector<vector<vector<double> > > (grid_i,temp_j);
	option_commit =  new vector<vector<vector<double> > > (grid_i,temp_j);
	option_temp =  new vector<vector<vector<double> > > (grid_i,temp_j);
	option_vector_isnt_setup=false;
	
}

void Optimizer::print_results(string filename,bool debug)
{
	int Width=20;
	ofstream graph;
	graph.open(filename.c_str());
	// check if the file is open or not
	if (graph.is_open() and !option_vector_isnt_setup)
	{
		if(debug)cout << "Print to file... " << endl;
		// write out stock vs option value
		for(int j = 0;j<grid_j;j = j + max(1,grid_j/50))
		{
			for(int i=0;i<grid_i;i = i + max(1,grid_i/50))
			{
				graph.setf(ios::scientific,ios::floatfield);
				graph << setprecision(8);
				graph << i*dx << " " << j*dy << " " 
					 << (*option_max_new)[i][j][0] << " "
					 << (*option_max_new)[i][j][grid_k-1]<< " "
					 << (*option_commit)[i][j][0]<< " "
					 << (*option_commit)[i][j][grid_k-1]<< endl;
			}
			graph << "\n";
		}
	}
	graph.close();
	if(debug)cout << "Printing Finished\n";
}

double Optimizer:: option_value_1(double value_stock_1,int j,int k,int accuracy)
{
	
	// 	cout << "XX :: " << endl;
	
	Methods option_value;
	
	vector<double> yy,gg;
	
	int option_int;
	double half_temp;
	half_temp = double(accuracy)/2. - double(accuracy/2);
	option_int = int( value_stock_1 / dx  + half_temp) - (accuracy-1)/2;
	option_int = max(0,min(grid_i-accuracy,option_int));
	for(int i=option_int;i<option_int+accuracy;i++){
		yy.push_back(i*dx);
		gg.push_back((*option_max_new)[i][j][k]);
		// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
	}
	return option_value.interpolate(value_stock_1,yy,gg);
	
}

double Optimizer:: option_value_2(double value_stock_1,double value_stock_2,int k,int accuracy)
{
	
	// 	cout << "YY :: " << endl;
	
	Methods option_value;
	
	vector<double> yy,gg;
	
	int option_int;
	double half_temp;
	half_temp = double(accuracy)/2. - double(accuracy/2);
	option_int = int( value_stock_2 / dy + half_temp) - (accuracy-1)/2;
	option_int = max(0,min(grid_j-accuracy,option_int));
	for(int j=option_int;j<option_int+accuracy;j++){
		yy.push_back(j*dy);
		gg.push_back(option_value_1(value_stock_1,j,k,accuracy));
		// 		    cout <<j << " " << yy[j-option_int] << " " << gg[j-option_int] << endl;
	}
	return option_value.interpolate(value_stock_2,yy,gg);
	
}

double Optimizer:: option_value_3(double value_stock_1,double value_stock_2,double value_stock_3,int accuracy)
{
	if(not we_have_option_max)return 0.;
	// 	cout << "QQ :: "<<accuracy << endl;
	
	Methods option_value;
	
	vector<double> yy,gg;
	
	int option_int;
	double half_temp;
	half_temp = double(accuracy)/2. - double(accuracy/2);
	option_int = int( value_stock_3 / dz+ half_temp) - (accuracy-1)/2;
	option_int = max(0,min(grid_k-accuracy,option_int));
	for(int k=option_int;k<option_int+accuracy;k++){
		yy.push_back(k*dz);
		gg.push_back(option_value_2(value_stock_1,value_stock_2,k,accuracy));
		// 		    cout << k << " " << yy[k-option_int] << " " << gg[k-option_int] << endl;
	}
	return option_value.interpolate(value_stock_3,yy,gg);
	
}

double Optimizer::option_average()
{
	if(not we_have_option_max)return 0.;
	
	double temp=0.;
	
	for(int i=0;i<grid_i;i++){
		for(int j=0;j<grid_j;j++){
			for(int k=0;k<grid_k;k++){
				// 				cout << i << " " << j << " " << k << " " << (*option_max_new)[i][j][k] << endl;
				temp+=(*option_max_new)[i][j][k];
			}
		}
	}
	// 	cout << temp << " " << (*option_max_new)[0][0][0] << endl;
	return temp/grid_i/grid_j/grid_k;
	
}

double Optimizer::rms_error(void)
{
	if(not we_have_option_max)return 0.;
	
	double temp=0.;
	
	for(int i=0;i<grid_i;i++){
		for(int j=0;j<grid_j;j++){
			for(int k=0;k<grid_k;k++){
				// 				cout << i << " " << j << " " << k << " " << (*option_max_new)[i][j][k] << endl;
				temp+=
				((*option_max_new)[i][j][k] - (*option_temp)[i][j][k] ) *
				((*option_max_new)[i][j][k] - (*option_temp)[i][j][k] );
			}
		}
	}
	// 	cout << temp << " " << (*option_max_new)[0][0][0] << endl;
	return pow(temp,0.5)/grid_i/grid_j/grid_k;
	
}

double Optimizer::estimate_perpetual_option(double time,double value)
{
	
	double interest_rate_per_hour;
	// interest rate per hour
	interest_rate_per_hour = all_parameters[0]->returnDouble(5)/24./365.;
	
	return value/time/interest_rate_per_hour  ;
	
}

void Optimizer::set_perpetual_option_conditions(double value)
{
	
	if(not we_have_option_max)return;
	
	for(int i=0;i<grid_i;i++){
		for(int j=0;j<grid_j;j++){
			for(int k=0;k<grid_k;k++){
				(*option_max_new)[i][j][k]=value;
				(*option_max_old)[i][j][k]=value;
			}
		}
	}
}

void Optimizer::adjust_perpetual_option_conditions(double value)
{
	
	if(not we_have_option_max)return;
	
	for(int i=0;i<grid_i;i++){
		for(int j=0;j<grid_j;j++){
			for(int k=0;k<grid_k;k++){
				(*option_max_new)[i][j][k]+=value;
				(*option_max_old)[i][j][k]+=value;
			}
		}
	}
}

void Optimizer::run_perpetual(void)
{
	
	for(int ii=0;ii<total_option_value.size();ii++)
	{
		set_base_parameters(ii);	
	}
	
	perpetual_pricer(true);
	
}

int Optimizer::force_commit(double tau)
{
  if(-24.5<tau and tau<-15.5){
	return 1;
  }
  else
	return total_option_value.size();
}

// Try and extrapolate to find the perpetual option value quickly
void Optimizer::perpetual_pricer(bool print_results)
{
	bool debug=false;
	double total_time,average,per_val,average_old,tolerance;
	int contracts_in_a_day;
	// tolerance for convergence
	tolerance = 1.e-6;
	
	// contracts = period/contract_length
	contracts_in_a_day = int( all_parameters[0]->returnDouble(31)
	/all_parameters[0]->returnDouble(22) +0.5);
	// if the option is not periodic method won't work
	if(24%int(all_parameters[0]->returnDouble(22)+0.5)!=0)return;
	// If running for the first time, run a single day option to estimate value
	if(run_starter){
	  // solve once for initial estimate
	  solve_multiple_time(contracts_in_a_day,string("test."),false,false);
	  // time = contracts * contract_length
	  total_time = contracts_in_a_day*all_parameters[0]->returnDouble(22);
	  // calculate average
	  average = option_average();
	  // perpetual = average/interest_rate
	  per_val = estimate_perpetual_option(total_time,average);
	  // debug messages
	  if(debug){
		cout << " Average value : " << average << " after " << total_time << " hours\n";
		cout << "Resulting perpetual option :: " << per_val << endl;
	  }
	  // set up new initial condition
	  //      w = per_val
	  set_perpetual_option_conditions(per_val);
	}
	run_starter = false;
	// run through a few iteration at first to let the solution settle into periodic solution
	average = per_val;
	for(int loop=0;loop<6;loop++){
		// average old allows us to calculate the value added in a day
		average_old = average;
		// solve for the day
		solve_multiple_time(contracts_in_a_day,string("perpetual."),false,true);
		average = option_average();
		// debug messages
		if(debug){
			cout << setprecision(8);
			cout << "it : " << loop << " :: average :: " << 
			average << " :: ratio :: " << 1.-average/average_old << endl;
		}
	}
	
	for(int big_loop=0;big_loop<200;big_loop++){
		// calculate value of average coupon payment each day
		total_time = contracts_in_a_day*all_parameters[0]->returnDouble(22);
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
		average = option_average();
		if(rms_error()<tolerance)break;
		
		// adjust initial conditions so that
		//      w = w + per_val
		adjust_perpetual_option_conditions(per_val);
		// loop again to settle solution
		for(int loop=0;loop<4;loop++){
			average_old = average;
			//
			
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

void Optimizer::grid_test_X_max_per(void)
{
	
	int ii;
	double i_max,x_value,y_value,q_empty,q_full;
	
	cout << " Grid Test..." << endl;
	
	for(int ii=0;ii<total_option_value.size();ii++)
	{
		set_base_parameters(ii);	
	}
	
	stringstream str_stream;
	string temp;
	ofstream x_max;
	
	x_value = all_parameters[0]->returnDouble(7);
	y_value = all_parameters[0]->returnDouble(10);
	q_empty = 0.;
	q_full = all_parameters[0]->returnDouble(3);
	
	
	x_max.open("x_max_data.dat");
	for(int iii=0;iii<total_option_value.size();iii++)
	{
		all_parameters[iii]->changeParameter(101,2);
		all_parameters[iii]->changeParameter(21,2);
		all_parameters[iii]->changeParameter(11,3);
	}
	for(int grid = 0;grid < 101;grid+=25)
	{
		str_stream.str("");
		str_stream << grid;
		temp = "test_grid." + str_stream.str() + ".dat";
		ii = 100+ grid + 1;
		i_max = 100. + grid;
		for(int iii=0;iii<total_option_value.size();iii++)
		{
			all_parameters[iii]->changeParameter(ii,1);// change nodes in X
			all_parameters[iii]->changeParameter(i_max,1); // change X_max
		}
		cout << " Print Parameters... " << endl;
		all_parameters[0]->print_all_params(&cout);
		perpetual_pricer(false);
		print_results(temp,true);
		x_max << grid_i << " " << grid_j <<  " " << grid_k << " " <<
		" " << dx << " " << dy << " " << dz << " Value :: " <<
		option_value_3(x_value,y_value,q_empty,4);
		x_max << " :: " << option_value_3(x_value,y_value,q_full,4) << endl;
	}
	
	x_max.close();
	
}

void Optimizer::increase_timesteps(int i)
{
  for(int ii=0;ii<total_option_value.size();ii++)
  {
	total_option_value[ii]->set_stability_factor(i);
  }
}