#include "simulator.h"
#include "newton_iteration.h"
#include "generic_SDE.h"
#include "optimizer.h"
#include "param_vecs.h"
#include "storage.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

Simulator::~Simulator()
{
	delete monte_X;
	delete monte_Y;
}

Simulator::Simulator()
{
   commit_strategy_setup=false;generate_distributions=false;
   monte_X = NULL;monte_Y=NULL;
}

Simulator::Simulator(int n,int m,int p,int q,int no_of_options,bool progress_bar,int no_price):Optimizer(no_of_options,progress_bar,no_price)
{
  commit_strategy_setup=true;generate_distributions=false;
  double contract_length;
  contract_length = int( 24./double(q) + 0.5 );
  for(int ii=0;ii<total_option_value.size();ii++)
  {
	set_base_parameters(ii);
	all_parameters[ii]->changeParameter(n,1);
 	all_parameters[ii]->changeParameter(m,2);
 	all_parameters[ii]->changeParameter(p,3);
	all_parameters[ii]->changeParameter(contract_length,22);
  }
  
  std::cout << " Make commit vector with n : " << n << " m : " << m << " p : " << p << " q : " << q << "\n";
  
  commit_strategy.resize(n);
  for(int i=0;i<n;i++)
  {
	commit_strategy[i].resize(m);
	for(int j=0;j<m;j++)
	{
	  commit_strategy[i][j].resize(p);
	  for(int k=0;k<p;k++)
	  {
		commit_strategy[i][j][k].resize(q);
	  }
	}
  }
  
  double alpha_X,average_X,amp_X,freq_X,phase_X,sigma_X;
  double alpha_Y,average_Y,amp_Y,freq_Y,phase_Y,sigma_Y;
  
  alpha_X = all_parameters[0]->returnDouble(6);
  average_X = all_parameters[0]->returnDouble(7);
  sigma_X = all_parameters[0]->returnDouble(8);
  alpha_Y = all_parameters[0]->returnDouble(9);
  average_Y = all_parameters[0]->returnDouble(10);
  sigma_Y = all_parameters[0]->returnDouble(11);
  amp_X = all_parameters[0]->returnDouble(27);
  phase_X = all_parameters[0]->returnDouble(29);
  amp_Y= all_parameters[0]->returnDouble(30);
  phase_Y= all_parameters[0]->returnDouble(32);
  
  X = total_option_value[0]->average_X_t(0);
  Y = total_option_value[0]->average_Y_t(0);
  no_of_points = 2400;
  no_of_runs = 10;
  dt = 24./no_of_points;

  monte_X = new Generic_SDE(no_of_runs,no_of_points,10,24.,X,1.,alpha_X,average_X,amp_X,phase_X,sigma_X);
  monte_Y = new Generic_SDE(no_of_runs,no_of_points,20,24.,Y,1.,alpha_Y,average_Y,amp_Y,phase_Y,sigma_Y);;

}

void Simulator::solve_multiple_time(int no_of_hours,string filename,bool display,bool reset)
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
		  write_strategy_to_vector(no_of_hours - i - 1);
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


void Simulator::run_perpetual(void)
{
	
  bool debug = false;
  
  if(debug){std::cout << " Using simulator pricer...\n";}
  perpetual_pricer(true);
	
}

void Simulator::write_strategy_to_vector(int hour)
{
  bool debug = false;
  if(debug){std::cout << "\n write commit to vector at hour "<< hour << "\n";
  std::cout << "\n grid size n : " << grid_i << " m : " << grid_j << " p : " << grid_k << "\n";}
  for(int i=0;i<grid_i;i++)
  {
	for(int j=0;j<grid_j;j++)
	{
	  for(int k=0;k<grid_k;k++)
	  {
		commit_strategy[i][j][k][hour] = (*option_commit)[i][j][k];
	  }
	}
  }
}

void Simulator::generate_simulated_day(double store_level,std::string filename)
{
  
  int i,j,k,hour,hour_temp;
  double dQ;
  Q = store_level;
  
  std::ofstream day_sim;
  day_sim.open(filename.c_str());
  for(int runs=0;runs<no_of_runs;runs++){
	Q = store_level;
	hour_temp = -1;
	for(int loop=0;loop<no_of_points;loop++)
	{
	  X = monte_X->X_t[runs][loop];
	  Y = monte_Y->X_t[runs][loop];
	  hour = max(0,int(dt*loop));
	  if(hour_temp<hour){
		i = min(int(X/dx + 0.5),grid_i);
		j = min(int(Y/dy + 0.5),grid_j);
		k = min(int(Q/dz + 0.5),grid_k);
		hour_temp = hour;
	  }
	  C = commit_strategy[i][j][k][hour];
	  dQ = total_option_value[0]->charge_rate(X,Q,C)*dt;
	  day_sim << loop*dt << " " << X << " " << Y << " " << Q << " " << C << " " << output() <<" " << average_x(loop*dt) <<" " << average_y(loop*dt) << "\n";
	  Q = max(0.,min((grid_k-1)*dz,Q + dQ));
	}
	day_sim << "\n\n";
  }
  
  day_sim.close();
}

double Simulator::average_x(double tau)
{
  return total_option_value[0]->average_X_t(tau);
}

double Simulator::average_y(double tau)
{
  return total_option_value[0]->average_Y_t(tau);
}

double Simulator::output(void)
{
  return total_option_value[0]->power_output(X,Q,C);
}

void Simulator::smoothing_store(void)
{
  double spread=0.,store_size=1.,max_charge=0.25;
  
  for(int ii=0;ii<total_option_value.size();ii++)
  {
	all_parameters[ii]->changeParameter(store_size,3);
 	all_parameters[ii]->changeParameter(max_charge,13);
 	all_parameters[ii]->changeParameter(max_charge,14);
	all_parameters[ii]->changeParameter(spread,35);
  }
}

void Simulator::arbitrage_store(void)
{
  double spread=0.,store_size=4.,max_charge=1.;
  
  for(int ii=0;ii<total_option_value.size();ii++)
  {
	all_parameters[ii]->changeParameter(store_size,3);
 	all_parameters[ii]->changeParameter(max_charge,13);
 	all_parameters[ii]->changeParameter(max_charge,14);
	all_parameters[ii]->changeParameter(spread,35);
  }
}

void Simulator::no_store(void)
{
  double spread=0.,store_size=1.,max_charge=0.;
  
  for(int ii=0;ii<total_option_value.size();ii++)
  {
	all_parameters[ii]->changeParameter(store_size,3);
 	all_parameters[ii]->changeParameter(max_charge,13);
 	all_parameters[ii]->changeParameter(max_charge,14);
	all_parameters[ii]->changeParameter(spread,35);
  }
}

void Simulator::test_case(void)
{
  
  for(int ii=0;ii<total_option_value.size();ii++)
  {
	all_parameters[ii]->changeParameter(1.,3);
	// no store
 	all_parameters[ii]->changeParameter(0.,13);
 	all_parameters[ii]->changeParameter(0.,14);
	// zero spread
	all_parameters[ii]->changeParameter(0.,35);
	// set price  variation
	all_parameters[ii]->changeParameter(0.4,9); // Mean Reversion in Y
	all_parameters[ii]->changeParameter(40.,10);// Long-term Average in Y
	all_parameters[ii]->changeParameter(0.075,11);// Variance in Y
	all_parameters[ii]->changeParameter(0.375,30);// Amplitude of sin wave around av_Y
	all_parameters[ii]->changeParameter(-14.,32);// Phase of sin wave in Y (hours)
	// set wind  variation
	all_parameters[ii]->changeParameter(0.1,6);
	all_parameters[ii]->changeParameter(7.7788148,7); // long-term average
	all_parameters[ii]->changeParameter(0.2,8); // variance in X
	all_parameters[ii]->changeParameter(0.375,27);
	// perfect store
	all_parameters[ii]->changeParameter(1.,12);
	// return rate 
	all_parameters[ii]->changeParameter(.1,5);
  }
}

void Simulator::changeParameter(double x,int parameter_ref)
{
  for(int ii=0;ii<total_option_value.size();ii++)
  {
	all_parameters[ii]->changeParameter(x,parameter_ref);
  }
}

double Simulator::hours_per_year(void)
{
  
  cout << "\n Option average is :: " << option_average() << endl; 
  cout << "\n Then yearly output is :: " << option_average()*all_parameters[0]->returnDouble(5)
 << " " << option_value_3(20.,0.,0.,4)<< endl; 
  cout << "\n And percentage output is :: " << option_average()*all_parameters[0]->returnDouble(5)/24./365.
 << endl; 
  return option_average()*all_parameters[0]->returnDouble(5)/24./365.;
}


void Simulator::vary_perpetual(int parameter_ref_1,double low_1,double high_1,int its_1,int parameter_ref_2,double low_2,double high_2,int its_2,string filename)
{
  
  bool file_setup=true,debug=false;
  int temp_int;
  double param_1,param_2,x_value,x_value_star,y_value,y_value_star,q_empty,q_full,dX,d_its_1;
  vector<double> grid_size(3);
  ofstream vary_param,log_file;
  string results_file;
  // Set base parameters
  //for(int ii=0;ii<total_option_value.size();ii++)
  //{
//	set_base_parameters(ii);	
  //}
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
		vary_param << all_parameters[0]->returnDouble(5)
		*option_value_3(x_value_star,y_value_star,q_empty,4) <<  " \n" ;
	  }
	  vary_param << " \n" ;
	  vary_param.flush();
	}
	vary_param.close();
	
  }
  
