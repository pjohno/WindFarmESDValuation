#include "mtrand.h"
#include "generic_SDE.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

Generic_SDE::~Generic_SDE()
{
	for(int i=0;i<no_of_runs;i++)
	{
	  delete [] normal_dist[i];
	  delete [] X_t[i];
	}
	delete [] normal_dist;
	delete [] X_t;
	delete [] time;
}

Generic_SDE::Generic_SDE()
{
	no_of_runs = 1;
	*normal_dist = NULL;
	*X_t = NULL;
	time=NULL;
	no_of_points = 1000;
	T = 24.;
	seeder = 0;
	beta = 1.;
	alpha_X = 0.1;
	average_X = 8;
	amp_X = 0.375;
	freq_X = 24./2./Pi;
	phase_X = -2.;
	sigma_X = 0.2;
	initial_conditions(10.);
}

Generic_SDE::Generic_SDE(int _points,int _seed,double _T,double _start,double _beta)
{
	no_of_runs = 1;
	*normal_dist = NULL;
	*X_t = NULL;
	time=NULL;
	no_of_points = _points;
	seeder = _seed;
	T = _T;
	beta = _beta;
	alpha_X = 0.1;
	average_X = 8;
	amp_X = 3;
	freq_X = 24./2./Pi;
	phase_X = -2.;
	sigma_X = 0.2;
	initial_conditions(_start);
}

Generic_SDE::Generic_SDE(int _runs,int _points,int _seed,double _T,double _start,double _beta,double _alpha_X,double _average_X,double _amp_X,double _phase_X,double _sigma_X)
{
	no_of_runs = _runs;
	normal_dist = NULL;
	X_t = NULL;
	time=NULL;
	no_of_points = _points;
	seeder = _seed;
	T = _T;
	beta = _beta;
	alpha_X = _alpha_X;
	average_X = _average_X;
	amp_X = _amp_X;
	freq_X = 2.*Pi/24.;
	phase_X = _phase_X;
	sigma_X = _sigma_X;
	initial_conditions(_start);
	generate_dist();
}



/*
 * For demand...
 * alpha = 0.25, theta = 40. + 15.*sin(2.*Pi/24.*(t+6)) and sigma = 0.75 seem to
 * produce a reasonable distribution
 * 
 * For wind...
 * alpha = 0.1, theta = 10. + 10.*sin(2.*Pi/24.*(t+6));and sigma = 0.2
*/
double Generic_SDE::alpha(double X,double t)
{
	return alpha_X;
}

double Generic_SDE::theta(double X,double t)
{
	return average_X * (1 + amp_X*sin(freq_X*(t-phase_X)));
}

double Generic_SDE::sigma(double X,double t)
{
	return sigma_X;
}

double Generic_SDE::forcing(double X,double t)
{
	return average_X*freq_X*amp_X*cos(freq_X*(t-phase_X));
}

void Generic_SDE::initial_conditions(double _start)
{
  bool debug = false;
  if(debug)std::cout << " initial conditions \n";
	X_start = _start;
	dt = T/double(no_of_points);
	normal_dist = new double*[no_of_runs];
	for(int i=0;i<no_of_runs;i++)
	{
	  normal_dist[i] = NULL;
	  gas_dev(seeder+i,i);
	}
}

void Generic_SDE::generate_dist(void)
{
  
  bool debug = false;
  if(debug)std::cout << " generate dist \n";
	if(X_t){
	  for(int runs=0;runs<no_of_runs;runs++)
	  {
		delete [] X_t[runs];
	  }
	}
	delete [] X_t;
	X_t = new double*[no_of_runs];
	for(int runs=0;runs<no_of_runs;runs++)
	{
	  X_t[runs] = new double[no_of_points+1];
	  X_t[runs][0] = X_start;
	}
	delete [] time;
	time = new double[no_of_points+1];
	time[0] = 0.;
	for(int i=0;i<no_of_points;i++)
	{
		time[i+1] = time[i] + dt;
		for(int runs=0;runs<no_of_runs;runs++)
		{
		  X_t[runs][i+1] = X_t[runs][i] + alpha(X_t[runs][i],time[i]) * (theta(X_t[runs][i],time[i]) - X_t[runs][i]) * dt 
		  + forcing(X_t[runs][i],time[i])*dt
		  + sigma(X_t[runs][i],time[i])*pow(X_t[runs][i],beta)*pow(dt,0.5)*normal_dist[runs][i];
		}
	}
}

double Generic_SDE::mean_dist(double *distribution,int n)
{
	double temp=0.;
	for(int i=0;i<no_of_points;i++)
	{
		temp += distribution[i];
	}
	return temp/double(n+1);
}

double Generic_SDE::variance_dist(double *distribution,int n)
{
	double temp=0.,temp_mean;
	temp_mean = mean_dist(distribution,n);
	for(int i=0;i<no_of_points;i++)
	{
		temp += (distribution[i] - temp_mean)*(distribution[i] - temp_mean);
	}
	return temp/double(n+1);
}

void Generic_SDE::print_dist(void)
{
	
	std::ofstream graph;
	graph.open("results.data");
	// check if the file is open or not
	if (graph.is_open())
	{
		graph << std::setprecision(8);
		graph.setf(std::ios::scientific,std::ios::floatfield);
		graph.fill(' ');
		std::cout << "Print to file... \n";
		// write out time vs SDE
		for(int i = 0;i<=no_of_points;i = i++){
			graph << time[i] << " " << X_t[0][i] << " " << theta(X_t[0][i],time[i]) <<" " << Power_Stack(X_t[0][i]) << std::endl;
		}
	}
	graph.close();
	std::cout << "Printing Finished\n";
}

double Generic_SDE:: Power_Stack(double yy)
{
	return yy;
}

void Generic_SDE::solve(void)
{
	generate_dist();
	print_dist();
	std::cout << " Mean :: " << mean_dist(normal_dist[0],no_of_points) 
		<< " Variance :: " << variance_dist(normal_dist[0],no_of_points) << std::endl;
	
}

void Generic_SDE::gas_dev(int seed,int runs)
{
	  bool debug = false;
  if(debug)std::cout << " gas dev ... \n";

	delete [] normal_dist[runs];
	MTRand drand(seed);
	int i=0;
	double v1,v2,rsq;
	int n=no_of_points;
	if(n==0)return;
	if(n%2 == 1)
	{
		n=+1;
	}
	normal_dist[runs] = new double[n];
	do
	{
		
		v1 = drand();
		v2 = drand();
		v1 = 2.*v1-1.;
		v2 = 2.*v2-1.;
		rsq = v1*v1+v2*v2;
		
		if(rsq<1.)
		{
			rsq=pow(-2.*log(rsq)/rsq , 0.5);
			normal_dist[runs][i] = v1*rsq;
			normal_dist[runs][i+1] = v2*rsq;
			i = i + 2;
		}
		
	}while(i<n);
}

