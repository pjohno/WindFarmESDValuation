/* Set up a class for the option value */

// #define DEBUG

#include "countdown.h"
#include "option.h"
#include "fixed_grid.h"
#include "storage.h"
#include "param_vecs.h"
#include "finite_diff_methods.h"
#include "nodes.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#define _NEW_MODEL_

using namespace std;

/* 
 * 
 * I nitial Conditions    *
 * 
 */

void Storage:: print_results_old(string filename)
{
  ofstream graph;
  graph.open(filename.c_str());
  // check if the file is open or not
  if (graph.is_open())
  {
    graph << setprecision(8);
    graph.setf(ios::scientific,ios::floatfield);
    graph.fill(' ');
    cout << "Print to file... " << endl;
    // write out stock vs option value
    for(int k = 0;k<kmax;k = k + max(1,kmax/50)){for(int i=0;i<n;i = i + max(1,n/50)){for(int j=0;j<m;j = j + max(1,m/50)){
      graph << X[i]<< " "  << Y[j] << " " << Q[k] << " " << option_new[i][j][k] <<endl;
    }graph << "\n";}graph << "\n";graph << "\n";}
  }
  graph.close();
  cout << "Printing Finished\n";
}

void Storage:: print_results(string filename)
{
  ofstream graph;
  graph.open(filename.c_str());
  double a,b,c;int adjust;
  // check if the file is open or not
  if (graph.is_open())
  {
    graph.setf(ios::scientific,ios::floatfield);
    graph << setprecision(8);
    graph.fill(' ');
    cout << "Print to file... " << endl;
    // write out stock vs option value
    for(int k = 0;k<kmax;k = k + kmax-1){for(int i=0;i<n;i = i + 1){for(int j=0;j<m;j ++){
      GG.grid_Q[i][k].return_stencil(&a,&b,&c,&adjust);
      graph << X[i] <<" "<< Y[j] << " " << Q[k]<< " " << option_new[i][j][k] << " " << b << " " << adjust << endl;
    }graph << "\n";}graph << "\n";graph << "\n";}
  }
  graph.close();
  cout << "Printing Finished\n";
}

void Storage:: print_results(string filename,int j)
{
  ofstream graph;
  graph.open(filename.c_str());
  // check if the file is open or not
  if (graph.is_open())
  {
    graph.setf(ios::scientific,ios::floatfield);
    graph << setprecision(8);
    graph.fill(' ');
    cout << "Print to file... " << endl;
    // write out stock vs option value
    for(int i=0;i<n;i++){
      graph << X[i] << " "  << Y[j] << " " << option_new[i][j][0] << " " << option_new[i][j][kmax-1] <<endl;
    }
  }
  graph.close();
  cout << "Printing Finished\n";
}

void Storage:: print_results(string filename,int j,double param,ofstream *graph)
{
  // check if the file is open or not
  if ((*graph).is_open())
  {
    (*graph).setf(ios::scientific,ios::floatfield);
    (*graph) << setprecision(8);
    (*graph).fill(' ');					
    cout << "Print to file... " << endl;
    // write out stock vs option value
    for(int i=0;i<n;i++){
      (*graph) << X[i] <<" "  << param<< " " << option_new[i][j][0] << " " << option_new[i][j][kmax-1] <<endl;
    }(*graph) << "\n";
  }
  cout << "Printing Finished\n";
}

void Storage:: return_vector_u(int j,int k,vector<double> *u)
{
  (*u).clear();
  for(int i=0; i<option_new.size_x();i++){(*u).push_back(option_new[i][j][k]);}
}

void Storage:: return_vector_uold(int j,int k,vector<double> *u)
{
  (*u).clear();
  for(int i=0; i<option_old.size_x();i++){(*u).push_back(option_old[i][j][k]);}
}


void Storage:: generate_default_params(void)
{
  option_name = "default_option";option_type = 1;
  n = 251;m = 21;kmax = 21;time_steps = 400;
  C_commit = 0.5;X_max = 100.;Y_max = 100.;Q_max = 1.;maturity = 2.;
  interest_rate=0.05;contract_length = 1.;
  mean_reversion_X=0.1;average_X=8.;sigma_X=0.2;
  mean_reversion_Y=0.04;average_Y=40.;sigma_Y=0.075;
  kappa = 0.7;X_C=0.25;X_D=0.25;lambda_C=0.2;lambda_D=0.2;
  tolerance=1.e-3;omega=1.;theta=0.5;rho=0.;	initial_results="default_ic.dat";
  final_results="default_fr.dat";
  plot_file = "default_plot.plt";
  power_cut_out=25.;power_max=1.;power_turn=10.;lambda_cut=2.;
  phase_X = -14.;phase_Y=-2.;
  amp_X=0.375;amp_Y=0.375;
  period_X=24.;period_Y=24.;
  maximum_commit = 1.;minimum_commit = -X_C;
  buy_sell_spread = 0.5;
  contract_rate = 1. - buy_sell_spread;
  penalty_charge = 1. + buy_sell_spread;
  stability_factor = 4;
}

void Storage:: default_parameters(void)
{
  generate_default_params();
}

void Storage:: default_parameters(Param_vecs *my_parameters)
{
  generate_default_params();
  assign_pointers(my_parameters);
}

void Storage:: default_parameters(ParameterMaps *my_parameters)
{
  generate_default_params();
  assign_pointers(my_parameters);
}

void Storage:: assign_pointers(Param_vecs *my_parameters)
{
  
  my_parameters->clear_vecs();
  my_parameters->assignParameter(&option_name ,"Option Name :: ");
  // ints
  my_parameters->assignParameter(&option_type ,"Option Type :: ");
  my_parameters->assignParameter(&n ,"No of X points :: ");
  my_parameters->assignParameter(&m ,"No of Y points :: ");
  my_parameters->assignParameter(&kmax ,"No of Q points :: ");
  my_parameters->assignParameter(&time_steps ,"No of Time Steps :: ");
  // reals
  my_parameters->assignParameter(&C_commit ,"Commit rate :: "); // 0
  my_parameters->assignParameter(&X_max ,"Maximum X value :: ");
  my_parameters->assignParameter(&Y_max ,"Maximum Y value :: ");
  my_parameters->assignParameter(&Q_max ,"Maximum Q value :: ");
  my_parameters->assignParameter(&maturity ,"Option Maturity :: ");
  my_parameters->assignParameter(&interest_rate ,"Interest Rate :: ");
  my_parameters->assignParameter(&mean_reversion_X ,"Mean Reversion in X :: ");
  my_parameters->assignParameter(&average_X ,"Long-term Average in X :: ");
  my_parameters->assignParameter(&sigma_X ,"Variance in X :: ");
  my_parameters->assignParameter(&mean_reversion_Y ,"Mean Reversion in Y :: ");
  my_parameters->assignParameter(&average_Y ,"Long-term Average in Y :: "); // 10
  my_parameters->assignParameter(&sigma_Y ,"Variance in Y :: ");
  my_parameters->assignParameter(&kappa ,"Energy Efficiency :: ");
  my_parameters->assignParameter(&X_C ,"Rate of Charge :: ");
  my_parameters->assignParameter(&X_D ,"Rate of Discharge :: ");
  my_parameters->assignParameter(&lambda_C ,"Smoothing on Charge :: ");
  my_parameters->assignParameter(&lambda_D ,"Smoothing on Discharge :: ");
  my_parameters->assignParameter(&contract_rate ,"Contract Rate :: ");
  my_parameters->assignParameter(&penalty_charge ,"Penalty Rate :: ");
  my_parameters->assignParameter(&power_cut_out ,"Power Cut Out Level :: ");
  my_parameters->assignParameter(&power_max ,"Maximum Power :: "); // 20
  my_parameters->assignParameter(&power_turn ,"Power Turn :: ");
  my_parameters->assignParameter(&contract_length ,"Contract Length :: "); 
  my_parameters->assignParameter(&tolerance ,"Tolerance :: ");
  my_parameters->assignParameter(&omega ,"Omega :: ");
  my_parameters->assignParameter(&theta ,"Differencing parameter (theta) :: ");
  my_parameters->assignParameter(&rho ,"Correlation (rho) :: ");
  my_parameters->assignParameter(&amp_X ," Amplitude of sin wave around av_X :: ");
  my_parameters->assignParameter(&period_X ," Period of sin wave in X (hours) :: ");
  my_parameters->assignParameter(&phase_X ," Phase of sin wave in X (hours) :: ");
  my_parameters->assignParameter(&amp_Y ," Amplitude of sin wave around av_Y :: "); // 30
  my_parameters->assignParameter(&period_Y ," Period of sin wave in Y (hours) :: ");
  my_parameters->assignParameter(&phase_Y ," Phase of sin wave in Y (hours) :: ");
  my_parameters->assignParameter(&maximum_commit ," Maximum Commitment (% power) :: ");
  my_parameters->assignParameter(&minimum_commit ," Minimum Commitment (% power) :: ");
  my_parameters->assignParameter(&buy_sell_spread ," Buy/Sell spread (% clearing price) :: ");
  
  // strings
  my_parameters->assignParameter(&initial_results ,"Initial Results File :: ");
  my_parameters->assignParameter(&final_results ,"Final Results File :: ");
  my_parameters->assignParameter(&plot_file ,"Plot File :: ");
  
}


void Storage:: assign_pointers(ParameterMaps *my_parameters)
{
  my_parameters->clear();
  my_parameters->assignParameter(&option_name,"option_name","Option Name :: ",option_name);
  // ints
  my_parameters->assignParameter(&option_type,"option_type" ,"Option Type :: ",option_type);
  my_parameters->assignParameter(&n ,"n","No of X points :: ",n,2,10000);
  my_parameters->assignParameter(&m ,"m","No of Y points :: ",m,2,10000);
  my_parameters->assignParameter(&kmax ,"kmax","No of Q points :: ",kmax,2,10000);
  my_parameters->assignParameter(&time_steps,"time_steps" ,"No of Time Steps :: ",time_steps,2,10000);
  // reals
  my_parameters->assignParameter(&C_commit,"C_commit" ,"Commit rate :: ",C_commit,-1.,2.); // 0
  my_parameters->assignParameter(&X_max,"X_max" ,"Maximum X value :: ",X_max);
  my_parameters->assignParameter(&Y_max,"Y_max" ,"Maximum Y value :: ",Y_max);
  my_parameters->assignParameter(&Q_max,"Q_max" ,"Maximum Q value :: ",Q_max);
  my_parameters->assignParameter(&maturity,"maturity" ,"Option Maturity :: ",1.);
  my_parameters->assignParameter(&interest_rate,"interest_rate" ,"Interest Rate :: ",0.1,0.,1.);
  my_parameters->assignParameter(&mean_reversion_X,"mean_reversion_X" ,"Mean Reversion in X :: ",0.1);
  my_parameters->assignParameter(&average_X,"average_X" ,"Long-term Average in X :: ",10.);
  my_parameters->assignParameter(&sigma_X,"sigma_X" ,"Variance in X :: ",0.2);
  my_parameters->assignParameter(&mean_reversion_Y,"mean_reversion_Y" ,"Mean Reversion in Y :: ",0.04);
  my_parameters->assignParameter(&average_Y,"average_Y" ,"Long-term Average in Y :: ",40.); // 10
  my_parameters->assignParameter(&sigma_Y,"sigma_Y" ,"Variance in Y :: ",0.075);
  my_parameters->assignParameter(&kappa,"kappa" ,"Energy Efficiency :: ",0.8366666666666,0.,1.);
  my_parameters->assignParameter(&X_C,"X_C" ,"Rate of Charge :: ",0.25,0.,1.);
  my_parameters->assignParameter(&X_D,"X_D" ,"Rate of Discharge :: ",0.25,0.,1.);
  my_parameters->assignParameter(&lambda_C,"lambda_C" ,"Smoothing on Charge :: ",0.2,0.,1.);
  my_parameters->assignParameter(&lambda_D,"lambda_D" ,"Smoothing on Discharge :: ",0.2,0.,1.);
  my_parameters->assignParameter(&contract_rate,"contract_rate" ,"Contract Rate :: ",0.5,0.,1.);
  my_parameters->assignParameter(&penalty_charge,"penalty_charge" ,"Penalty Rate :: ",1.5);
  my_parameters->assignParameter(&power_cut_out,"power_cut_out" ,"Power Cut Out Level :: ",25.);
  my_parameters->assignParameter(&power_max,"power_max" ,"Maximum Power :: ",1.); // 20
  my_parameters->assignParameter(&power_turn,"power_turn" ,"Power Turn :: ",10.);
  my_parameters->assignParameter(&contract_length,"contract_length" ,"Contract Length :: ",1.); 
  my_parameters->assignParameter(&tolerance,"tolerance" ,"Tolerance :: ",0.001);
  my_parameters->assignParameter(&omega,"omega" ,"Omega :: ",1.);
  my_parameters->assignParameter(&theta,"theta" ,"Differencing parameter (theta) :: ",0.5);
  my_parameters->assignParameter(&rho,"rho" ,"Correlation (rho) :: ",0.);
  my_parameters->assignParameter(&amp_X,"amp_X" ," Amplitude of sin wave around av_X :: ",0.375);
  my_parameters->assignParameter(&period_X,"period_X" ," Period of sin wave in X (hours) :: ",24.);
  my_parameters->assignParameter(&phase_X,"phase_X" ," Phase of sin wave in X (hours) :: ",-2.);
  my_parameters->assignParameter(&amp_Y,"amp_Y" ," Amplitude of sin wave around av_Y :: ",0.375); // 30
  my_parameters->assignParameter(&period_Y,"period_Y" ," Period of sin wave in Y (hours) :: ",24.);
  my_parameters->assignParameter(&phase_Y,"phase_Y" ," Phase of sin wave in Y (hours) :: ",-14.);
  my_parameters->assignParameter(&maximum_commit,"maximum_commit" ," Maximum Commitment (% power) :: ",0.);
  my_parameters->assignParameter(&minimum_commit,"minimum_commit" ," Minimum Commitment (% power) :: ",1.);
  my_parameters->assignParameter(&buy_sell_spread,"buy_sell_spread" ," Buy/Sell spread (% clearing price) :: ",0.5);
  
  // strings
  my_parameters->assignParameter(&initial_results,"initial_results" ,"Initial Results File :: ","default_i.dat");
  my_parameters->assignParameter(&final_results,"final_results" ,"Final Results File :: ","default_f.dat");
  my_parameters->assignParameter(&plot_file,"plot_file" ,"Plot File :: ","default_plot.plt");
  
}

void Storage:: assign_grid(void)
{
  //std::cout << "\n assign grid \n";
  
  dQ = Q_max * power_max / double(kmax-1);
  dX = X_max / double(n-1);
  dY = Y_max / double(m-1);
  
  X_D = X_C; // set charge and discharge to the same rate
  if(max(X_C,X_D)==0.){
    dt = dX;
  }else{
    dt = min(dX,dQ/max(X_C,X_D)/power_max);
  }
  time_steps = int(stability_factor*contract_length/dt+0.5);
  time_steps = max(2,time_steps);
  dt = contract_length / double(time_steps);
  
  // set up buy/sell spread
  contract_rate = 1. - buy_sell_spread;penalty_charge = 1. + buy_sell_spread;
  
  lambda_C = max(0.,min(1.,lambda_C)); // lambda_C is the percentage smoothing in Q
  lambda_D = lambda_C;
  
  //	C_commit = max(C_commit,0.01);
  
  double power_curve_zero;
  power_curve_zero = root_finder(C_commit,0.,power_cut_out- 0.5*lambda_cut);
  n_zero_1 = int( power_curve_zero / dX + 0.5);
  // 	cout << "n_zero :: " << n_zero_1 << " Zero :: " << Power_Curve(power_curve_zero) - C_commit << endl;
  
  power_curve_zero = power_cut_out - 0.5*lambda_cut + lambda_cut*(1. - C_commit);
  n_zero_2 = int( power_curve_zero / dX + 0.5);
  // 	cout << "n_zero :: " << n_zero_2*dX << " Zero :: " << Power_Curve(n_zero_2*dX) - C_commit << endl;
  
  if(C_commit<=tolerance)
  {
    n_zero_1 = -1;
    n_zero_2 = n;
  }
  if(C_commit>=1.)
  {
    n_zero_1 = n;
    n_zero_2 = n;
  }
  
  // set up period for mean reversion
  
  freq_X=2.*Pi/period_X;freq_Y=2.*Pi/period_Y;
}

void Storage:: assign_initial_conditions(void)
{
  #ifdef DEBUG
  std::cout << "\n Update grids...\n";
  #endif
  X.setup_x(n,0.,X_max);
  Y.setup_x(m,0.,Y_max);
  Q.setup_x(kmax,0.,Q_max*power_max);
  
  #ifdef DEBUG
  std::cout << "\n Update options...\n";
  #endif
  option_old.update(&X,&Y,&Q);
  option_new.update(&X,&Y,&Q);
  option_old.change_accuracy_z(2);
  option_new.change_accuracy_z(2);
  // assign initial conditions
  option_old.assign_value(&payoff);
  option_new.assign_value(&payoff);
  
  #ifdef DEBUG
  std::cout << "\n Update methods...\n";
  #endif
  
  ADI_time_march.update(X,Y,Q);
  GG.update(X,Y,Q);
  
}


void Storage::reset_initial_conditions(void)
{
  assign_initial_conditions();
}

void Storage::reset_initial_conditions(double ***U)
{
  assign_initial_conditions();
  option_old = U;
  option_new = U;
}

void Storage:: reset_initial_conditions(Option3& U)
{
  assign_initial_conditions();
  option_old = U;
  option_new = U;
}

void Storage:: reset_initial_conditions(double ***u_old,double ***u_new)
{
  assign_initial_conditions();
  option_old = u_old;
  option_new = u_new;
}

void Storage:: reset_initial_conditions(vector<vector<vector<double> > > *u_old,vector<vector<vector<double> > > *u_new)
{
  option_old.clear();option_new.clear();
  assign_initial_conditions();
  for(int i=0;i<n;i++){for(int j=0;j<m;j++){for(int k=0;k<kmax;k++){option_old[i][j][k] = (*u_old)[i][j][k];}}}
  for(int i=0;i<n;i++){for(int j=0;j<m;j++){for(int k=0;k<kmax;k++){option_new[i][j][k] = (*u_new)[i][j][k];}}}
}
void Storage::runSolver(vector<vector<vector<double> > > *u_old,vector<vector<vector<double> > > *u_new,double tau)
{
  
  assign_grid();
  reset_initial_conditions(u_old,u_new);
  set_time(tau);
  main_solver();
  
}
void Storage:: print_initial_conditions(void)
{
  
  ofstream power_funcs;
  power_funcs.open("power_functions.dat");
  // check if the file is open or not
  if (power_funcs.is_open())
  {
    power_funcs.setf(ios::scientific,ios::floatfield);
    power_funcs << setprecision(8);
    power_funcs.fill(' ');
    if(show_progress)cout << "Print initial conditions... " << endl;
    for(int k=0;k<kmax;k = k + max(1,kmax/50)){
      for(int i=0;i<max(0,n_zero_1);i++){
	power_funcs << X[i] << " " << Q[k] << " " << L_D(X[i],Q[k]) << " " << penalties(X[i],average_Y,Q[k])<< " " << payments(average_Y) << " " << Power_Curve(X[i]) << " " << Power_Stack(average_Y) << endl;
      }
      for(int i=max(0,n_zero_1);i<min(n,n_zero_1+1);i++){
	// 				cout << "test";
	power_funcs << X[i]<< " " << Q[k]<< " " << 0.<< " " << penalties(X[i],average_Y,Q[k])<< " " << payments(average_Y)<< " " << Power_Curve(X[i]) << " " << Power_Stack(average_Y) << endl;
      }
      for(int i=max(0,n_zero_1+1);i<min(n-1,n_zero_2);i++){
	power_funcs << X[i]<< " " << Q[k]<< " "  << L_C(X[i],Q[k])<< " " << penalties(X[i],average_Y,Q[k])<< " " << payments(average_Y)<< " " << Power_Curve(X[i]) << " " << Power_Stack(average_Y) << endl;
      }
      for(int i=n_zero_2;i<min(n,n_zero_2+1);i++){
	// 			  cout << "test";
	power_funcs << X[i]<< " " << Q[k]<< " " << 0.<< " " << penalties(X[i],average_Y,Q[k])<< " " << payments(average_Y)<< " " << Power_Curve(X[i]) << " " << Power_Stack(average_Y) << endl;
      }
      for(int i=min(n-1,n_zero_2)+1;i<n;i++){
	power_funcs << X[i]<< " " << Q[k]<< " " << L_D(X[i],Q[k])<< " " << penalties(X[i],average_Y,Q[k])<< " " << payments(average_Y)<< " " << Power_Curve(X[i])<< " " << Power_Stack(average_Y) << endl;
      }power_funcs << "\n";
    }power_funcs << "\n";power_funcs << "\n";
    for(int i=0;i<n;i++){
      power_funcs << X[i]<< " " ;
      power_funcs << Power_Curve(X[i])<< endl  ;
    }power_funcs << "\n";power_funcs << "\n";
    for(int j=0;j<m;j++){
      power_funcs << Y[j]<< " " ;
      power_funcs << Power_Stack(Y[j])<< endl  ;
    }
  }
  power_funcs.close();
  if(show_progress)cout << "Printing Finished\n";
}

double Storage:: root_finder(double cc,double a,double b)
{
  
  double bottom,middle,top;
  bottom = Power_Curve(a) - cc; // bottom < 0
  top = Power_Curve(b) - cc; // top > 0
  middle = Power_Curve((a+b)/2.) - cc;
  if(fabs(middle)<tolerance)return (a+b)/2.;
  for(int loop = 1;loop<1000;loop++)
  {
    
    if( ( middle < 0. and bottom > 0. ) or ( middle > 0. and bottom < 0. ) )
    {
      top = middle;
      b = (a+b)/2.;
      middle = Power_Curve((a+b)/2.) - cc;
      if(fabs(middle)<tolerance)return (a+b)/2.;
    }
    else if ( ( middle < 0. and top > 0. ) or ( middle > 0. and top < 0. ) )
    {
      bottom = middle;
      a = (a+b)/2.;
      middle = Power_Curve((a+b)/2.) - cc;
      if(fabs(middle)<tolerance)return (a+b)/2.;
    }
    
  }
  
}

double Storage:: average_X_t(double tau)
{
  
  return average_X * (1 + amp_X*sin(freq_X*(tau-phase_X)));
  
}

double Storage:: average_X_t_dash(double tau)
{
  return average_X * amp_X*freq_X*cos(freq_X*(tau-phase_X));
}

double Storage:: average_Y_t(double tau)
{
  
  return average_Y * (1 + amp_Y*sin(freq_Y*(tau-phase_Y)));
  
}

double Storage:: average_Y_t_dash(double tau)
{
  
  return average_Y * amp_Y*freq_Y*cos(freq_Y*(tau-phase_Y));
  
}

double Storage:: Power_Stack(double yy)
{
  return yy;
}

double Storage:: Power_Curve(double xx)
{
  if(xx<tolerance)
  {
    return 0.;
  }
  else if(xx<power_cut_out-2.)
  {
    return 0.5*(erf(0.25*( xx - power_turn ) ) + 1. );
  }
  else if(xx<power_cut_out-0.5*lambda_cut)
  {
    return 1.;
  }
  else if(xx<power_cut_out+0.5*lambda_cut)
  {
    return 1. - 1./lambda_cut*(xx-power_cut_out+0.5*lambda_cut);
  }
  else
  {
    return 0.;
  }
}

// energy lost as it passes into the store 
double Storage:: L_C(double xx,double qq)
{
  return kappa * power_max*min(Power_Curve(xx)-C_commit,min(X_C,X_C/lambda_C/Q_max*(Q_max-qq/power_max)));
}

double Storage:: L_C(double xx,double qq,double cc)
{
  return kappa * power_max*min(Power_Curve(xx)-cc,min(X_C,X_C/lambda_C/Q_max*(Q_max-qq/power_max)));
}

double Storage:: L_D(double xx,double qq)
{
  return - power_max*min((-Power_Curve(xx)+C_commit)/kappa,min(X_D,X_D/lambda_D/Q_max/power_max*qq));
}

double Storage:: L_D(double xx,double qq,double cc)
{
  return - power_max*min((-Power_Curve(xx)+cc)/kappa,min(X_D,X_D/lambda_D/Q_max/power_max*qq));
}

// sell at a discount below the market price
// buy at a premium above the market price
double Storage:: penalties(double xx,double yy,double qq)
{
  // when we sell electricity not allowed to sell more than the maximum commit
  return penalty_charge * Power_Stack(yy) * min( power_max*(Power_Curve(xx) - C_commit) - kappa * L_D(xx,qq) , 0.) 
  + contract_rate * Power_Stack(yy) * max( min(power_max*(Power_Curve(xx) - C_commit) - L_C(xx,qq)/kappa,power_max*(maximum_commit-C_commit)) , 0.)
#ifdef _NEW_MODEL_
  + Power_Stack(yy) * power_max * C_commit
#endif
  ;
  
}

double Storage::charge_rate(double xx,double qq,double cc)
{
  if(Power_Curve(xx)<cc)
    return kappa * L_D(xx,qq,cc);
  else
    return L_C(xx,qq,cc)/kappa;
  
}

// power output from the farm
double Storage::power_output(double xx,double qq,double cc)
{
  if(Power_Curve(xx)<cc)
    return power_max*Power_Curve(xx) - kappa * L_D(xx,qq,cc);
  else
    return power_max*Power_Curve(xx) - L_C(xx,qq,cc)/kappa;
}

// payed at the clearing rate for the hour ahead for chosen commitment
double Storage:: payments(double yy)
{
  #ifdef _NEW_MODEL_
return 0.;
#endif
return contract_length * Power_Stack(yy) * power_max * C_commit;
}

/* 
 * 
 * S olver...             *
 * 
 */

void Storage:: setup_stencils(double tau)
{
  
  double aa_,bb_,cc_,g_,adjust_;
  
  double r_over_sigXsquared,root_r_over,scnd_coeff_X,frst_coeff_X,scnd_coeff_Y,frst_coeff_Y;
  double interest_rate_hour,sigma_X_hour,sigma_Y_hour,temp;
  
  
  // adjust parameters units
  interest_rate_hour = interest_rate / 365. / 24.;
  sigma_X_hour = sigma_X ; // quote sigma x in hours
  sigma_Y_hour = sigma_Y ; // quote sigma price in hours
  
  r_over_sigXsquared = 2. * interest_rate_hour / sigma_X_hour / sigma_X_hour;
  root_r_over = pow(r_over_sigXsquared,0.5);
  scnd_coeff_X = 0.5*sigma_X_hour*sigma_X_hour/dX/dX;
  scnd_coeff_Y = 0.5*sigma_Y_hour*sigma_Y_hour/dY/dY;
  
  
  GG.grid_const = -interest_rate_hour;
  
  frst_coeff_X = (mean_reversion_X*average_X_t(tau) + average_X_t_dash(tau) )/ dX;
  aa_ = 1.*frst_coeff_X + interest_rate_hour;
  bb_ = -1.*frst_coeff_X; 
  cc_ = 0.*frst_coeff_X;
  
  GG.grid_X[0](aa_,bb_,cc_,1);
  
  for(int i=1;i<n-1;i++){
    frst_coeff_X = ( mean_reversion_X * (average_X_t(tau)-X[i]) + average_X_t_dash(tau) )/ dX;
    aa_ = -scnd_coeff_X*X[i]*X[i] + 0.5*frst_coeff_X;
    bb_ = 2.*scnd_coeff_X*X[i]*X[i] + interest_rate_hour;
    cc_ = -scnd_coeff_X*X[i]*X[i] - 0.5*frst_coeff_X;
    if(aa_>0.)
    {
      aa_ = -scnd_coeff_X*X[i]*X[i];
      bb_ = 2.*scnd_coeff_X*X[i]*X[i] + frst_coeff_X + interest_rate_hour;
      cc_ = -scnd_coeff_X*X[i]*X[i] - frst_coeff_X;
    }
    else if(cc_>0.)
    {
      aa_ = -scnd_coeff_X*X[i]*X[i] + frst_coeff_X ;
      bb_ = 2.*scnd_coeff_X*X[i]*X[i] - frst_coeff_X + interest_rate_hour;
      cc_ = -scnd_coeff_X*X[i]*X[i];
      // std::cout << "\n X = " << X[i] << " aa_ = " << aa_ << " cc_ = " << cc_ << " at time " << tau << "\n";
    }
    GG.grid_X[i](aa_,bb_,cc_,0);
  }
  
  frst_coeff_X = ( mean_reversion_X * (average_X_t(tau)-X[n-1]) + average_X_t_dash(tau) )/ dX;
  aa_ = -0.*frst_coeff_X;
  bb_ = 1.*frst_coeff_X;
  cc_ = -1.*frst_coeff_X + interest_rate_hour;
  
  GG.grid_X[n-1](aa_,bb_,cc_,-1);
  
  frst_coeff_Y = (mean_reversion_Y*average_Y_t(tau) + average_Y_t_dash(tau)) / dY; 
  aa_ = 1.*frst_coeff_Y;
  bb_ = -1.*frst_coeff_Y;
  cc_ = 0.*frst_coeff_Y;
  
  GG.grid_Y[0](aa_,bb_,cc_,1);
  
  for(int j=1;j<m;j++){
    frst_coeff_Y = (mean_reversion_Y*(average_Y_t(tau) - Y[j]) + average_Y_t_dash(tau)) / dY; 
    aa_ = -scnd_coeff_Y*Y[j]*Y[j] + 0.5*frst_coeff_Y;
    bb_ = 2.*scnd_coeff_Y*Y[j]*Y[j];
    cc_ = -scnd_coeff_Y*Y[j]*Y[j] - 0.5*frst_coeff_Y;
    GG.grid_Y[j](aa_,bb_,cc_,0);
  }
  
  frst_coeff_Y = (mean_reversion_Y*(average_Y_t(tau) - Y[m-1]) + average_Y_t_dash(tau)) / dY; 
  aa_ = -0.*frst_coeff_Y;
  bb_ = 1.*frst_coeff_Y;
  cc_ = -1.*frst_coeff_Y;
  
  GG.grid_Y[m-1](aa_,bb_,cc_,-1);
  
}

void Storage:: setup_scheme(void)
{
  
  //	GG.setup_grid(n,m,kmax);
  
  setup_stencils(theta*dt);
  
  double aa_,bb_,cc_,g_,adjust_,temp;
  
  double ddq;
  ddq=0.5*dQ;
  
  for(int i=0;i<max(0,min(n_zero_1,n));i++){
    for(int j=0;j<m;j++){GG[i][j][0].set_g(penalties(X[i],Y[j],Q[0]));}
    // middle points type = 8 X<0
    
    for(int k=1;k<kmax;k++){			
      cc_ = 0.;
      aa_ = L_D(X[i],Q[k]-ddq)/dQ;
      bb_ = -L_D(X[i],Q[k]-ddq)/dQ;
      if(fabs(L_D(X[i],Q[k]-ddq))>tolerance){
	for(int j=0;j<m;j++){GG[i][j][k].set_g(penalties(X[i],Y[j],Q[k]-ddq));}
	GG.grid_Q[i][k](aa_,bb_,cc_,-1);
      }
      else{
	for(int j=0;j<m;j++){GG[i][j][k].set_g(penalties(X[i],Y[j],Q[k]));}
	GG.grid_Q[i][k](0.,0.,0.,0);
      }
    }
  }
  
  for(int i=max(0,n_zero_1);i<max(0,min(n_zero_1+1,n));i++){
    for(int k=0;k<kmax;k++){for(int j=0;j<m;j++){GG[n_zero_1][j][k].set_g(penalties(X[i],Y[j],Q[k]));}}
  }
  
  for(int i=max(0,n_zero_1+1);i<min(n,n_zero_2);i++){
    // middle points type = 2 X>0
    for(int k=0;k<kmax-1;k++){
      bb_ = L_C(X[i],Q[k]+ddq)/dQ;
      aa_ = -L_C(X[i],Q[k]+ddq)/dQ;
      cc_ = 0.;
      if(fabs(L_C(X[i],Q[k]+ddq))>tolerance){
	for(int j=0;j<m;j++){GG[i][j][k].set_g(penalties(X[i],Y[j],Q[k]+ddq));}
	GG.grid_Q[i][k](aa_,bb_,cc_,1);
      }
      else{
	for(int j=0;j<m;j++){GG[i][j][k].set_g(penalties(X[i],Y[j],Q[k]));}
	GG.grid_Q[i][k](0.,0.,0.,0);
      }
    }
    for(int j=0;j<m;j++){GG[i][j][kmax-1].set_g(penalties(X[i],Y[j],Q[kmax-1]));}
  }
  
  for(int i=min(n,n_zero_2);i<min(n,n_zero_2+1);i++){
    for(int k=0;k<kmax;k++){for(int j=0;j<m;j++){GG[n_zero_2][j][k].set_g(penalties(X[i],Y[j],Q[k]));}}
  }
  
  for(int i=min(n,n_zero_2+1);i<n;i++){
    for(int j=0;j<m;j++){GG[i][j][0].set_g(penalties(X[i],Y[j],Q[0]));}
    // middle points type = 8 X<0
    for(int k=1;k<kmax;k++){		
      cc_ = 0.;
      aa_ = L_D(X[i],Q[k]-ddq)/dQ;
      bb_ = -L_D(X[i],Q[k]-ddq)/dQ;
      if(fabs(L_D(X[i],Q[k]-ddq))>tolerance){
	for(int j=0;j<m;j++){GG[i][j][k].set_g(penalties(X[i],Y[j],Q[k]-ddq));}
	GG.grid_Q[i][k](aa_,bb_,cc_,-1);
      }
      else{
	for(int j=0;j<m;j++){GG[i][j][k].set_g(penalties(X[i],Y[j],Q[k]));}
	GG.grid_Q[i][k](0.,0.,0.,0);
      }
    }
  }
  
}


void Storage:: adjust_scheme_Y(void)
{
  
  int adjust_x,adjust_y,adjust_q;
  double aa_x,bb_x,cc_x;
  double bb_y,aa_y,cc_y;
  double bb_q,aa_q,cc_q,r;
  
  double ZZ_;
  
  
  for(int i=0;i<n;i++){GG.grid_X[i](&aa_x,&bb_x,&cc_x,&adjust_x);
    for(int k=0;k<kmax;k++){GG.grid_Q[i][k](&aa_q,&bb_q,&cc_q,&adjust_q);
      if(adjust_q==0){
	for(int j=0;j<m;j++){
	  ZZ_ = GG[i][j][k]() 
	  - (aa_x*option_new[i+adjust_x-1][j][k+adjust_q] + 
	  bb_x*option_new[i+adjust_x][j][k+adjust_q] + 
	  cc_x*option_new[i+adjust_x+1][j][k+adjust_q] )
	  + 2./dt*option_new[i][j][k+adjust_q];
	  GG[i][j][k].set_Z(ZZ_);
	}
      }else{
	for(int j=0;j<m;j++){GG.grid_Y[j](&aa_y,&bb_y,&cc_y,&adjust_y);					
	  ZZ_ = GG[i][j][k]() 
	  - (aa_x*option_new[i+adjust_x-1][j][k+adjust_q] + 
	  bb_x*option_new[i+adjust_x][j][k+adjust_q] + 
	  cc_x*option_new[i+adjust_x+1][j][k+adjust_q] )
	  - 0.5*(aa_y*option_old[i][j-1+adjust_y][k] + 
	  bb_y*option_old[i][j+adjust_y][k] + 
	  cc_y*option_old[i][j+1+adjust_y][k] )
	  - 2.*aa_q*option_new[i][j][k+adjust_q]
	  + (1./dt - 2.*bb_q)*option_old[i][j][k] ;
	  GG[i][j][k].set_Z(ZZ_);
	}}}}
}



void Storage:: adjust_scheme_X(void)
{
  
  int adjust_x,adjust_y,adjust_q;
  double aa_x,bb_x,cc_x,r;
  double bb_y,aa_y,cc_y;
  double bb_q,aa_q,cc_q;
  
  double ZZ_;
  
  
  // 	cout << "\n Adjust X \n";
  //	r = GG.grid_const;
  GG.grid_X[0](&aa_x,&bb_x,&cc_x,&adjust_x);
  
  for(int j=0;j<m;j++){GG.grid_Y[j](&aa_y,&bb_y,&cc_y,&adjust_y);
    for(int k=0;k<kmax;k++){
      for(int i=0;i<n;i++){GG.grid_Q[i][k](&aa_q,&bb_q,&cc_q,&adjust_q);
	
	// This hints that option_old must have an appropriate value at every point in the grid. Need to make sure we write the correct value to each point.
	if(adjust_q==0)
	{
	  ZZ_ = GG[i][j][k]()
	  - (	aa_y*option_new[i][j+adjust_y-1][k] + 
	  bb_y*option_new[i][j+adjust_y][k] + 
	  cc_y*option_new[i][j+adjust_y+1][k])
	  + 2./dt*option_new[i][j][k];
	  GG[i][j][k].set_Z(ZZ_);
	}else{
	  ZZ_ = GG[i][j][k]()
	  - 0.5*(	aa_y*option_new[i][j+adjust_y-1][k] + 
	  bb_y*option_new[i][j+adjust_y][k] + 
	  cc_y*option_new[i][j+adjust_y+1][k])
	  - 0.5*(	aa_y*option_old[i][j+adjust_y-1][k] + 
	  bb_y*option_old[i][j+adjust_y][k] + 
	  cc_y*option_old[i][j+adjust_y+1][k])
	  - (2.*aa_q + 1./dt)*option_new[i][j][k]
	  + 1./dt*option_old[i][j][k];
	  GG[i][j][k].set_Z(ZZ_);
	}
      }
    }
  }
}

void Storage:: smooth_initial_conditions(void)
{
  double a,b,c,q;
  int adjust;
  for(int k=0;k<kmax;k++){
    for(int i=0;i<n;i++){
      adjust = GG.grid_Q[i][k]();
      q = Q[k]+0.5*Q.dx()*adjust;
      for(int j=0;j<m;j++)option_old[i][j][k] = option_new(i,j,q);
    }
  }
  #ifdef DEBUG
  std::cout << " smooth initial conditions...\n";
  #endif
}

void Storage:: main_solver(void)
{
  
  /* 
   *	 H ere we must solve the problem::
   *	 dV/dt + 1/2 sigma_x^2 x^2 d^2V/dx^2 + 1/2 sigma_y^2 y^2 d^2V/dy^2 
   *	 + alpha_x(theta_x-x)dV/dx + alpha_y(theta_y-y)dV/dy  + L dV/dQ  
   *		= rV + Payments + Penalties
   */
  
  if(show_progress)cout << "Solving...\n";
  
  Countdown screen_output;
  
  // setup grid 
  setup_scheme();
  
  int iter,temp_int,maturity_steps;
  double tau;
  string temp_str;
  
  // set time
  tau = global_time;
  
  smooth_initial_conditions();
  
  // calculate number of contracts in option
  maturity_steps = int(maturity/contract_length +0.5);
  for(int maturity_loop = 0;maturity_loop<maturity_steps;maturity_loop++){
    
    // hourly contract
    for(int time_loop = 1;time_loop<=time_steps;time_loop++){
      // write progress to screen
      if(show_progress)screen_output.update(time_loop + maturity_loop*time_steps,maturity_steps*time_steps);
      //finished writing progress
      //if(C_commit<=tolerance)smooth_initial_conditions();
      
      //solve at time tau
      setup_stencils(tau - theta*dt);
      adjust_scheme_Y();
      ADI_time_march.ADI_imp_Y(dt,option_new(),&GG);
      adjust_scheme_X();
      option_old = option_new;
      ADI_time_march.ADI_imp_X(dt,option_new(),&GG);
      tau = tau-dt;
    }
    
    // receive payment at beginning of the hour at the spot rate
    for(int i=0;i<n;i++){
      for(int j = 0;j<m;j++){
	for(int k = 0;k<kmax;k++){
	  option_new[i][j][k] += payments(Y[j]);option_old[i][j][k] += payments(Y[j]);
	}
	
      }
      
    }
    
  }
  
  if(show_progress){
    cout << "\nSolution found.\n";
  }
  
}

void Storage::generate_pointers(Param_vecs *my_parameters)
{
  assign_pointers(my_parameters);
}

void Storage::generate_pointers(ParameterMaps *my_parameters)
{
  assign_pointers(my_parameters);
}

void Storage::runSolver()
{
  
  assign_grid();
  assign_initial_conditions();
  set_time(0.);
  main_solver();
  
}

void Storage::test_initial_conditions()
{
  
  assign_grid();
  reset_initial_conditions();
  set_time(0.);
  print_initial_conditions();
  
}

void Storage::runSolver(double ***U,double tau)
{
  
  assign_grid();
  reset_initial_conditions(U);
  set_time(tau);
  main_solver();
  
}

void Storage::runSolver(double ***u_old,double ***u_new,double tau)
{
  
  assign_grid();
  reset_initial_conditions(u_old,u_new);
  set_time(tau);
  main_solver();
  
}

void Storage::runSolver(Option3& U,double tau)
{
  assign_grid();
  reset_initial_conditions(U);
  set_time(tau);
  main_solver();
}

void Storage::set_time(double tau)
{
  global_time = tau;
}

void Storage::set_stability_factor(int i)
{
  if(i>0)stability_factor = i;
}

/*
 * 
 * C onstructors...      *
 * 
 */

Storage::Storage()
{
  
  default_parameters();
  show_progress = true;
}


Storage::Storage(Param_vecs *my_parameters,bool progress_bar)
{
  
  default_parameters(my_parameters);
  show_progress = progress_bar;
}


Storage::Storage(ParameterMaps *my_parameters,bool progress_bar)
{
  default_parameters(my_parameters);
  show_progress = progress_bar;
}


