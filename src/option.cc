// option class, containing interpolation and other benifits
#include "param_vecs.h"
#include "param_maps.h"
#include "option.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

// #define DEBUG

double Interpolation::interpolate(double x,std::vector<double> y,std::vector<double> g)
{
  double temp=0.;
  //check vectors
  if(y.size()!=g.size()){return 0.;}
  //
  // interpolate = sum{i=1,n} g(i) Pi{j=1,n,j/=i} (x - x_j) / (x_i - x_j)
  //
  for(int i=0;i<y.size();i++){
	double int_temp;
	int_temp = g[i];
	for(int j=0;j<y.size();j++){
	  if(j==i){continue;}
	  int_temp *= ( x - y[j] )/( y[i] - y[j] );
	}
	temp += int_temp;
  }
  return temp;
};

void Option::update(int n)
{
  if(n==no_of_elements and u_setup)return;
  if(u_setup)delete [] U;
  U = new double[n];
  u_setup = true;
  no_of_elements=n;
  accuracy_x = std::min(4,no_of_elements);
};

double Option::operator()(double x){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_x)/2. - double(accuracy_x/2);
  option_int = int( (x-x_min) / dx  + half_temp) - (accuracy_x-1)/2;
  option_int = std::max(0,std::min(no_of_elements-accuracy_x,option_int));
  for(int i=option_int;i<option_int+accuracy_x;i++){
	yy.push_back(x_min+i*dx);
	gg.push_back(U[i]);
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(x,yy,gg);
}

Option& Option::max(Option& a,Option& b)
{
  for(int i=0;i<no_of_elements;i++)U[i]=std::max(a[i],b[i]);
  return *this;
}

void Print::print(Option &u)
{
  double x,y;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  #endif
  for(int i=0;i<n;i++)
  {
	x = x_min + i*dx;
	file <<  x << " " << u(x) << "\n";
  }
  file.close();
}

void Print::link_printer_config(Param_vecs *print_config)
{
  print_config->assignParameter(&precision," Precision ");
  print_config->assignParameter(&n," No of points in X ");
  print_config->assignParameter(&x_min," Minimum value of X ");
  print_config->assignParameter(&x_max," Maximum value of X ");
}

void Print::link_printer_config(ParameterMaps *print_config)
{
  print_config->assignParameter(&precision,"precision"," Precision ");
  print_config->assignParameter(&n,"n"," No of points in X ");
  print_config->assignParameter(&x_min,"x_min"," Minimum value of X ");
  print_config->assignParameter(&x_max,"x_max"," Maximum value of X ");
}

Option2::Option2(Fixed_Grid *x,Fixed_Grid *y)
{
  default_values();
  update(x,y);
  assign_value();
};

Option2::Option2(Fixed_Grid *x,Fixed_Grid *y,Payoff2 *vanilla)
{
  default_values();
  update(x,y);
  assign_value(vanilla);
};

void Option2::default_values(void)
{
  accuracy_x = 4;
  accuracy_y = 4;
  U = NULL;
  u_setup = false;
  n = 4;
  x_min = 0.;
  x_max = 1.;
  m = 4;
  y_min = 0.;
  y_max = 1.;
}

void Option2::assign_value(Payoff2 *vanilla)
{
  for(int i=0;i<n;i++){
	for(int j=0;j<m;j++){
	  U[i][j]=(*vanilla)(x_min+i*dx,y_min+j*dy);
	}
  }
}

void Option2::update(Fixed_Grid *x,Fixed_Grid *y)
{
  setup_x(x->size(),x->min(),x->max());
  setup_y(y->size(),y->min(),y->max());
  if(u_setup)
  {
	if(x->size()==n and y->size()==m)return;
	else clear();
  }
  n=x->size();m=y->size();
  U = new double*[n];
  for(int i=0;i<n;i++){U[i]= new double[m];}
  u_setup = true;
  accuracy_x = std::min(n,accuracy_x);
  accuracy_y = std::min(m,accuracy_y);
};

double Option2::operator()(double x,double y){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_x)/2. - double(accuracy_x/2);
  option_int = int( (x-x_min) / dx  + half_temp) - (accuracy_x-1)/2;
  option_int = std::max(0,std::min(n-accuracy_x,option_int));
  for(int i=option_int;i<option_int+accuracy_x;i++){
	yy.push_back(x_min+i*dx);
	gg.push_back((*this)(i,y));
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(x,yy,gg);
}

double Option2::operator()(int i,double y){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_y)/2. - double(accuracy_y/2);
  option_int = int( (y-y_min) / dy  + half_temp) - (accuracy_y-1)/2;
  option_int = std::max(0,std::min(m-accuracy_y,option_int));
  for(int j=option_int;j<option_int+accuracy_y;j++){
	yy.push_back(y_min+j*dy);
	gg.push_back(U[i][j]);
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(y,yy,gg);
}

Option2& Option2::max(Option2& a,Option2& b)
{
  for(int i=0;i<n;i++)for(int j=0;j<m;j++)U[i][j]=std::max(a[i][j],b[i][j]);
  return *this;
}

void Print2::calc_grid(void)
{
  dx = (x_max - x_min)/double(n-1);
  if(n==1)dx=0.;
  dy = (y_max - y_min)/double(m-1);
  if(m==1)dy=0.;
}

void Print2::print(Option2 &u)
{
  double x,y;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  std::cout << " dx : " << dx << " dy : " << dy << "\n";
  #endif
  for(int i=0;i<n;i++)
  {
	x = x_min + i*dx;
	for(int j=0;j<m;j++)
	{
	  y = y_min + j*dy;
	  file <<  x << " " << y << " " << u(x,y) << "\n";
	}
	file << "\n";
  }
  file.close();
}

void Print2::print_y_x(Option2 &u)
{
  double x,y;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  std::cout << " dx : " << dx << " dy : " << dy << "\n";
  #endif
  for(int j=0;j<m;j++)
  {
    y = y_min + j*dy;
    for(int i=0;i<n;i++)
    {
      x = x_min + i*dx;
      file <<  x << " " << y << " " << u(x,y) << "\n";
    }
    file << "\n";
  }
  file.close();
}

void Print2::print_U_and_V(Option2 &u,Option2 &v)
{
  double x,y;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  std::cout << " dx : " << dx << " dy : " << dy << "\n";
  #endif
  for(int i=0;i<n;i++)
  {
    x = x_min + i*dx;
    for(int j=0;j<m;j++)
    {
      y = y_min + j*dy;
      file <<  x << " " << y << " " << u(x,y) << " " << v(x,y) << "\n";
    }
    file << "\n";
  }
  file.close();
}

void Print2::print_y_x_U_and_V(Option2 &u,Option2 &v)
{
  double x,y;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  std::cout << " dx : " << dx << " dy : " << dy << "\n";
  #endif
  for(int j=0;j<m;j++)
  {
    y = y_min + j*dy;
    for(int i=0;i<n;i++)
    {
      x = x_min + i*dx;
      file <<  x << " " << y << " " << u(x,y) << " " << v(x,y) <<"\n";
    }
    file << "\n";
  }
  file.close();
}

void Print2::link_printer_config2(Param_vecs *print_config)
{
  print_config->assignParameter(&m," No of points in Y ");
  print_config->assignParameter(&y_min," Minimum value of Y ");
  print_config->assignParameter(&y_max," Maximum value of Y ");
}

void Print2::link_printer_config2(ParameterMaps *print_config)
{
  print_config->assignParameter(&m,"m"," No of points in Y ");
  print_config->assignParameter(&y_min,"y_min"," Minimum value of Y ");
  print_config->assignParameter(&y_max,"y_max"," Maximum value of Y ");
}

// 
//
//  3D option value with interpolation
//
//

Option3::Option3(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z)
{
  default_values();
  update(x,y,z);
  assign_value();
};

Option3::Option3(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,Payoff3 *vanilla)
{
  default_values();
  update(x,y,z);
  assign_value(vanilla);
};

void Option3::default_values(void)
{
  accuracy_x = 4;
  accuracy_y = 4;
  accuracy_z = 4;
  U = NULL;
  u_setup = false;
  n = 4;
  x_min = 0.;
  x_max = 1.;
  m = 4;
  y_min = 0.;
  y_max = 1.;
  p = 4;
  z_min = 0.;
  z_max = 1.;
}

void Option3::assign_value(Payoff3 *vanilla)
{
  for(int i=0;i<n;i++){
	for(int j=0;j<m;j++){
	  for(int k=0;k<p;k++){
		U[i][j][k]=(*vanilla)(x_min+i*dx,y_min+j*dy,z_min+k*dz);
	  }
	}
  }
}

void Option3::setup_x(int n_,double x_min_,double x_max_){
  x_min = x_min_;
  x_max = x_max_;
  dx = (x_max - x_min)/double(n_-1);
}

void Option3::setup_y(int n_,double x_min_,double x_max_){
  y_min = x_min_;
  y_max = x_max_;
  dy = (y_max - y_min)/double(n_-1);
}

void Option3::setup_z(int n_,double x_min_,double x_max_){
  z_min = x_min_;
  z_max = x_max_;
  dz = (z_max - z_min)/double(n_-1);
};

void Option3::update(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z)
{
  setup_x(x->size(),x->min(),x->max());
  setup_y(y->size(),y->min(),y->max());
  setup_z(z->size(),z->min(),z->max());
  if(u_setup){
	if(x->size()==n and y->size()==m and z->size()==p)return;
	else clear();
  }
  n=x->size();m=y->size();p=z->size();
  U = new double**[n];
  for(int i=0;i<n;i++)
  {
	U[i]= new double*[m];
	for(int j=0;j<m;j++){U[i][j]= new double[p];}
  }
  u_setup = true;
  accuracy_x = std::min(n,accuracy_x);
  accuracy_y = std::min(m,accuracy_y);
  accuracy_z = std::min(p,accuracy_z);
};

double Option3::operator()(double x,double y,double z){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_x)/2. - double(accuracy_x/2);
  option_int = int( (x-x_min) / dx  + half_temp) - (accuracy_x-1)/2;
  option_int = std::max(0,std::min(n-accuracy_x,option_int));
  for(int i=option_int;i<option_int+accuracy_x;i++){
	yy.push_back(x_min+i*dx);
	gg.push_back((*this)(i,y,z));
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(x,yy,gg);
}

double Option3::operator()(int i,double y,double z){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_y)/2. - double(accuracy_y/2);
  option_int = int( (y-y_min) / dy  + half_temp) - (accuracy_y-1)/2;
  option_int = std::max(0,std::min(m-accuracy_y,option_int));
  for(int j=option_int;j<option_int+accuracy_y;j++){
	yy.push_back(y_min+j*dy);
	gg.push_back((*this)(i,j,z));
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(y,yy,gg);
}

double Option3::operator()(int i,int j,double z){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_z)/2. - double(accuracy_z/2);
  option_int = int( (z-z_min) / dz  + half_temp) - (accuracy_z-1)/2;
  option_int = std::max(0,std::min(p-accuracy_z,option_int));
  for(int k=option_int;k<option_int+accuracy_z;k++){
	yy.push_back(z_min+k*dz);
	gg.push_back(U[i][j][k]);
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(z,yy,gg);
}

double rms(Option3& a,Option3& b)
{
  double sum;
  sum = 0.;
  for(int i=0;i<a.size_x();i++)
	for(int j=0;j<a.size_y();j++)
	  for(int k=0;k<a.size_z();k++)
	  {
		sum+=(a[i][j][k] - b[i][j][k])*(a[i][j][k] - b[i][j][k]);
	  }
  return std::sqrt(sum)/a.size_x()/a.size_y()/a.size_z();
}

Option3& Option3::max(Option3& a,Option3& b)
{
  for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)U[i][j][k]=std::max(a[i][j][k],b[i][j][k]);
  return *this;
}


void Print3::calc_grid(void)
{
  if(n==1)dx=0.;
  else dx = (x_max - x_min)/double(n-1);
  if(m==1)dy=0.;
  else dy = (y_max - y_min)/double(m-1);
  if(p==1)dz=0.;
  else dz = (z_max - z_min)/double(p-1);
}

void Print3::print(Option3 &u)
{
  double x,y,z;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  #endif
  for(int i=0;i<n;i++)
  {
	x = x_min + i*dx;
	for(int j=0;j<m;j++)
	{
	  y = y_min + j*dy;
	  for(int k=0;k<p;k++)
	  {
		z = z_min + k*dz;
		file <<  x << " " << y << " " << z << " " << u(x,y,z) << "\n";
	  }
	  file << "\n";
	}
	file << "\n";	
  }
  file.close();
}

void Print3::print_x_y(Option3 &u)
{
  double x,y,z;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  #endif
  for(int k=0;k<p;k++)
  {
	z = z_min + k*dz;
	for(int j=0;j<m;j++)
	{
	  y = y_min + j*dy;
	  for(int i=0;i<n;i++)
	  {
		x = x_min + i*dx;
		file <<  x << " " << y << " " << z << " " << u(x,y,z) << "\n";
	  }
	  file << "\n";
	}
	file << "\n";	
  }
  file.close();
}

void Print3::print_x_z(Option3 &u)
{
  double x,y,z;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  #endif
  for(int j=0;j<m;j++)
  {
	y = y_min + j*dy;
	for(int k=0;k<p;k++)
	{
	  z = z_min + k*dz;
	  for(int i=0;i<n;i++)
	  {
		x = x_min + i*dx;
		file <<  x << " " << y << " " << z << " " << u(x,y,z) << "\n";
	  }
	  file << "\n";
	}
	file << "\n";	
  }
  file.close();
}

void Print3::link_printer_config3(Param_vecs *print_config)
{
  print_config->assignParameter(&p," No of points in Z ");
  print_config->assignParameter(&z_min," Minimum value of Z ");
  print_config->assignParameter(&z_max," Maximum value of Z ");
}

void Print3::link_printer_config3(ParameterMaps *print_config)
{
  print_config->assignParameter(&p,"p"," No of points in Z ");
  print_config->assignParameter(&z_min,"z_min"," Minimum value of Z ");
  print_config->assignParameter(&z_max,"z_max"," Maximum value of Z ");
}


// 
//
//  4D option value with interpolation
//
//

Option4::Option4(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,Fixed_Grid *w)
{
  default_values();
  update(x,y,z,w);
  assign_value();
};

Option4::Option4(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,Fixed_Grid *w,Payoff4 *vanilla)
{
  default_values();
  update(x,y,z,w);
  assign_value(vanilla);
};

void Option4::default_values(void)
{
  accuracy_x = 4;
  accuracy_y = 4;
  accuracy_z = 4;
  accuracy_w = 4;
  U = NULL;
  u_setup = false;
  n = 4;
  x_min = 0.;
  x_max = 1.;
  m = 4;
  y_min = 0.;
  y_max = 1.;
  p = 4;
  z_min = 0.;
  z_max = 1.;
  q = 4;
  w_min = 0.;
  w_max = 1.;
}

void Option4::assign_value(Payoff4 *vanilla)
{
  for(int i=0;i<n;i++){
	for(int j=0;j<m;j++){
	  for(int k=0;k<p;k++){
		for(int l=0;l<q;l++){
		  U[i][j][k][l]=(*vanilla)(x_min+i*dx,y_min+j*dy,z_min+k*dz,w_min+l*dw);
		}
	  }
	}
  }
}

void Option4::setup_x(int size,double min,double max)
{
  x_min = min;
  x_max = max;
  dx = (max - min)/double(size-1);
}

void Option4::setup_y(int size,double min,double max)
{
  y_min = min;
  y_max = max;
  dy = (max - min)/double(size-1);
}

void Option4::setup_z(int size,double min,double max)
{
  z_min = min;
  z_max = max;
  dz = (max - min)/double(size-1);
}

void Option4::setup_w(int size,double min,double max)
{
  w_min = min;
  w_max = max;
  dw = (max - min)/double(size-1);
}
	
void Option4::update(Fixed_Grid *x,Fixed_Grid *y,Fixed_Grid *z,Fixed_Grid *w)
{
  setup_x(x->size(),x->min(),x->max());
  setup_y(y->size(),y->min(),y->max());
  setup_z(z->size(),z->min(),z->max());
  setup_w(w->size(),w->min(),w->max());
  if(u_setup){
	if(n==x->size() and m==y->size() and p==z->size() and q==w->size())return;
	else clear();
  }
  n=x->size();m=y->size();p=z->size();q=w->size();
  U = new double***[n];
  for(int i=0;i<n;i++)
  {
	U[i]= new double**[m];
	for(int j=0;j<m;j++){
	  U[i][j]= new double*[p];
	  for(int k=0;k<p;k++){
		U[i][j][k]= new double[q];
	  }
	}
  }
  u_setup = true;
  accuracy_x = std::min(n,accuracy_x);
  accuracy_y = std::min(m,accuracy_y);
  accuracy_z = std::min(p,accuracy_z);
  accuracy_w = std::min(q,accuracy_w);
};

double Option4::operator()(double x,double y,double z,double w){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_x)/2. - double(accuracy_x/2);
  option_int = int( (x-x_min) / dx  + half_temp) - (accuracy_x-1)/2;
  option_int = std::max(0,std::min(n-accuracy_x,option_int));
  for(int i=option_int;i<option_int+accuracy_x;i++){
	yy.push_back(x_min+i*dx);
	gg.push_back((*this)(i,y,z,w));
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(x,yy,gg);
}

double Option4::operator()(int i,double y,double z,double w){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_y)/2. - double(accuracy_y/2);
  option_int = int( (y-y_min) / dy  + half_temp) - (accuracy_y-1)/2;
  option_int = std::max(0,std::min(m-accuracy_y,option_int));
  for(int j=option_int;j<option_int+accuracy_y;j++){
	yy.push_back(y_min+j*dy);
	gg.push_back((*this)(i,j,z,w));
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(y,yy,gg);
}

double Option4::operator()(int i,int j,double z,double w){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_z)/2. - double(accuracy_z/2);
  option_int = int( (z-z_min) / dz  + half_temp) - (accuracy_z-1)/2;
  option_int = std::max(0,std::min(p-accuracy_z,option_int));
  for(int k=option_int;k<option_int+accuracy_z;k++){
	yy.push_back(z_min+k*dz);
	gg.push_back((*this)(i,j,k,w));
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(z,yy,gg);
}

double Option4::operator()(int i,int j,int k,double w){
  std::vector<double> yy,gg;
  int option_int;
  double half_temp;
  half_temp = double(accuracy_w)/2. - double(accuracy_w/2);
  option_int = int( (w-w_min) / dw  + half_temp) - (accuracy_w-1)/2;
  option_int = std::max(0,std::min(q-accuracy_w,option_int));
  for(int l=option_int;l<option_int+accuracy_w;l++){
	yy.push_back(w_min+l*dw);
	gg.push_back(U[i][j][k][l]);
	// 		    cout << i << " " << yy[i-option_int] << " " << gg[i-option_int] << endl;
  }
  return interpolate(w,yy,gg);
}

double rms(Option4& a,Option4& b)
{
  double sum;
  sum = 0.;
  for(int i=0;i<a.size_x();i++)
	for(int j=0;j<a.size_y();j++)
	  for(int k=0;k<a.size_z();k++)
		for(int l=0;l<a.size_w();l++)
		{
		  sum+=(a[i][j][k][l] - b[i][j][k][l])*(a[i][j][k][l] - b[i][j][k][l]);
		}
  return std::sqrt(sum)/a.size_x()/a.size_y()/a.size_z()/a.size_w();
}

void Print4::calc_grid(void)
{
  if(n==1)dx=0.;
  else dx = (x_max - x_min)/double(n-1);
  if(m==1)dy=0.;
  else dy = (y_max - y_min)/double(m-1);
  if(p==1)dz=0.;
  else dz = (z_max - z_min)/double(p-1);
  if(q==1)dw=0.;
  else dw = (w_max - w_min)/double(q-1);
}

void Print4::print(Option4 &u)
{
  double x,y,z,w;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  #endif
  if(n*m*p*q>1.e6){std::cout << " file would be huge aborting...\n";return;}
  for(int i=0;i<n;i++)
  {
	x = x_min + i*dx;
	for(int j=0;j<m;j++)
	{
	  y = y_min + j*dy;
	  for(int k=0;k<p;k++)
	  {
		z = z_min + k*dz;
		for(int l=0;l<q;l++)
		{
		  w = w_min + l*dw;
		  file <<  x << " " << y << " " << z << " " << w << " " << u(x,y,z,w) << "\n";
		}
		file << "\n";
	  }
	  file << "\n";	
	}
	file << "#########################################################\n\n";
  }
  file.close();
}


void Print4::print_y_z(Option4 &u)
{
  double x,y,z,w;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  #endif
  if(n*m*p*q>1.e6){std::cout << " file would be huge aborting...\n";return;}
  for(int i=0;i<n;i++)
  {
	x = x_min + i*dx;
	for(int l=0;l<q;l++)
	{
	  w = w_min + l*dw;
	  for(int k=0;k<p;k++)
	  {
		z = z_min + k*dz;
		for(int j=0;j<m;j++)
		{
		  y = y_min + j*dy;
		  file <<  x << " " << y << " " << z << " " << w << " " << u(x,y,z,w) << "\n";
		}
		file << "\n";
	  }
	  file << "\n";	
	}
	file << "#########################################################\n\n";
  }
  file.close();
}


void Print4::print_x_y(Option4 &u)
{
  double x,y,z,w;
  calc_grid();
  if(file.is_open()){file.close();}
  file.open(filename.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  file << std::setprecision(precision);
  #ifdef DEBUG
  std::cout << " printing to file " << filename << "\n";
  #endif
  if(n*m*p*q>1.e6){std::cout << " file would be huge aborting...\n";return;}
  for(int l=0;l<q;l++)
  {
	w = w_min + l*dw;
	for(int k=0;k<p;k++)
	{
	  z = z_min + k*dz;  
	  for(int i=0;i<n;i++)
	  {
		x = x_min + i*dx;
		for(int j=0;j<m;j++)
		{
		  y = y_min + j*dy;
		  file <<  x << " " << y << " " << z << " " << w << " " << u(x,y,z,w) << "\n";
		}
		file << "\n";
	  }
	  file << "\n";	
	}
	file << "#########################################################\n\n";
  }
  file.close();
}


void Print4::link_printer_config4(Param_vecs *print_config)
{
  print_config->assignParameter(&q," No of points in W ");
  print_config->assignParameter(&w_min," Minimum value of W ");
  print_config->assignParameter(&w_max," Maximum value of W ");
}

void Print4::link_printer_config4(ParameterMaps *print_config)
{
  print_config->assignParameter(&q,"q"," No of points in W ");
  print_config->assignParameter(&w_min,"w_min"," Minimum value of W ");
  print_config->assignParameter(&w_max,"w_max"," Maximum value of W ");
}

