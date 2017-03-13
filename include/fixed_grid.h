#ifndef _fixed_grid_h_included_
#define _fixed_grid_h_included_

#include <iostream>
#include <cmath>

class Fixed_Grid{

  // check the grid is setup
  bool grid_is_setup;
  // number of nodes in the grid
  int n;
  // size of the grid and variable storage 
  double x_min,x_max,dX,*X;
  // reset the value of dX
  void set_dX(int n_){if(n_==1)dX=0.;else dX=( x_max - x_min )/double(n_-1);};
  // resize grid
  void resize_grid(int n_)
  {
    // std::cout << " resize grid \n";
	if(n_==n and n_>2)return;
    if(grid_is_setup)delete [] X;
	X = new double[n_];
    //    std::cout << " new grid size \n";
    grid_is_setup = true;
  };
  // input constant grid
  void setup_grid(int n_)
  {
    //    std::cout << " recalculate grid \n" ;
    if(!grid_is_setup)return;
    if(n_<1){std::cout << " can't have grid size < 2, using default\n " ; return;}
    set_dX(n_);
    for(int i=0;i<n_;i++)
      {
	X[i] = x_min + i*dX;
      }
  };
  // default values, a two point grid over [0:1]
  void default_values(void)
  {
    grid_is_setup=false;X=NULL;
  };
  
public:
  
  // default constructor
  Fixed_Grid(){default_values();setup_x(2,0.,1.);};
  // specify grid parameters in constructor
  Fixed_Grid(int n_,double min,double max)
  {
    default_values();
    setup_x(n_,min,max);
  }
  // delete array x on exit
  ~Fixed_Grid(){delete [] X;};
  // overload = so that it copies all values from Y
  Fixed_Grid& operator=(Fixed_Grid& Y)
  {
	setup_x(Y.size(),Y.min(),Y.max());
	return *this;
  }
  // overload [] so that it returns the value in X[]
  double operator[](int index)
  {
    //    if(grid_is_setup and index<n and index>=0)return X[index];
    //    else{std::cout << " error returning value Grid[]\n"; abort();}
    return X[index];
  };
  void setup_x(int n_,double min,double max)
  {
	change_n(n_);
	set_range(min,max);
  }
  // change the number of nodes in the grid
  void change_n(int n_){
    resize_grid(n_);
    setup_grid(n_);
    n=n_;
  };
  // change the minimum value of x
  void set_range(double min,double max)
  {
    x_min = min;
	x_max = max;
    setup_grid(n);
  };
  // return current size of the array
  int size(void){return n;};
  double dx(void){return dX;};
  double min(void){return x_min;};
  double max(void){return x_max;};
  double last(void){return X[n-1];};
  
};

#endif

