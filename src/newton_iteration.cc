#include "newton_iteration.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

double Newton_iteration::iteration(void)
{
  if(func_dash_exists)newtons_method();
  else secant_method();
  
  return x_curr;
}
	// secant method
void Newton_iteration::secant_method(void)
{
  x_old = x_start;
  func_old = func(x_old);
  x_curr = x_old + dX;
  for(counter = 0;counter<100;counter++){
	func_curr = func(x_curr);
	dX = - func_curr*dX/(func_curr-func_old);
	x_old = x_curr;
	func_old = func_curr;
	x_curr = x_old + dX;
	if(std::abs(func_curr)<tolerance)break;
  }
  if(counter>=100)std::cout << "Not converged";
  
}
	// newtons method
void Newton_iteration::newtons_method(void)
{
  x_curr = x_start;
  for(counter = 0;counter<100;counter++){
	func_curr = func(x_curr);
	x_curr = x_curr - func_curr/func_dash(x_curr);
	//std::cout << x_curr << " " << func_curr << "\n";
	if(std::abs(func_curr)<tolerance)break;
  }
  if(counter>=100)std::cout << "Not converged";
  
}