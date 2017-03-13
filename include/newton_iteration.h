#ifndef _newton_iteration_h_included_
#define _newton_iteration_h_included_

class Newton_iteration
{
  
  private:

	// counter for iterations
	int counter;
	// initial guess, final value and tolerance
	double x_start,x_curr,dX;
	// store old func value if don't know f'
	bool func_dash_exists;
	double func_old,func_curr,x_old;
	// secant method
	void secant_method(void);
	// newtons method
	void newtons_method(void);
	// input functions
	virtual double func(double x){return x;};
	virtual double func_dash(double x){return 1.;};

  protected:
	
	double tolerance;
	
  public:

	// default  constructor
	Newton_iteration(){counter=0;x_start=1.;dX=0.1;tolerance=1.e-10;};
	// Allow for initial guess
	Newton_iteration(double x_start_){x_start = x_start_;func_dash_exists=true;tolerance=1.e-10;};
	Newton_iteration(double x_start_,double dX_){x_start=x_start_;dX=dX_;func_dash_exists=false;tolerance=1.e-10;};
	// choose iteration method and return x_final
	double iteration(void);
	
};

#endif