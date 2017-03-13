// implementation of finite difference stencil, not really using this anymore...
// try another change...
#include "diff_grid.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

void Stencil:: assign_adjust(int adjust_)
{
	
	adjust = adjust_;
	
}

void Stencil:: setup_stencil(double aa_,double bb_,double cc_,int adjust_)
{
	
	aa = aa_;
	bb = bb_;
	cc = cc_;
	adjust = adjust_;
	
}

void Stencil:: return_stencil(double *aa_,double *bb_,double *cc_,int *adjust_)
{
	
	*aa_ = aa;
	*bb_ = bb;
	*cc_ = cc;
	*adjust_ = adjust;
	
}

void Stencil:: print_stencil(void)
{
	cout << "grid X " << aa << " " << bb << " " << cc <<endl;
}

int Stencil:: adjuster(void)
{
	return adjust;
}

Stencil:: Stencil()
{
	setup_stencil(0.,0.,0.,0);
}

