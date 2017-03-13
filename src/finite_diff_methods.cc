/* Set up a class for the option value */

#include "finite_diff_methods.h"
#include "fixed_grid.h"
#include "nodes.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

void Methods::setup_storage(int n_,int m_,int p_)
{
  resize_X(n_);
  resize_Y(m_);
  p = p_;
};

void Methods::resize_X(int n_)
{
  #ifdef DEBUG
  std::cout << " Methods :: resize x grid \n";
  #endif // DEBUG
  if(n_==n and n_>2)return;
  if(x_setup)clear_X();
  assign_X(n_);
   #ifdef DEBUG
  std::cout << " new grid size for x \n";
  #endif // DEBUG
  x_setup = true;
  n=n_;
}

void Methods::resize_Y(int m_)
{
  #ifdef DEBUG
  std::cout << " Methods :: resize y grid \n";
  #endif // DEBUG
  if(m_==m and m_>2)return;
  if(y_setup)clear_Y();
  assign_Y(m_);
  #ifdef DEBUG
  std::cout << " new grid size for y \n";
  #endif // DEBUG
  y_setup = true;
  m=m_;
}


void Methods::clear_X()
{
  delete [] a_X;
  delete [] b_X;
  delete [] c_X;
  delete [] d_X;
  delete [] w_X;
}

void Methods::clear_Y()
{
  delete [] a_Y;
  delete [] b_Y;
  delete [] c_Y;
  delete [] d_Y;
  delete [] w_Y;
}

void Methods::assign_X(int n_)
{
  a_X = new double[n_];
  b_X = new double[n_];
  c_X = new double[n_];
  d_X = new double[n_];
  w_X = new double[n_];
}

void Methods::assign_Y(int m_)
{
  a_Y = new double[m_];
  b_Y = new double[m_];
  c_Y = new double[m_];
  d_Y = new double[m_];
  w_Y = new double[m_];
}

void Methods::default_values()
{
  x_setup=false;y_setup=false;
  n=2;m=2;p=2;
}

Methods:: Methods()
{
  default_values();
  setup_storage(n,m,p);
}
// constructor with initial setup for storage
Methods:: Methods(Fixed_Grid &X,Fixed_Grid &Y,Fixed_Grid &Z)
{
  default_values();
  update(X,Y,Z);
}
// update vector sizes
void Methods:: update(Fixed_Grid &X,Fixed_Grid &Y,Fixed_Grid &Z)
{
  setup_storage(X.size(),Y.size(),Z.size());
}

void Methods:: tridag_solver(vector<double> *w,vector<double> aa,
			     vector<double>bb,vector<double>cc,
				      vector<double> dd)
{
  /*
  ! tridiagonal solver uses Gaussian elimination
  ! note that no checks are made for zero diagonal elements
  !  A(1) B(1) C(1)                     D(1)
  !  A(2) B(2) C(2)                     D(2)
  !       A(3) B(3) C(3)                D(3)
  !               ----------            ||
  !        ----- A(i) B(i) C(i) ---     D(i)
  !               ----------            ||
  !         A(n-2) B(n-2) C(n-2)        D(n-2)
  !                A(n-1) B(n-1) C(n-1) D(n-1)
  !                A(n)   B(n)   C(n)   D(n)
  */
  
  int n = aa.size();
  
  bb[1]=bb[1]-bb[0]*aa[1]/aa[0]; // elim A(2]
  cc[1]=cc[1]-cc[0]*aa[1]/aa[0];
  dd[1]=dd[1]-dd[0]*aa[1]/aa[0];
  for(int j = 2;j<=n-2;j++)
  {
    bb[j]=bb[j]-cc[j-1]*aa[j]/bb[j-1];  // elim A(3]-A(n-1]
    dd[j]=dd[j]-dd[j-1]*aa[j]/bb[j-1];
  }
  bb[n-1]=bb[n-1]-cc[n-3]*aa[n-1]/bb[n-3]; // elim A(n]
  dd[n-1]=dd[n-1]-dd[n-3]*aa[n-1]/bb[n-3];
  cc[n-1]=cc[n-1]-cc[n-2]*bb[n-1]/bb[n-2]; // elim B(n]
  dd[n-1]=dd[n-1]-dd[n-2]*bb[n-1]/bb[n-2];
  // calculate solution 
  (*w)[n-1]=dd[n-1]/cc[n-1];
  for(int j=n-2;j>=1;j--)
  {
    (*w)[j]=(dd[j]-cc[j]*(*w)[j+1])/bb[j];
  }
  (*w)[0]=(dd[0]-cc[0]*(*w)[2]-bb[0]*(*w)[1])/aa[0];
}

void Methods::tridag_solver_X()
{ 
  /*
  ! tridiagonal solver uses Gaussian elimination
  ! note that no checks are made for zero diagonal elements
  !  A(1) B(1) C(1)                     D(1)
  !  A(2) B(2) C(2)                     D(2)
  !       A(3) B(3) C(3)                D(3)
  !               ----------            ||
  !        ----- A(i) B(i) C(i) ---     D(i)
  !               ----------            ||
  !         A(n-2) B(n-2) C(n-2)        D(n-2)
  !                A(n-1) B(n-1) C(n-1) D(n-1)
  !                A(n)   B(n)   C(n)   D(n)
  */
  
  b_X[1]=b_X[1]-b_X[0]*a_X[1]/a_X[0]; // elim A(2]
  c_X[1]=c_X[1]-c_X[0]*a_X[1]/a_X[0];
  d_X[1]=d_X[1]-d_X[0]*a_X[1]/a_X[0];
  for(int j = 2;j<=n-2;j++)
  {
    b_X[j]=b_X[j]-c_X[j-1]*a_X[j]/b_X[j-1];  // elim A(3]-A(n-1]
    d_X[j]=d_X[j]-d_X[j-1]*a_X[j]/b_X[j-1];
  }
  b_X[n-1]=b_X[n-1]-c_X[n-3]*a_X[n-1]/b_X[n-3]; // elim A(n]
  d_X[n-1]=d_X[n-1]-d_X[n-3]*a_X[n-1]/b_X[n-3];
  c_X[n-1]=c_X[n-1]-c_X[n-2]*b_X[n-1]/b_X[n-2]; // elim B(n]
  d_X[n-1]=d_X[n-1]-d_X[n-2]*b_X[n-1]/b_X[n-2];
  // calculate solution 
  w_X[n-1]=d_X[n-1]/c_X[n-1];
  for(int j=n-2;j>=1;j--)
  {
    w_X[j]=(d_X[j]-c_X[j]*w_X[j+1])/b_X[j];
  }
  w_X[0]=(d_X[0]-c_X[0]*w_X[2]-b_X[0]*w_X[1])/a_X[0];
}


void Methods::tridag_solver_Y()
{
  /*
  ! tridiagonal solver uses Gaussian elimination
  ! note that no checks are made for zero diagonal elements
  !  A(1) B(1) C(1)                     D(1)
  !  A(2) B(2) C(2)                     D(2)
  !       A(3) B(3) C(3)                D(3)
  !               ----------            ||
  !        ----- A(i) B(i) C(i) ---     D(i)
  !               ----------            ||
  !         A(n-2) B(n-2) C(n-2)        D(n-2)
  !                A(n-1) B(n-1) C(n-1) D(n-1)
  !                A(n)   B(n)   C(n)   D(n)
  */
  
  b_Y[1]=b_Y[1]-b_Y[0]*a_Y[1]/a_Y[0]; // elim A(2]
  c_Y[1]=c_Y[1]-c_Y[0]*a_Y[1]/a_Y[0];
  d_Y[1]=d_Y[1]-d_Y[0]*a_Y[1]/a_Y[0];
  for(int j = 2;j<=m-2;j++)
  {
    b_Y[j]=b_Y[j]-c_Y[j-1]*a_Y[j]/b_Y[j-1];  // elim A(3]-A(n-1]
    d_Y[j]=d_Y[j]-d_Y[j-1]*a_Y[j]/b_Y[j-1];
  }
  b_Y[m-1]=b_Y[m-1]-c_Y[m-3]*a_Y[m-1]/b_Y[m-3]; // elim A(m]
  d_Y[m-1]=d_Y[m-1]-d_Y[m-3]*a_Y[m-1]/b_Y[m-3];
  c_Y[m-1]=c_Y[m-1]-c_Y[m-2]*b_Y[m-1]/b_Y[m-2]; // elim B(m]
  d_Y[m-1]=d_Y[m-1]-d_Y[m-2]*b_Y[m-1]/b_Y[m-2];
  // calculate solutiom 
  w_Y[m-1]=d_Y[m-1]/c_Y[m-1];
  for(int j=m-2;j>=1;j--)
  {
    w_Y[j]=(d_Y[j]-c_Y[j]*w_Y[j+1])/b_Y[j];
  }
  w_Y[0]=(d_Y[0]-c_Y[0]*w_Y[2]-b_Y[0]*w_Y[1])/a_Y[0];
}

void Methods::ADI_imp_Y(double dt,double ***w, Grid *GG)
{

	double r;
	double aa_,bb_,cc_;
	double dd_,ee_,ff_;
	int adjust_;
	
	bool comment=false;
	
	for(int k=0;k<p;k++)
	{
		for(int i=0;i<n;i++)
		{
			
			adjust_=(*GG).grid_Q[i][k]();
			if(adjust_==0)
			{
				r = 1.;
			}
			else
			{
				r = 0.5;
			}
			(*GG).grid_Y[0].return_stencil(&aa_,&bb_,&cc_,&adjust_);
			
			a_Y[0] = r*(aa_ + 2./dt);
			b_Y[0] = r*bb_ ;
			c_Y[0] = r*cc_ ;
			d_Y[0] = (*GG)[i][0][k].ZZ();
			for(int j=1;j<m-1;j++)
			{
				(*GG).grid_Y[j].return_stencil(&aa_,&bb_,&cc_,&adjust_);
				a_Y[j] = r*aa_ ;
				b_Y[j] = r*(bb_ + 2./dt);
				c_Y[j] = r*cc_ ;
				d_Y[j] = (*GG)[i][j][k].ZZ();
			}
			(*GG).grid_Y[m-1].return_stencil(&aa_,&bb_,&cc_,&adjust_);
			a_Y[m-1] = r*aa_  ;
			b_Y[m-1] = r*bb_ ;
			c_Y[m-1] = r*(cc_ + 2./dt) ;
			d_Y[m-1] = (*GG)[i][m-1][k].ZZ();
			
			tridag_solver_Y();
			
				if(comment)
				{cout << i << " i:k " << k << " dt "<<dt<<" aa "<<a_Y[m/2]<< " bb "<<b_Y[m/2]<< " cc "<<c_Y[m/2]<< " dd "<<d_Y[m/2]<< " ee " << ee_ << " ww "<<w_Y[m/2]<<endl;}
				for(int j=0;j<m;j++)
				{
					w[i][j][k] = w_Y[j];
					
				}
				
				
		}
		if(comment){cin >> adjust_;}
	}
	
}

void Methods::ADI_imp_X(double dt,double ***w, Grid *GG)
{
	
	double aa_,bb_,cc_,flip;
	double dd_,ee_,ff_;
	int adjust_;
	
	bool comment=false;
	
	for(int k=0;k<p;k++)
	{
		for(int j=0;j<m;j++)
		{
			
			(*GG).grid_X[0](&aa_,&bb_,&cc_,&adjust_);
			(*GG).grid_Q[0][k](&dd_,&ee_,&ff_,&adjust_);
			
			if(adjust_==0)
			{
				flip = 2./dt;
			}
			else
			{
				flip = 2.*ee_;
			}
			a_X[0] = flip + aa_ ;
			b_X[0] = bb_;
			c_X[0] = cc_ ;
			d_X[0] = (*GG)[0][j][k].ZZ();
			for(int i=1;i<n-1;i++)
			{
				(*GG).grid_X[i](&aa_,&bb_,&cc_,&adjust_);
				(*GG).grid_Q[i][k](&dd_,&ee_,&ff_,&adjust_);
				if(adjust_==0)
				{
					flip = 2./dt;
				}
				else
				{
					flip = 2.*ee_;
				}
				a_X[i] = aa_ ;
				b_X[i] = flip + bb_;
				c_X[i] = cc_ ;
				d_X[i] = (*GG)[i][j][k].ZZ();
			}
			(*GG).grid_X[n-1](&aa_,&bb_,&cc_,&adjust_);
			(*GG).grid_Q[n-1][k](&dd_,&ee_,&ff_,&adjust_);
			if(adjust_==0)
			{
				flip = 2./dt;
			}
			else
			{
				flip = 2.*ee_;
			}
			a_X[n-1] = aa_ ;
			b_X[n-1] = bb_ ;
			c_X[n-1] = flip + cc_  ;
			d_X[n-1] = (*GG)[n-1][j][k].ZZ();
			
			tridag_solver_X();
			
			if(comment and j==m/2){cout << "Solve at j="<<j<<" k="<<k<<endl;}
			for(int i=0;i<n;i++)
			{
				if(comment and j==m/2)
				{cout << i << " aa "<<a_X[i]<< " bb "<<b_X[i]<< " cc "<<c_X[i]<< " dd "<<d_X[i]<< " ww "<<w_X[i]<<endl;}
				w[i][j][k] = w_X[i];
				
			}
			
			if(comment and j==m/2){cin >> adjust_;}
			
		}
	}
	
}


double Methods::interpolate(double x,vector<double> y,vector<double> g)
{
  
  double temp=0.;
  //check vectors
  if(y.size()!=g.size()){return 0.;}
  //
  // interpolate = sum{i=1,n} g(i) Pi{j=1,n,j/=i} (x - x_j) / (x_i - x_j)
  //
  for(int i=0;i<y.size();i++)
  {
    
    double int_temp;
    int_temp = g[i];
    for(int j=0;j<y.size();j++)
    {
      if(j==i){continue;}
      int_temp *= ( x - y[j] )/( y[i] - y[j] );
    }
    
    temp += int_temp;
    
  }
  
  return temp;
  
}

void Methods::polynomial(vector<double> y,vector<double> g,vector<double> *coefficients)
{
	
	(*coefficients).clear();
	
// 	cout << "\n Generate coefficients... \n";
	
// 	for(int i=0;i<y.size();i++)
// 	{
// 		cout << y[i] << " " << g[i] << endl;
// 	}
// 	
// 	cout << "\n \n";
	
	double temp=0.;
	//check vectors
	if(y.size()!=4){return;}
	//
	// interpolate = sum{i=1,n} g(i) Pi{j=1,n,j/=i} (x - x_j) / (x_i - x_j)
	//
	//// coefficient[0] = sum{i=1,n} g(i) Pi{j=1,n,j/=i} -x_j / (x_i - x_j)
	//
	for(int i=0;i<y.size();i++)
	{
		
		double int_temp;
		int_temp = g[i];
		for(int j=0;j<y.size();j++)
		{
			if(j==i){continue;}
			int_temp *= - y[j] /( y[i] - y[j] );
		}
		
		temp += int_temp;
		
	}
	
	(*coefficients).push_back(temp);
	
	temp=0.;
	// 	coefficient[1] = sum{i=1,n} g(i) Pi{j=1,n,j/=i} ( sum{k=1,n,k/=i} sum{l=k,n,l/=i,l/=k} (-x_k)(-x_l)  / (x_i - x_j) )
	// 	
	for(int i=0;i<y.size();i++)
	{
		
		double int_temp,sum_temp = 0.;
		int_temp = g[i];
		for(int j=0;j<y.size();j++)
		{
			if(j==i){continue;}
			int_temp *= 1. /( y[i] - y[j] );
		}
		
		for(int k=0;k<y.size();k++)
		{
			if(k==i){continue;}
			for(int l=k;l<y.size();l++)
			{
				if(l==i){continue;}	
				if(l==k){continue;}	
				
				
				sum_temp += y[k] * y[l] ;
			}
		}
		
		int_temp *= sum_temp;
		
		temp += int_temp;
		
	}
	
	(*coefficients).push_back(temp);
	
	temp=0.;
	
	// coefficient[2] = sum{i=1,n} g(i) Pi{j=1,n,j/=i} 1 / (x_i - x_j) sum{j=1,n,j/=i} -x_j
	//
	for(int i=0;i<y.size();i++)
	{
		
		double int_temp,sum_temp=0.;
		int_temp = g[i];
		for(int j=0;j<y.size();j++)
		{
			if(j==i){continue;}
			int_temp *= 1. /( y[i] - y[j] );
		}
		
		for(int k=0;k<y.size();k++)
		{
			if(k==i){continue;}
			sum_temp += - y[k] ;
		}
		
		temp += int_temp*sum_temp;
		
	}
	
	(*coefficients).push_back(temp);
	
	temp=0.;
	// coefficient[3] = sum{i=1,n} g(i) Pi{j=1,n,j/=i}  1 / (x_i - x_j)
	//
	for(int i=0;i<y.size();i++)
	{
		
		double int_temp;
		int_temp = g[i];
		for(int j=0;j<y.size();j++)
		{
			if(j==i){continue;}
			int_temp *= 1. /( y[i] - y[j] );
		}
		
		temp += int_temp;
		
	}
	
	(*coefficients).push_back(temp);
	
	
}

double Methods:: poly_eval(double x,vector<double> a_n)
{
	int n = a_n.size();
	double temp = a_n[n-1];

	for(int i=n-2;i>=0;i--)
	{
		temp = a_n[i] + temp * x;
	}
	return temp;
	
}

double Methods:: poly_max(double x_min,double x_max,vector<double> a_n)
{
	
	if(a_n.size()!=4){cout << "\ninconsistent input poly_max\n";return 0.;}
	
	double temp,temp_pm;
	
	temp = 4.*a_n[2]*a_n[2] - 12.*a_n[1]*a_n[3];
	if(temp<0.){cout << "\n no maximum exists... \n";return 0.;};
// 	cout << -2.*a_n[2] << " " << pow(temp,0.5) << "  " << 2.*3.*a_n[3] <<endl;
	temp_pm = ( -2.*a_n[2] + pow(temp,0.5) ) / 2./3./a_n[3];
	if(temp_pm>x_min and temp_pm<x_max){return temp_pm;};
	temp_pm = ( -2.*a_n[2] - pow(temp,0.5) ) / 2./3./a_n[3];
	if(temp_pm>x_min and temp_pm<x_max){return temp_pm;};
	
}	

double Methods:: richardson_extrap(double y_1,double y_2,int n_1,int n_2,int p)
{
	return ( pow(double(n_2),p) * y_2 - pow(double(n_1),p) * y_1 ) / ( pow(double(n_2),p) - pow(double(n_1),p) );
}
