
#include<cmath>
#include<iostream>
#include<iomanip>
#include "mtrand.h"
#include "random_nd.h"

// deconstructor
Random_ND::~Random_ND()
{
  delete randoms;
  delete [] normal_distribution;
}
// default constructor
Random_ND::Random_ND()
{
  seed = 0;
  randoms = NULL;
  normal_distribution = NULL;
  resize(2);
  randoms = new MTRand(seed);
  gaussian_distribution(seed);
}
// assign the size of the normal distribution array
Random_ND::Random_ND(int size)
{
  seed = 0;
  randoms = NULL;
  normal_distribution = NULL;
  resize(size);
  randoms = new MTRand(seed);
  gaussian_distribution(seed);
}
// assign the size and seed of the normal distribution array
Random_ND::Random_ND(int size,int _seed)
{
  seed = _seed;
  randoms = NULL;
  normal_distribution = NULL;
  resize(size);
  randoms = new MTRand(seed);
  gaussian_distribution(seed);
}
// could choose to change the size of the array, will need to regenerate numbers
void Random_ND::resize(int size)
{
  delete [] normal_distribution;
  if(size%2 == 1)
  {
	size+=1;
  }
  size_of_array = size;
  normal_distribution = new double[size_of_array];
}

// generate the numbers
void Random_ND::gaussian_distribution(int new_seed)
{
  bool debug = false;
  if(debug)std::cout << " gas dev ... \n";
  
  int i=0;
  double v1,v2,rsq;
  int n=size_of_array;
  if(n==0)return;
  if(new_seed!=seed)
  {
	std::cout << "\n reset seed \n" ;
	delete randoms;
	randoms = new MTRand(new_seed);
  }
  do
  {
	
	v1 = (*randoms)();
	v2 = (*randoms)();
	v1 = 2.*v1-1.;
	v2 = 2.*v2-1.;
	rsq = v1*v1+v2*v2;
	
	if(rsq<1.)
	{
	  rsq=std::pow(-2.*log(rsq)/rsq , 0.5);
	  normal_distribution[i] = v1*rsq;
	  normal_distribution[i+1] = v2*rsq;
	  i = i + 2;
	}
	
  }while(i<n);
  seed = new_seed;
}
// regenerate new numbers;
void Random_ND::regenerate()
{
  gaussian_distribution(seed);
}
// set new seed;
void Random_ND::set_new_seed(int _seed)
{
  gaussian_distribution(_seed);
  seed = _seed;
}
// change size of array
void Random_ND::change_size(int size,bool reset)
{
  resize(size);
  if(reset)reset_seed();
  else gaussian_distribution(seed);
}
// reset the seed 
void Random_ND::reset_seed()
{
  int temp_seed = seed;
  seed += 1;
  gaussian_distribution(temp_seed);
}
