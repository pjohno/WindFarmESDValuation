// Take a uniform random distribution and convert to normal distribution
#ifndef _random_ND_h_included_
#define _random_ND_h_included_

#include<cmath>
#include<iostream>
#include<iomanip>

class MTRand;

class Random_ND
{
  // have randoms number been generated
  bool numbers_generated;
  // create a pointer to the random number generator
  MTRand *randoms;
  // the size of the array and the seed
  int size_of_array,seed;
  // this will a pointer to the array of numbers that are normally distributed
  double *normal_distribution;
  // generate the numbers
  void gaussian_distribution(int new_seed);
  // could choose to change the size of the array, will need to regenerate numbers
  void resize(int size);
  // reset seed
  void reset_seed();
  
  public:
	
	// deconstructor
	~Random_ND();
	// default constructor
	Random_ND();
	// assign the size of the normal distribution array
	Random_ND(int size);
	// assign the size and seed of the normal distribution array
	Random_ND(int size,int _seed);
	// overload [] so that it returns the appropriate value in normal_distribution
	inline double operator[](int index){ return normal_distribution[index]; };
	// regenerate new numbers;
	void regenerate();
	// set new seed;
	void set_new_seed(int _seed);
	// return current size of the array;
	inline int size(){return size_of_array;};
	// change size of the array
	void change_size(int size,bool reset);
};

#endif
