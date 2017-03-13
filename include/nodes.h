#ifndef _nodes_h_included_
#define _nodes_h_included_

#include "fixed_grid.h"
#include "diff_grid.h"
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Grid
{
  
  private:
	
	// size of dimensions
	int n,m,p;
	bool x_setup,y_setup,z_setup,nodes_setup;
	
  protected:

	// assign default values
	void default_values();
	// setup the grid
	void setup_storage(int n_,int m_,int p_);
	// clear storage
	void clear_nodes(),clear_stencil_Q();
	// Assign storage 
	void assign_stencil_Q(int n_,int p_),assign_nodes(int n_,int m_,int p_);
	// resize elements
	void resize_X(int n_),resize_Y(int m_),resize_Q(int n_,int p_),resize_nodes(int n_,int m_,int p_);
	
  public:
	
	// deconstructor
	~Grid(){clear_nodes();clear_stencil_Q();delete [] grid_X;delete [] grid_Y;}
	// constructor
	Grid();
	// some constant on the grid
	double grid_const;
	// stencils in X and Y
	Stencil *grid_X,*grid_Y;
	// stencil in Q depends on X and Q
	Stencil **grid_Q;
	// node returns value of gg
	Node ***grid;
	// overload access to grid
	Node** operator[](int index){return grid[index];};
	// update storage in grid
	void update(Fixed_Grid &X,Fixed_Grid &Y,Fixed_Grid &Z);
	// reset all values
	Grid& operator=(double a){
	  for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++)grid[i][j][k]=0.;
	  for(int i=0;i<n;i++)grid_X[i]=0.;
	  for(int j=0;j<m;j++)grid_Y[j]=0.;
	  for(int i=0;i<n;i++)for(int k=0;k<p;k++)grid_Q[i][k]=0.;
	  return *this;
	}
};

#endif
