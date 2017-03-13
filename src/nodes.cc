#include "diff_grid.h"
#include "nodes.h"
#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

void Grid:: setup_storage(int n_,int m_,int p_)
{
  resize_X(n_);
  resize_Y(m_);
  resize_Q(n_,p_);
  resize_nodes(n_,m_,p_);
  n=n_;
  m=m_;
  p=p_;
  // reset all values to zero
  *this=0.;
}

void Grid::resize_X(int n_)
{
  #ifdef DEBUG
  std::cout << " Grid:: resize x grid \n";
  #endif // DEBUG
  if(n_==n and n_>2)return;
  if(x_setup)delete [] grid_X;
  grid_X = new Stencil[n_];
  #ifdef DEBUG
  std::cout << " new grid size for x \n";
  #endif // DEBUG
  x_setup = true;
}

void Grid::resize_Y(int m_)
{
  #ifdef DEBUG
  std::cout << " Grid:: resize y grid \n";
  #endif // DEBUG
  if(m_==m and m_>2)return;
  if(y_setup)delete [] grid_Y;
  grid_Y = new Stencil[m_];
  #ifdef DEBUG
  std::cout << " new grid size for y \n";
  #endif // DEBUG
  y_setup = true;
}

void Grid::resize_Q(int n_,int p_)
{
  #ifdef DEBUG
  std::cout << " Grid:: resize q grid \n";
  #endif // DEBUG
  if(n_==n and n_>2 and p_==p and p_>2)return;
  if(z_setup)clear_stencil_Q();
  assign_stencil_Q(n_,p_);
  #ifdef DEBUG
  std::cout << " new grid size for q \n";
  #endif // DEBUG
  z_setup = true;
}
//
void Grid::resize_nodes(int n_,int m_,int p_)
{
  #ifdef DEBUG
  std::cout << " Grid:: resize nodes \n";
  #endif // DEBUG
  if(n_==n and n_>2 and m_==m and m_>2 and p_==p and p_>2)return;
  if(nodes_setup)clear_nodes();
  assign_nodes(n_,m_,p_);
  #ifdef DEBUG
  std::cout << " new grid size for nodes \n";
  #endif // DEBUG
  nodes_setup = true;
}
// clear grid nodes
void Grid::clear_nodes()
{
  for(int i=0;i<n;i++){
	for(int j=0;j<m;j++)delete [] grid[i][j];
	delete [] grid[i];
  }
  delete [] grid;
}
// clear Q stencil
void Grid::clear_stencil_Q()
{
  for(int i=0;i<n;i++)delete [] grid_Q[i];
  delete [] grid_Q;
}
// 
void Grid::assign_stencil_Q(int n_,int p_)
{
  grid_Q = new Stencil*[n_];
  for(int i=0;i<n_;i++)grid_Q[i] = new Stencil[p_];
}
void Grid::assign_nodes(int n_,int m_,int p_)
{
  grid = new Node**[n_];
  for(int i=0;i<n_;i++)
  {
	grid[i]= new Node*[m_];
	for(int j=0;j<m_;j++){grid[i][j]= new Node[p_];}
  }
}
//
void Grid::default_values()
{
  x_setup=false;y_setup=false;z_setup=false;nodes_setup=false;
  n=2;m=2;p=2;
}

void Grid::update(Fixed_Grid &X,Fixed_Grid &Y,Fixed_Grid &Z)
{
  setup_storage(X.size(),Y.size(),Z.size());
}

Grid::Grid()
{
  default_values();
  setup_storage(n,m,p);
}
