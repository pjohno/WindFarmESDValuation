// this stencil stores the three piont central differencing stencil
#include "fixed_grid.h"
#include "stencil.h"

using namespace PJLib;

// constructor
GenericStencil::GenericStencil()
{
  // current storage is zero
  n=0;
}

GenericStencil::~GenericStencil()
{
  clear_storage();
//   std::cout << " Stencil cleared\n";
}

// clear storage
void GenericStencil::clear_storage()
{
  if(n!=0)
  {
    // delete storage from stencil
    delete [] aa;
    delete [] bb;
    delete [] cc;
    // reset size of storage to 0
    n=0;
  }
}

// assign storage
void GenericStencil::assign_storage()
{
  // if storage does not need to reassigned don't bother
  if(n==X->size())return;
  // clear current storage
  clear_storage();
  // set up stencil storage to the same size as the grid
  n=X->size();
  aa = new double[n];
  bb = new double[n];
  cc = new double[n];
}

// create stencils
void GenericStencil::create_stencils(double tau)
{
  // boundary conditions
  double a_coeff,b_coeff;
  for(int i=0;i<1;i++)
  {
    a_coeff = sde->alpha((*X)[i],tau)/((*X)[i+1]-(*X)[i]);
    aa[i] = -a_coeff;
    bb[i] = a_coeff;
    cc[i] = 0.;
  }
  // midpoints
  for(int i=1;i<X->size()-1;i++)
  {
    a_coeff = sde->alpha((*X)[i],tau);
    b_coeff = sde->beta((*X)[i],tau);
    b_coeff = 0.5*b_coeff*b_coeff/((*X)[i+1]-(*X)[i])/((*X)[i]-(*X)[i-1]);
    aa[i] = b_coeff;
    bb[i] = -2*b_coeff;
    cc[i] = b_coeff;
    if(a_coeff/((*X)[i+1]-(*X)[i-1])>b_coeff)
    {
      bb[i] = bb[i]-a_coeff/((*X)[i+1]-(*X)[i]);
      cc[i] = cc[i]+a_coeff/((*X)[i+1]-(*X)[i]);
    }
    else if(a_coeff/((*X)[i+1]-(*X)[i-1])<-b_coeff)
    {
      aa[i] = aa[i]-a_coeff/((*X)[i]-(*X)[i-1]);
      bb[i] = bb[i]+a_coeff/((*X)[i]-(*X)[i-1]);
    }
    else
    {
      aa[i] = aa[i]-a_coeff/((*X)[i+1]-(*X)[i-1]);
      cc[i] = cc[i]+a_coeff/((*X)[i+1]-(*X)[i-1]);
    }
      
  }
  // boundary condition
  for(int i=X->size()-1;i<X->size();i++)
  {
    a_coeff = sde->alpha((*X)[i],tau)/((*X)[i]-(*X)[i-1]);
    aa[i] = 0.;
    bb[i] = -a_coeff;
    cc[i] = a_coeff;
  }
  
}

// create stencils
void CharStencil::create_stencils(int k,double tau)
{
  for(int i=0;i<n;i++)
  {
    bb[i] = sde->alpha(i,tau);
  }
}

CharStencil::~CharStencil()
{
//   std::cout << " Char Stencil cleared\n";
}


// create stencils
void SecondOrderStencil::create_stencils(int k,double theta)
{
  // boundary conditions
  double a_coeff,b_coeff;
  for(int i=0;i<1;i++)
  {
    a_coeff = sde->alpha((*X)[i],k,theta)/((*X)[i+1]-(*X)[i]);
    aa[i] = -a_coeff;
    bb[i] = a_coeff;
    cc[i] = 0.;
  }
  // midpoints
  for(int i=1;i<X->size()-1;i++)
  {
    a_coeff = sde->alpha((*X)[i],k,theta);
    b_coeff = sde->beta((*X)[i],k,theta);
    b_coeff = 0.5*b_coeff*b_coeff/((*X)[i+1]-(*X)[i])/((*X)[i]-(*X)[i-1]);
    aa[i] = b_coeff;
    bb[i] = -2*b_coeff;
    cc[i] = b_coeff;
    aa[i] = aa[i]-a_coeff/((*X)[i+1]-(*X)[i-1]);
    cc[i] = cc[i]+a_coeff/((*X)[i+1]-(*X)[i-1]);
  }
  // boundary condition
  for(int i=X->size()-1;i<X->size();i++)
  {
    a_coeff = sde->alpha((*X)[i],k,theta)/((*X)[i]-(*X)[i-1]);
    aa[i] = 0.;
    bb[i] = -a_coeff;
    cc[i] = a_coeff;
  }
  
}
