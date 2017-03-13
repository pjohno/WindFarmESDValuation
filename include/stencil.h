#ifndef _GENERIC_STENCIL_h_included_
#define _GENERIC_STENCIL_h_included_

#include "fixed_grid.h"
#include "sde.h"

namespace PJLib {
  
  class GenericStencil
  {
    bool check_linked;
    protected:
      // size of current storage
      int n;
      // link the associated grid
      Fixed_Grid *X;
      // link the associated SDE
      SDE *sde;
      // storage for the stencils
      double *aa,*bb,*cc;
      // clear storage
      void clear_storage();
      
    public:
      // constructor
      GenericStencil();
      ~GenericStencil();
      // assign storage
      void assign_storage();
      // link up to a model
      void link(Fixed_Grid *X_,SDE *sde_){X=X_;sde=sde_;};
      // create stencils
      virtual void create_stencils(double tau);
      // public access for the stencils
      inline double a(int index){return aa[index];};
      inline double b(int index){return bb[index];};
      inline double c(int index){return cc[index];};
  };
  
  class CharStencil: public GenericStencil
  {
    
    public:
      
      CharStencil(){};
      ~CharStencil();
      // create stencils
      void create_stencils(int k,double tau);
      
  };
  
  class SecondOrderStencil: public GenericStencil
  {
    
    public:
      
      // create stencils
      void create_stencils(int k,double theta);
      
  };
  
}//namespace

#endif
