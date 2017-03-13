#ifndef _no_price_h_included_
#define _no_price_h_included_
#include <param_vecs.h>
#include <storage.h>

class No_Price: public Storage
{
  
	
  
  public:
	
	No_Price(PJLib::Param_vecs *my_parameters,bool progress_bar);
	
	virtual double Power_Stack(double yy);
	
  
};

class Constant_Contract: public Storage
{
  
	
  
  public:
	
	Constant_Contract(PJLib::Param_vecs *my_parameters,bool progress_bar);
	
	virtual double payments(double yy);
	
  
};

#endif
