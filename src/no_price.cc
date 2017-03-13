#include "param_vecs.h"
#include "storage.h"
#include "no_price.h"
#include <iostream>

No_Price::No_Price(Param_vecs *my_parameters,bool progress_bar):Storage(my_parameters,progress_bar)
{

}
	
double No_Price::Power_Stack(double yy)
{
  return 1.;
}

Constant_Contract::Constant_Contract(Param_vecs *my_parameters,bool progress_bar):Storage(my_parameters,progress_bar)
{
// 	std::cout << " I pay constant contracts... ";
}

// wind farm receives payment at average price of electricity
double Constant_Contract::payments(double yy)
{
  return contract_length * average_Y * power_max * C_commit;
}
