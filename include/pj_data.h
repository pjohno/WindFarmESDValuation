// pj data holds a pointer to the parameter, its name, along with a default vaue and min/max values.
#ifndef _pj_data_h_included_
#define _pj_data_h_included_

#include <iostream>
#include <string>

// #define DEBUG

namespace PJLib {

template <class param>
class PJ_data {
  // desccription of data
  std::string name;
  // pointer to data
  param* ptr_x;
	
  public:
	// Constructor
	PJ_data(){ptr_x=NULL;};
	PJ_data(std::string name_,param* ptr_x_){set_pointer(ptr_x_);set_name(name_);};
	// assign pointer
	void set_pointer(param* ptr_x_){ptr_x=ptr_x_;};
	// change value
	virtual void change_x(param x_) {*ptr_x=x_;};
	// set name
	void set_name(std::string name_){name=name_;};
	// print value
	void print_data(std::ostream *output){(*output) << *ptr_x;};
	// print name
	void print_name(std::ostream *output){(*output) << name;};
	// print name and value
	virtual void print(std::ostream *output){print_name(output);print_data(output);};
	// return name
	std::string return_name(void){return name;};
	// return position
	param* return_pointer(void){return ptr_x;};
	// return value
	param return_value(void){return *ptr_x;};
};

template <class param2>
class PJ_data2:public PJ_data<param2>
{
  // don't allow assignment of value outside min/max
  bool constraint,default_available;
  // include a default value, minimum and a maximum
  param2 default_x,min_x,max_x;
  public:
	// constructors
	PJ_data2():PJ_data<param2>(){default_available=false;};
	PJ_data2(std::string name_,param2* ptr_x_):PJ_data<param2>(name_,ptr_x_){constraint = false;default_available=false;};
	PJ_data2(std::string name_,param2* ptr_x_,param2 x_):PJ_data<param2>(name_,ptr_x_){set_default(x_);constraint = false;};
	PJ_data2(std::string name_,param2* ptr_x_,param2 x_,param2 min,param2 max):PJ_data<param2>(name_,ptr_x_){set_default(x_);set_min(min);set_max(max);};
	// check on type of parameter
	bool is_default(){return default_available;};
	bool is_bounded(){return constraint;};
	// reset value to default
	void reset(){if(default_available)change_x(default_x);};
	// set default, min and max
	void set_default(param2 x_){default_x=x_;default_available=true;};
	void set_min(param2 x){min_x=x;};
	void set_max(param2 x){max_x=x;constraint = true;};
	// return min,max and default
	param2 return_default(){return default_x;};
	param2 return_min(){return min_x;};
	param2 return_max(){return max_x;};
	// print default, min and max
	void print_default(std::ostream *output){(*output) << " " << default_x;};
	void print_min(std::ostream *output){(*output) << " " << min_x;};
	void print_max(std::ostream *output){(*output) << " " << max_x;};
	void print(std::ostream *output){this->print_name(output);this->print_data(output);print_default(output);print_min(output);print_max(output);};
	// overload change_x to check against bounds
	void change_x(param2 x_){
	  if(constraint)
	  {
		#ifdef DEBUG
		std::cout << " &x " << this->return_pointer() <<  " x " 
		<< *this->return_pointer() << " : " << x_ << "\n";
		#endif
		x_ = std::max(x_,min_x);
		x_ = std::min(x_,max_x);
		*this->return_pointer()=x_;
	  }
	  else 
	  {
		#ifdef DEBUG
		std::cout << " &x " << *this->return_pointer()<< " x " << this->return_pointer() << "\n";
		#endif
		*this->return_pointer()=x_;
	  }
	};
	
};

} // namespace

#endif
