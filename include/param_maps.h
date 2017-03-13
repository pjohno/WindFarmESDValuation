#ifndef _param_maps_h_included_
#define _param_maps_h_included_

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "pj_data.h"

// #define DEBUG

namespace PJLib {

class ParameterMaps
{
  
  protected:
    
    // create a container for all pointers to the integer parameters
    std::map<std::string,PJ_data2<int> > the_ints;
    std::map<std::string,PJ_data2<int> >::iterator ii;
    // create a container for all pointers to the real parameters
    std::map<std::string,PJ_data2<double> > the_reals;
    std::map<std::string,PJ_data2<double> >::iterator jj;
    // create a container for all pointers to the string parameters
    std::map<std::string,PJ_data2<std::string> > the_strings;
    std::map<std::string,PJ_data2<std::string> >::iterator kk;
    
  public:
    
    //deconstructor
    ~ParameterMaps(){clear();};
    // constructor
    ParameterMaps(){};
    // empty all parameters
    void clear(void);
    // empty all parameters
    void reset(void);
    // set up new parameter - int
    void assignParameter(int *xx,std::string ref,std::string comment);
    // - double
    void assignParameter(double *xx,std::string ref,std::string comment);
    // - string
    void assignParameter(std::string *xx,std::string ref,std::string comment);
    // set up new parameter - int
    void assignParameter(int *xx,std::string ref,std::string comment,int reset_val);
    // - double
    void assignParameter(double *xx,std::string ref,std::string comment,double reset_val);
    // - string
    void assignParameter(std::string *xx,std::string ref,std::string comment,std::string reset_val);
    // set up new parameter - int
    void assignParameter(int *xx,std::string ref,std::string comment,int reset_val,int min,int max);
    // - double
    void assignParameter(double *xx,std::string ref,std::string comment,double reset_val,double min,double max);
    // change parameter value - int
    bool changeParameter(int xx,std::string ref);
    // - double
    bool changeParameter(double xx,std::string ref);
    // - string
    bool changeParameter(std::string xx,std::string ref);
    // erase a parameter from the map
    bool eraseParameter(std::string ref);
    // check whether parameter is there or not
    bool check(std::string ref);
    // return value - int
    int returnInteger(std::string ref);
    // - double
    double returnDouble(std::string ref);
    // - string
    std::string returnString(std::string ref);
    // return pointers - int
    int* returnPtrInteger(std::string ref);
    // - double
    double* returnPtrDouble(std::string ref);
    // - string
    std::string* returnPtrString(std::string ref);
    // return the parameter description
    std::string return_name(std::string ref);
    // check for default
    bool is_default(std::string ref);
    // check for bounds
    bool is_bounded(std::string ref);
    // print all parameters to stream
    virtual void print_all_params(std::ostream *output);
    // - print parameter by refernce
    void print(std::string ref,std::ostream *output);
    // read all parameters to stream
    bool read_all_params(std::string filename);
    // create vectors of all possible parameters
    bool create_vector_refs(std::vector<std::string> &my_ints,std::vector<std::string> &my_doubles,std::vector<std::string> &my_strings);
    
    ParameterMaps& operator=(ParameterMaps& rhs){
      for(ii=the_ints.begin();ii!=the_ints.end();ii++)
      {
	// check if the 
	if(rhs.check((*ii).first))
	{
// 	  std::cout << " change parameter " << (*ii).first << " to " << rhs.returnInteger((*ii).first) << "\n";
	  the_ints[(*ii).first].change_x(rhs.returnInteger((*ii).first));
	}
      }
      for(jj=the_reals.begin();jj!=the_reals.end();jj++)
      {
	// check if the 
	if(rhs.check((*jj).first))
	{
// 	  std::cout << " change parameter " << (*jj).first << " to " << rhs.returnDouble((*jj).first) << "\n";
	  the_reals[(*jj).first].change_x(rhs.returnDouble((*jj).first));
	}
      }
      for(kk=the_strings.begin();kk!=the_strings.end();kk++)
      {
	// check if the 
	if(rhs.check((*kk).first))
	{
// 	  std::cout << " change parameter " << (*kk).first << " to " << rhs.returnString((*kk).first) << "\n";
	  the_strings[(*kk).first].change_x(rhs.returnString((*kk).first));
	}
      };
      return *this;
    };
    
};

}// namespace

#endif
