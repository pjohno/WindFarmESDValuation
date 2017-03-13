
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <string.h>
#include "pj_data.h"
#include "param_vecs.h"

// #define DEBUG
using namespace PJLib;

void Param_vecs::clear_vecs(void)
{
	
	the_ints.clear();the_reals.clear();the_strings.clear();
	
}

void Param_vecs::assignParameter(int *xx,std::string name)
{
	
	
	PJ_data<int> integer_parameter;
	integer_parameter.set_name(name);
	integer_parameter.set_pointer(xx);
	
	the_ints.push_back(integer_parameter);
	
}


void Param_vecs::assignParameter(double *xx,std::string name)
{
	
	PJ_data<double> double_parameter;
	double_parameter.set_name(name);
	double_parameter.set_pointer(xx);
	
	the_reals.push_back(double_parameter);
	
}


void  Param_vecs::assignParameter(std::string *xx,std::string name)
{
	
	PJ_data<std::string> string_parameter;
	string_parameter.set_name(name);
	string_parameter.set_pointer(xx);
	
	the_strings.push_back(string_parameter);
	
}

void Param_vecs::changeParameter(int xx,int parameter_ref)
{
	
	if(parameter_ref<0 or unsigned(parameter_ref)>=the_ints.size())return;
	the_ints[parameter_ref].change_x(xx);
	
}

void Param_vecs::changeParameter(double xx,int parameter_ref)
{
	
  if(parameter_ref<0 or unsigned(parameter_ref)>=the_reals.size())return;
	the_reals[parameter_ref].change_x(xx);
	
}

void Param_vecs::changeParameter(std::string xx,int parameter_ref)
{
	
  if(parameter_ref<0 or unsigned(parameter_ref)>=the_strings.size())return;
	the_strings[parameter_ref].change_x(xx);
	
}

void Param_vecs::print_all_params(std::ostream *output)
{
	(*output) << std::endl;
	for(unsigned i=0;i<the_ints.size();i++)
	{
		(*output) << " # (" << i << ") # ";
		the_ints[i].print(output);
		(*output) << std::endl;
	}
	for(unsigned i=0;i<the_reals.size();i++)
	{
		(*output) << " # (" << i << ") # ";
		the_reals[i].print(output);
		(*output) << std::endl;
	}
	for(unsigned i=0;i<the_strings.size();i++)
	{
		(*output) << " # (" << i << ") # ";
		the_strings[i].print(output);
		(*output) << std::endl;
	}
	
}

int Param_vecs::returnInteger(int parameter_ref)
{
  if(parameter_ref<0 or unsigned(parameter_ref)>=the_ints.size())return 0;
	return the_ints[parameter_ref].return_value();
}
double Param_vecs::returnDouble(int parameter_ref)
{
  if(parameter_ref<0 or unsigned(parameter_ref)>=the_reals.size())return 0.;
	return the_reals[parameter_ref].return_value();

}

std::string Param_vecs::returnString(int parameter_ref)
{
  if(parameter_ref<0 or unsigned(parameter_ref)>=the_strings.size())return std::string(" ") ;
  return the_strings[parameter_ref].return_value();
  
}


void Param_vecs2::print_all_params(std::ostream *output)
{
  (*output) << " # Parameter Configuration File #\n";
  (*output) << " # <ints> " << the_ints.size();
  (*output) << " <doubles> " << the_reals.size();
  (*output) << " <strings> " << the_strings.size();
  (*output) << std::endl;
  for(unsigned i=0;i<the_ints.size();i++)
  {
	(*output) << " # (" << i << ") # ";
	the_ints[i].print_name(output);
	(*output) << " <int> ";
	the_ints[i].print_data(output);
	(*output) << std::endl;
  }
  for(unsigned i=0;i<the_reals.size();i++)
  {
	(*output) << " # (" << i << ") # ";
	the_reals[i].print_name(output);
	(*output) << " <double> ";
	the_reals[i].print_data(output);
	(*output) << std::endl;
  }
  for(unsigned i=0;i<the_strings.size();i++)
  {
	(*output) << " # (" << i << ") # ";
	the_strings[i].print_name(output);
	(*output) << " <string> ";
	the_strings[i].print_data(output);
	(*output) << std::endl;
  }
  (*output) << "\n# End of File # \n";
}

// read all parameters to stream
bool Param_vecs2::read_all_params(std::string filename)
{
  
  std::FILE *pfile;
  pfile = fopen (filename.c_str(),"r");
  
  int size_int,size_real,size_string,i_int,i_real,i_string,error;
  size_int = the_ints.size();size_real = the_reals.size();size_string = the_strings.size();
  i_int=0;i_real=0;i_string=0;error=0;
  
  if (pfile == NULL) std::cerr << "Error opening file " << filename << std::endl;
  else
  {
	while(! std::feof (pfile))
	{
	  int temp_int;
	  float temp_real=0.;
	  std::string temp_string;
	  char char_string[100];
	  
	  error = fscanf (pfile, "%s", char_string);
	  //std::cout << char_string << "\n";
	  if(!strcmp(char_string, "<ints>"))
	  {
		error = fscanf (pfile, "%i", &temp_int);
		if(temp_int!=size_int)
		{
		  error=0;
		  #ifdef DEBUG
		  std::cout << " Incompatible data file...\n";
		  #endif
		  break;
		}
		#ifdef DEBUG
		std::cout << " no of ints : " << temp_int;
		#endif
	  }
	  if(!strcmp(char_string, "<doubles>"))
	  {
		error = fscanf (pfile, "%i", &temp_int);
		#ifdef DEBUG
		std::cout << " no of doubles : " << temp_int;
		#endif
		if(temp_int!=size_real)
		{
		  	error=0;
			#ifdef DEBUG
			std::cout << "  Incompatible data file...\n";
			#endif
		  break;
		}
	  }
	  if(!strcmp(char_string, "<strings>"))
	  {
		error = fscanf (pfile, "%i", &temp_int);
		#ifdef DEBUG
		std::cout << " no of strings : " << temp_int;
		#endif
		if(temp_int!=size_string)
		{
		  error=0;
		  #ifdef DEBUG
		  std::cout << " Incompatible data file ...\n";
		  #endif
		  break;
		}
	  }
	  
	  if(!strcmp(char_string, "<int>"))
	  {
		error = fscanf (pfile, "%i", &temp_int);
		#ifdef DEBUG
		std::cout << " int parameter : " << temp_int;
		#endif
		changeParameter(temp_int,i_int);
		i_int++;
	  }
	  if(!strcmp(char_string, "<double>"))
	  {
		error = fscanf (pfile, "%f", &temp_real);
		#ifdef DEBUG
		std::cout << " real parameter : " << temp_real;
		#endif
		changeParameter(temp_real,i_real);
		i_real++;
	  }
	  if(!strcmp(char_string, "<string>"))
	  {
		error = fscanf (pfile, "%s", char_string);
		#ifdef DEBUG
		std::cout << " string parameter : " << char_string;
		#endif
		changeParameter(std::string(char_string),i_string);
		i_string++;
	  }
	}
	#ifdef DEBUG
	std::cout << " \n test :: " << error << " " << i_int << " " <<size_int << " " <<i_real << " " <<size_real  << " " << i_string << " " <<size_string;
	#endif
	fclose (pfile);
  }
  if(error and i_int==size_int and i_real==size_real and i_string==size_string)
  {
	#ifdef DEBUG
	std::cout << " \n Successful read of data...\n";
	#endif
	return true;
  }
  else return false;
  
}

