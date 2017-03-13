
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <string.h>
#include "pj_data.h"
#include "param_maps.h"

//#define DEBUG
using namespace PJLib;  

void ParameterMaps::clear(void)
{
	the_ints.clear();the_reals.clear();the_strings.clear();
}

void ParameterMaps::reset()
{
  for(ii=the_ints.begin();ii!=the_ints.end();ii++)
  {
	(*ii).second.reset();
  }
  for(jj=the_reals.begin();jj!=the_reals.end();jj++)
  {
	(*jj).second.reset();
  }
  for(kk=the_strings.begin();kk!=the_strings.end();kk++)
  {
	(*kk).second.reset();
  }
}

bool ParameterMaps::eraseParameter(std::string ref)
{
  if(the_ints.erase(ref))return true;
  if(the_reals.erase(ref))return true;
  if(the_strings.find(ref)==the_strings.end())
  {
    #ifdef DEBUG
    std::cout << " could not find a parameter for < " << ref << " > \n";
    #endif
    return false;
  }
  else {the_strings.erase(ref);return true;}
}

bool ParameterMaps::check(const std::string ref)
{
  if(the_ints.find(ref)!=the_ints.end())return true;
  if(the_reals.find(ref)!=the_reals.end())return true;
  if(the_strings.find(ref)!=the_strings.end())return true;
  return false;
}

void ParameterMaps::assignParameter(int *xx,std::string ref,std::string comment)
{
	PJ_data2<int> integer_parameter(comment,xx);
	the_ints[ref] = integer_parameter;
}

void ParameterMaps::assignParameter(int *xx,std::string ref,std::string comment,int reset_value)
{
	PJ_data2<int> integer_parameter(comment,xx,reset_value);
	the_ints[ref] = integer_parameter;
}

void ParameterMaps::assignParameter(int *xx,std::string ref,std::string comment,int reset_value,int min,int max)
{
	PJ_data2<int> integer_parameter(comment,xx,reset_value,min,max);
	the_ints[ref] = integer_parameter;
}

void ParameterMaps::assignParameter(double *xx,std::string ref,std::string comment)
{
	PJ_data2<double> double_parameter(comment,xx);
	the_reals[ref] = double_parameter;
}

void ParameterMaps::assignParameter(double *xx,std::string ref,std::string comment,double reset_value)
{
	PJ_data2<double> double_parameter(comment,xx,reset_value);
	the_reals[ref] = double_parameter;
}

void ParameterMaps::assignParameter(double *xx,std::string ref,std::string comment,double reset_value,double min,double max)
{
	PJ_data2<double> double_parameter(comment,xx,reset_value,min,max);
	the_reals[ref] = double_parameter;
}

void  ParameterMaps::assignParameter(std::string *xx,std::string ref,std::string comment)
{
	PJ_data2<std::string> string_parameter(comment,xx);
	the_strings[ref] = string_parameter;
}

void  ParameterMaps::assignParameter(std::string *xx,std::string ref,std::string comment,std::string reset_value)
{
	PJ_data2<std::string> string_parameter(comment,xx,reset_value);
	the_strings[ref] = string_parameter;
}

bool ParameterMaps::changeParameter(int xx,std::string ref)
{
	
  if(the_ints.find(ref)==the_ints.end())
  {
	#ifdef DEBUG
	std::cout << " could not find parameter < " << ref << " > \n";
	#endif
	return false;
  }
  the_ints[ref].change_x(xx);return true;
	
}

bool ParameterMaps::changeParameter(double xx,std::string ref)
{
	
  if(the_reals.find(ref)==the_reals.end())
  {
	#ifdef DEBUG
	std::cout << " could not find parameter < " << ref << " > \n";
	#endif
	return false;
  }
  the_reals[ref].change_x(xx);return true;
	
}

bool ParameterMaps::changeParameter(std::string xx,std::string ref)
{
	
  if(the_strings.find(ref)==the_strings.end())
  {
	#ifdef DEBUG
	std::cout << " could not find parameter < " << ref << " > \n";
	#endif
	return false;
  }
  the_strings[ref].change_x(xx);return true;
	
}

int ParameterMaps::returnInteger(std::string ref)
{
  if(the_ints.find(ref)==the_ints.end())return 0;
  return the_ints[ref].return_value();
}
double ParameterMaps::returnDouble(std::string ref)
{
  if(the_reals.find(ref)==the_reals.end())return 0.;
  return the_reals[ref].return_value();

}

std::string ParameterMaps::returnString(std::string ref)
{
  if(the_strings.find(ref)==the_strings.end())return std::string(" ") ;
  return the_strings[ref].return_value();
}

int* ParameterMaps::returnPtrInteger(std::string ref)
{
  if(the_ints.find(ref)==the_ints.end())return NULL;
  return the_ints[ref].return_pointer();
}
double* ParameterMaps::returnPtrDouble(std::string ref)
{
  if(the_reals.find(ref)==the_reals.end())return NULL;
  return the_reals[ref].return_pointer();
}

std::string* ParameterMaps::returnPtrString(std::string ref)
{
  if(the_strings.find(ref)==the_strings.end())return NULL ;
  return the_strings[ref].return_pointer();
}

void ParameterMaps::print(std::string ref,std::ostream *output)
{
  ii = the_ints.find(ref);
  if(ii!=the_ints.end())(*ii).second.print(output);
  jj=the_reals.find(ref);
  if(jj!=the_reals.end())(*jj).second.print(output);
  kk = the_strings.find(ref);
  if(kk!=the_strings.end())(*kk).second.print(output);
}

std::string ParameterMaps::return_name(std::string ref)
{
  ii = the_ints.find(ref);
  if(ii!=the_ints.end())return (*ii).second.return_name();
  jj=the_reals.find(ref);
  if(jj!=the_reals.end())return (*jj).second.return_name();
  kk = the_strings.find(ref);
  if(kk!=the_strings.end())return (*kk).second.return_name();
  return " ";
}

bool ParameterMaps::is_default(std::string ref)
{
  ii = the_ints.find(ref);
  if(ii!=the_ints.end())return (*ii).second.is_default();
  jj=the_reals.find(ref);
  if(jj!=the_reals.end())return (*jj).second.is_default();
  kk = the_strings.find(ref);
  if(kk!=the_strings.end())return (*kk).second.is_default();
  return false;
}

bool ParameterMaps::is_bounded(std::string ref)
{
  ii = the_ints.find(ref);
  if(ii!=the_ints.end())return (*ii).second.is_bounded();
  jj=the_reals.find(ref);
  if(jj!=the_reals.end())return (*jj).second.is_bounded();
  kk = the_strings.find(ref);
  if(kk!=the_strings.end())return (*kk).second.is_bounded();
  return false;
}

void ParameterMaps::print_all_params(std::ostream *output)
{
  (*output) << " # Parameter Configuration File #\n#\n";
  (*output) << " # <ints> " << the_ints.size();
  (*output) << " <doubles> " << the_reals.size();
  (*output) << " <strings> " << the_strings.size();
  (*output) << std::endl;
  for(ii=the_ints.begin();ii!=the_ints.end();ii++)
  {
	(*output) << " #  ";
	(*ii).second.print_name(output);
	(*output) << " <int>< ";
	(*output) << (*ii).first << " > ";
	(*ii).second.print_data(output);
	(*output) << std::endl;
  }
  for(jj=the_reals.begin();jj!=the_reals.end();jj++)
  {
	(*output) << " #  ";
	(*jj).second.print_name(output);
	(*output) << " <double>< ";
	(*output) << (*jj).first << " > ";
	(*jj).second.print_data(output);
	(*output) << std::endl;
  }
  for(kk=the_strings.begin();kk!=the_strings.end();kk++)
  {
	(*output) << " #  ";
	(*kk).second.print_name(output);
	(*output) << " <string>< ";
	(*output) << (*kk).first << " > ";
	(*kk).second.print_data(output);
	(*output) << std::endl;
  }
  (*output) << "#\n # End of File # \n";
}

// read all parameters to stream
bool ParameterMaps::read_all_params(std::string filename)
{
  std::FILE *pfile;
  pfile = fopen (filename.c_str(),"r");
  
  std::string::iterator string_it;

  int size_int,size_real,size_string,i_int,i_real,i_string,error;
  size_int=0;size_real=0;size_string=0;
  i_int=0;i_real=0;i_string=0;error=0;
  
  if (pfile == NULL) perror ("Error opening file");
  else
  {
	#ifdef DEBUG
	std::cout << " Read parameters from file :: " << filename << "\n";
	#endif
	while(! std::feof (pfile))
	{
	  int temp_int;
	  float temp_real=0.;
	  std::string temp_string;
	  char char_string[100];
	  char name_string[100];
	  
	  error = fscanf (pfile, "%s", char_string);
	  //std::cout << char_string << "\n";
	  if(!strcmp(char_string, "<ints>"))
	  {
		if(!fscanf (pfile, "%i", &size_int))
		{
		  #ifdef DEBUG
		  std::cout << " no of ints : " << size_int;
		  #endif
		  error=0;
		  #ifdef DEBUG
		  std::cout << " Incompatible data file...\n";
		  #endif
		  break;
		}
	  } // number of ints
	  if(!strcmp(char_string, "<doubles>"))
	  {
		if(!fscanf (pfile, "%i", &size_real))
		{
		  #ifdef DEBUG
		  std::cout << " no of doubles : " << size_real;
		  #endif
		  error=0;
		  #ifdef DEBUG
		  std::cout << "  Incompatible data file...\n";
		  #endif
		  break;
		}
	  } // number of doubles
	  if(!strcmp(char_string, "<strings>"))
	  {
		if(!fscanf (pfile, "%i", &size_string))
		{
		  #ifdef DEBUG
		  std::cout << " no of strings : " << size_string;
		  #endif
		  error=0;
		  #ifdef DEBUG
		  std::cout << " Incompatible data file ...\n";
		  #endif
		  break;
		}
	  } // number of strings
	  
	  if(!strcmp(char_string, "<int><"))
	  {
		error = fscanf (pfile, "%s", name_string);
		error = fscanf (pfile, "%s", char_string);
		if(!fscanf (pfile, "%i", &temp_int))i_int--;
		#ifdef DEBUG
		std::cout << " << " << name_string << " >> ";
		std::cout << " int parameter : " << temp_int;
		#endif
		if(changeParameter(temp_int,name_string))i_int++;
	  } // int
	  if(!strcmp(char_string, "<double><"))
	  {
		error = fscanf (pfile, "%s", name_string);
		error = fscanf (pfile, "%s", char_string);
		if(!fscanf (pfile, "%f", &temp_real))i_real--;
		#ifdef DEBUG
		std::cout << " << " << name_string << " >> ";
		std::cout << " real parameter : " << temp_real;
		#endif
		if(changeParameter(temp_real,name_string))i_real++;
	  } // double
	  if(!strcmp(char_string, "<string><"))
	  {
		error = fscanf (pfile, "%s", name_string);
		error = fscanf (pfile, "%s", char_string);
		error = fscanf (pfile, "%s", char_string);
		#ifdef DEBUG
		std::cout << " << " << name_string << " >> ";
		std::cout << " string parameter : " << char_string;
		#endif
		if(changeParameter(std::string(char_string),std::string(name_string)))i_string++;
	  } // string
	} // end of file
	#ifdef DEBUG
	std::cout << " \n test :: " << error << " " << i_int << " " <<size_int << " " <<i_real << " " <<size_real  << " " << i_string << " " <<size_string;
	#endif
	fclose (pfile);
  } // file open
  if(error and i_int==size_int and i_real==size_real and i_string==size_string)
  {
	#ifdef DEBUG
	std::cout << " \n Successful read of data...\n";
	#endif
	return true;
  }
  else 
  {
	std::cout << " file not read...\n";
	return false;
  }
}

bool ParameterMaps::create_vector_refs(std::vector<std::string> &my_ints,std::vector<std::string> &my_doubles,std::vector<std::string> &my_strings)
{
  my_ints.clear();my_doubles.clear();my_doubles.clear();
  for(ii=the_ints.begin();ii!=the_ints.end();ii++)
  {
	my_ints.push_back((*ii).first);
  }
  for(jj=the_reals.begin();jj!=the_reals.end();jj++)
  {
	my_doubles.push_back((*jj).first);
  }
  for(kk=the_strings.begin();kk!=the_strings.end();kk++)
  {
	my_strings.push_back((*kk).first);
  }
}

