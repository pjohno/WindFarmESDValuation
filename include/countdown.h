#ifndef _countdown_h_included_
#define _countdown_h_included_

#include <iostream>
#include <iomanip>
#include <sstream>

class Countdown{

  char *buffer;
  char *percent;
  int ii,jj,i_step;
	
  public:
	
	~Countdown(){
	  delete [] buffer;
	  delete [] percent;
	}
	Countdown(){
	  buffer = new char[45];
	  percent = new char[10];
	  sprintf(percent, "%3.2f%%", 0.0);
	  buffer[0] = '[';i_step = 0;buffer[42]='\0';
	}
	void update(int curr,int finish)
	{
	  if(curr>=finish){
		for(int k = 1; k < 41; k++)buffer[k] = '-';
		buffer[41] = ']';
		sprintf(percent, "%3.2f%%", 100.0);
		printf("%s%s", buffer, percent);
		return;
	  }
	  ii = int(double(curr)/double(finish) * 40.);jj = curr % 4;
	  if(jj == 0)buffer[ii + 1] = '\\';
	  else if(jj == 1)buffer[ii + 1] = '|';
	  else if(jj == 2)buffer[ii + 1] = '/';
	  else buffer[ii + 1] = '-';
	  if(i_step<ii){for(int k=1;k<ii+1;k++)buffer[k] = '-';i_step = ii;}
	  for(int k = ii + 2; k < 41; k++)buffer[k] = ' ';buffer[41] = ']';
	  std::sprintf(percent, "%3.2f%%", double(curr) / double(finish) * 100.0);
	  std::printf("%s%s\r", buffer, percent);
	  std::cout.flush();
	}

};

#endif

