#include "vary_params.h"
#include <iostream>

int main(int argc, const char* argv[])
{
  
  if(argc>=2)
  {
    Vary_Params test {std::string(argv[1])};
    test.runIterator();
  }
  else
  {
    Vary_Params test;
    test.runIterator();
  }
  
}

