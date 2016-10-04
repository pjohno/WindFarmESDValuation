#include "DelayedWindFarm.hpp"
#include <iostream>
#include <string>
using namespace std;

namespace DelayedWindFarm
{
  int testProject(void)
  {
    string home_location=getenv("HOME");
    int exit_status = 0;
    cout << " My new project:: DelayedWindFarm " << endl;
    cout << " Home location " << home_location << endl;
    return exit_status;
  }
}

