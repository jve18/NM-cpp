#include <iostream>
#include <cmath>
#include "../../source/NM-cpp.h"
using namespace std;

double y(double x)
{
  return x - 2.3;
}

int main()
{

  double x_L = 0; // Lower bound guess
  double x_U = 4; // Upper boung guess
  double eps_y = pow(10.0,-4.0); // Error tolerance
  int i_max = pow(10.0,6.0); // Max number of iterations
  double x_root;

  x_root = fzero_bisect(y,x_L,x_U,eps_y,i_max);

  cout << x_root << endl;

  return 0;
  
  

}
