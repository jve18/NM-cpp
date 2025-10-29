#include <iostream>
#include <cmath>
#include "../../../source/NM-cpp.h"

double y(double x)
{
  return pow((x+2.0),2.0) - 3;
}

int main()
{

  double x_L = -2.0; // Lower bound guess
  double x_U = 0.0; // Upper bound guess
  double eps_y = pow(10.0,-4.0); // Error tolerance
  int i_max = 10000; // Max number of iterations
  double x_root;

  x_root = fzero_bisect(y, x_L, x_U, eps_y, i_max);

  std::cout << x_root << std::endl;

  return 0;
  
  

}
