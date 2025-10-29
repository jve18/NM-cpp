#include <iostream>
#include <cmath>
#include "../../../source/NM-cpp.h"

double y(double x)
{
  return pow((x+2.0),2.0) - 3;
}

int main()
{

  double x_1 = -1.0; // First guess
  double x_2 = 0.0; // Second guess
  double eps_y = pow(10.0,-4.0); // Error tolerance
  int i_max = 10000; // Max number of iterations
  double x_root;

  x_root = fzero_secant(y, x_1, x_2, eps_y, i_max);

  std::cout << x_root << std::endl;

  return 0;
  
  

}
