#include <iostream>
#include <cmath>
#include "../../../source/NM-cpp.h"

double y(double x)
{
  return pow((x+2.0),2.0) - 3.0;
}

double y_p(double x)
{
  return 2.0*(x + 2.0);
}

double y_pp(double x)
{
  return 2.0;
}



int main()
{

  double x_0 = 0.0; // Initial guess
  double eps_y = pow(10.0,-4.0); // Error tolerance
  int i_max = 10000; // Max number of iterations
  double x_root;

  x_root = fzero_Halley(y, y_p, y_pp, x_0, eps_y, i_max);

  std::cout << x_root << std::endl;

  return 0;
  
  

}
