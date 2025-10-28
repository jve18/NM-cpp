#include <cmath>
#include <vector>

std::vector <double> ode_RK4(double (*dy)(double, double), double x0, double y0, double xmax, double dx);