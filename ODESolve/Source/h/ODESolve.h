#include <cmath>
#include <vector>

using namespace std;

vector <double> ODESolve_RK4(double (*dy)(double, double), double x0, double y0, double xmax, double dx);