#include <cmath>

using namespace std;

double fzero_bisect(double (*func)(double), double bound_lower, double bound_upper, double err_tol, int iter_max);
double findRoot_secant(double (*function)(double), double guessInitial1, double guessInitial2, double errorTolerance, int iterMax);
double findRoot_regulafalsi(double (*function)(double), double boundLower, double boundUpper, double errorTolerance, int iterMax);
double findRoot_NewtonRaphson(double (*function)(double), double (*functionDerivative)(double), double guessInitial, double errorTolerance, int iterMax);
double findRoot_Halley(double (*function)(double), double (*functionDerivative1)(double), double (*functionDerivative2)(double), double guessInitial, double errorTolerance, int iterMax);