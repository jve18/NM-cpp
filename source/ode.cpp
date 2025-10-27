#include "ode.h"
#include <iostream>
#include <complex>
#include <vector>

using namespace std;

vector<double> ode_RK4(double (*dy)(double, double), double x0, double y0, double xmax, double dx)
{
    double y1 = y0;
    double x1 = x0;
    double k1 = dx*dy(x1, y1);
    double k2 = dx*dy(x1 + dx/2, y1 + k1/2);
    double k3 = dx*dy(x1 + dx/2, y1 + k2/2);
    double k4 = dx*dy(x1 + dx, y1 + k3);
    
    vector<double> vect_x1;
    vector<double> vect_y1;
    
    int nIter = 0;
    
    while(x1 < xmax){
            k1 = dx*dy(x1, y1);
            k2 = dx*dy(x1 + dx/2, y1 + k1/2);
            k3 = dx*dy(x1 + dx/2, y1 + k2/2);
            k4 = dx*dy(x1 + dx, y1 + k3);
            
            y1 = y1 + (1.0/6)*(k1 + 2*k2 + 2*k3 + k4);
            x1 = x1 + dx;
            
            vect_x1.push_back(x1);
            vect_y1.push_back(y1);

            nIter++;   
    }
    return vect_y1;
    
}