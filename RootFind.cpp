
#include "RootFind.h"

using namespace std;

double findRoot_bisection(double (*function)(double), double boundLower, double boundUpper, double errorTolerance){
    double x0 = boundLower;
    double f0 = function(x0);
    double x2 = boundUpper;
    double f2 = function(x2);
    double x1 = (x0 + x2)/2;
    double f1 = function(x1);
    
    double errorCurrent = 10;
    
    while(abs(errorCurrent) > errorTolerance)
    {
        f0 = function(x0);
        
        f2 = function(x2);
        
        x1 = (x0 + x2)/2;
        f1 = function(x1);
        
        if(f0 * f1 > 0){
            x0 = x1;
        } else{
            x2 = x1;
        }
        
        errorCurrent = f1;
    }
    
    return x1;
}

double findRoot_secant(double (*function)(double), double guessInitial1, double guessInitial2, double errorTolerance, int iterMax){
    double x0 = guessInitial1;
    double f0 = function(x0);
    double x1 = guessInitial2;
    double f1 = function(x1);
    double x2 = (x0*f1 - x1*f0)/(f1 - f0);
    double f2 = function(x2);
    
    double errorCurrent = 10;
    
    int iterCurrent = 0;
    
    while(abs(errorCurrent) > errorTolerance && iterCurrent <= iterMax){
        f0 = function(x0);
        f1 = function(x1);
        x2 = (x0*f1 - x1*f0)/(f1 - f0);
        f2 = function(x2);
        x0 = x1;
        x1 = x2;
        errorCurrent = f2;
        iterCurrent++;
        
    }
    return x2;
}

double findRoot_regulafalsi(double (*function)(double), double boundLower, double boundUpper, double errorTolerance){
    double x0 = boundLower;
    double f0 = function(x0);
    double x2 = boundUpper;
    double f2 = function(x2);
    double x1 = x0 - ((f0*(x2-x0))/(f2-f0));
    double f1 = function(x1);
    
    double errorCurrent = 10;
    
    while(abs(errorCurrent) > errorTolerance)
    {
        f0 = function(x0);
        f2 = function(x2);
        x1 = x0 - ((f0*(x2-x0))/(f2-f0));
        f1 = function(x1);
        
        if(f0 * f1 > 0){
            x0 = x1;
        } else{
            x2 = x1;
        }
        
        errorCurrent = f1;
    }
    
    return x1;   
}

double findRoot_NewtonRaphson(double (*function)(double), double (*functionDerivative)(double), double guessInitial, double errorTolerance, int iterMax){
    double x0 = guessInitial;
    double f0 = function(x0);
    double fp0 = functionDerivative(x0);
    
    double errorCurrent = f0;
    
    int iterCurrent = 0;
    
    while(abs(errorCurrent) > errorTolerance && iterCurrent <= iterMax){
        f0 = function(x0);
        fp0 = functionDerivative(x0);
        x0 = x0 - f0/fp0;
        errorCurrent = f0;
        iterCurrent++;
    }
    return x0;
}

double findRoot_Halley(double (*function)(double), double (*functionDerivative1)(double), double (*functionDerivative2)(double), double guessInitial, double errorTolerance, int iterMax){
    double x0 = guessInitial;
    double f0 = function(x0);
    double fp0 = functionDerivative1(x0);
    double fpp0 = functionDerivative2(x0);
    
    double errorCurrent = f0;
    
    int iterCurrent = 0;
    
    while(abs(errorCurrent) > errorTolerance && iterCurrent <= iterMax){
        f0 = function(x0);
        fp0 = functionDerivative1(x0);
        fpp0 = functionDerivative2(x0);
        x0 = x0 - (2*f0*fp0)/(2*fp0*fp0-f0*fpp0);
        errorCurrent = f0;
        iterCurrent++;
    }
    return x0;
}