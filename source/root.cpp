
#include "root.h"
#include <iostream>
#include <complex>

// FUNCTIONS

double fzero_bisect(double (*func)(double), double bound_lower, double bound_upper, double err_tol, int iter_max){
    /*

    DESCRIPTION: Function finds the root/zero of a function using the root bisection method.
     
    INPUTS:
    *  func: The function that the user would like to find the root/zero of.
    *  bound_lower: The lower bound of the domain the root/zero is known to exist within.
    *  bound_upper: The upper bound of the domain the root/zero is known to exist within.
    *  err_tol: The amount by which the function evaluated at the best estimate of the root/zero is allowed to differ from 0 by; f(x1); how far from 0 can f(x1) be.
    *  iter_max: Maximum number of iterations for the bisection algorithm before it stops.
    
    OUTPUTS:
    *  x1: The approximated root of the function
    
    SOURCES:
    *  https://en.wikipedia.org/wiki/Bisection_method
    *  TAMU AERO 220 Summer 2016 Lecture Notes
    *  Numerical Methods in Engineering with MATLAB by Jaan Kiusalaas (3rd ed.)
    */

    // Error handling for entering a lower bound and upper bound that are the same
    try{
        if(bound_lower == bound_upper){
            throw 10001;
        }
        
    }
    catch(int err10001){
        std::cout << "Error: " << "Input lower bound equal to input upper bound. Invalid input. Disregard return value." << endl;
        return 0;
    }
    
    // Error handling for entering a lower bound greater than the upper bound.
    try{
        if(bound_lower > bound_upper){
            throw 10002;
        }
        
    }
    catch(int err10002){
        std::cout << "Warning: " << "Input lower bound greater than input upper bound. Evaluate reasonability of return value." << endl;
        
    }
    
    double x0 = bound_lower;                     //x0, current lower bound, initialized to user-defined lower bound
    double f0 = func(x0);                    //f0, function evaluated at current lower bound
    double x2 = bound_upper;                     //x2, current upper bound, initialized to user-defined upper bound
    double f2 = func(x2);                    //f2, function evaluated at current lower bound
    double x1 = (x0 + x2)/2;                     //x1, midpoint based on current lower and upper bounds
    double f1 = func(x1);                    //f1, function evaluated at current midpoint
    
    double err = f1;
    int iter = 0;
    
    while(abs(err) > err_tol && iter <= iter_max)
    {
        f0 = func(x0);                      //re-evaluate f0 based on new x0
        f2 = func(x2);                      //re-evaluate f2 based on new x2
        x1 = (x0 + x2)/2;                       //re-calculate midpoint based on new x0, x2
        f1 = func(x1);                      //re-evaluate f1 based on new midpoint, x1
        
        if(f0 * f1 > 0){                        //Is the sign of f0 and f1 the same?
            x0 = x1;                            //If so, the lower bound can be moved up to x1, so set x0 = x1
        } else{
            x2 = x1;                            //Otherwise, f1 and f2 must be of the same size, the upper bound can be moved down to x1, set x1 = x2
        }
        
        err = f1;                               //set error as function evaluated at current estimate of root
        iter++;                                 //increment iteration
    }
    
    // Error handling for if max number of iterations  
    try{
        if(iter >= iter_max){
            throw 10000;
        }
        
    }
    catch(int err10000){
        std::cout << "Error: " << "Maximum iteration limit exceeded" << endl;
        std::cout << "Disregard return value." << endl;
        std::cout << "Possible sources of error: " << endl;
        std::cout << "Root is not bracketed. Root must be bracketed for bisection." << endl;
        std::cout << "Not enough iterations were alloted. Bisection is one of the slower root finding techniques. A high number of iterations may be necessary ~ log2((boundUpper-boundLower)/errorTolerance)." << endl;
        return 0;
    }
   
    std::cout << "iter: " << iter << endl;
    return x1;
    
}

double findRoot_secant(double (*function)(double), double guessInitial1, double guessInitial2, double errorTolerance, int iterMax){
     /*
     INPUTS:
     *  function: The function that the user would like to find the root (0) of
     *  guessInitial1: First initial guess near the root. DOES NOT have to bracket the root. DOES NOT have to be less than guessInitial2
     *  guessInitial1: Second initial guess near the root. DOES NOT have to bracket the root. DOES NOT have to be greater than guessInitial1
     *  errorTolerance: The amount by which the function evaluated at the best estimate of the root is allowed to differ from 0 by; f(x1); how far from 0 can f(x1) be
     *  iterMax: Maximum number of iterations for the secant algorithm before it stops
     
     OUTPUTS:
     *  x2: The approximated root of the function
     
     SOURCES:
     *  https://en.wikipedia.org/wiki/Secant_method
     *  TAMU AERO 220 Summer 2016 Lecture Notes
     *  Numerical Methods in Engineering with MATLAB by Jaan Kiusalaas (3rd ed.)
     */
    
     // Error handling for entering a lower bound and upper bound that are the same
    try{
        if(guessInitial1 == guessInitial2){
            throw 10001;
        }
        
    }
    catch(int err10001){
        std::cout << "Error: " << "Input first initial guess equal to second initial guess. Invalid input. Disregard return value." << endl;
        return 0;
    }
    
    
    double x0 = guessInitial1;                          //x0, current first guess, initialized to user-defined first guess
    double f0 = function(x0);                           //f0, function evaluated at current first guess
    double x1 = guessInitial2;                          //x1, current second guess, initialized to user-defined first guess
    double f1 = function(x1);                           //f1, function evaluated at current second guess
    double x2 = (x0*f1 - x1*f0)/(f1 - f0);              //x1, current third guess (will replace second guess, x2 -> x1, after the second guess replaces the first guess, x1 -> x0
    double f2 = function(x2);                           //x2, function evaluated at third guess
    
    double errorCurrent = 10;
    
    int iterCurrent = 0;
    
    while(abs(errorCurrent) > errorTolerance && iterCurrent <= iterMax){
        f0 = function(x0);                              //re-evaluate f0 based on new x0
        f1 = function(x1);                              //re-evaluate f1 based on new x1
        x2 = (x0*f1 - x1*f0)/(f1 - f0);                 //re-calculate x2 based on new x0, f0, x1, f1
        f2 = function(x2);                              //re-evaluate f2 based on new x2
        x0 = x1;                                        //set next "first" guess equal to current "second" guess
        x1 = x2;                                        //set next "second" guess equal to current "third" guess
        errorCurrent = f2;
        iterCurrent++;
        
    }
    
    // Error handling for if max number of iterations  
    try{
        if(iterCurrent >= iterMax){
            throw 10000;
        }
        
    }
    catch(int err10000){
        std::cout << "Error: " << "Maximum iteration limit exceeded." << endl;
        std::cout << "Disregard return value." << endl;                      
        std::cout << "Possible sources of error: " << endl;
        std::cout << "Trying converge function extrema; secant method may not converge on function extrema." << endl;
        std::cout << "Inappropriate initial guesses; some initial guesses can cause secant method to DIVERGE." << endl;
        return 0;
    }
    
    std::cout << "iter: " << iterCurrent << endl;
    return x2;
    
}

double findRoot_regulafalsi(double (*function)(double), double boundLower, double boundUpper, double errorTolerance, int iterMax){
    /*
     INPUTS:
     *  function: The function that the user would like to find the root (0) of
     *  boundLower: The lower bound of the domain the root is known to exist within
     *  boundUpper: The upper bound of the domain the root is known to exist within
     *  errorTolerance: The amount by which the function evaluated at the best estimate of the root is allowed to differ from 0 by; f(x1); how far from 0 can f(x1) be
     *  iterMax: Maximum number of iterations for the regula falsi algorithm before it stops
     
     OUTPUTS:
     *  x1: The approximated root of the function
     
     SOURCES:
     *  https://en.wikipedia.org/wiki/False_position_method
     *  TAMU AERO 220 Summer 2016 Lecture Notes
     *  Numerical Methods in Engineering with MATLAB by Jaan Kiusalaas (3rd ed.)
     */
    
    // Error handling for entering a lower bound and upper bound that are the same
    try{
        if(boundLower == boundUpper){
            throw 10001;
        }
        
    }
    catch(int err10001){
        std::cout << "Error: " << "Input lower bound equal to input upper bound. Invalid input. Disregard return value." << endl;
        return 0;
    }
    
    // Error handling for entering a lower bound greater than the upper bound.
    try{
        if(boundLower > boundUpper){
            throw 10002;
        }
        
    }
    catch(int err10002){
        std::cout << "Warning: " << "Input lower bound greater than input upper bound. Evaluate reasonability of return value." << endl;
        
    }
    
    double x0 = boundLower;                                 //x0, current lower bound, initialized to user-defined lower bound
    double f0 = function(x0);                               //f0, function evaluated at current lower bound
    double x2 = boundUpper;                                 //x2, current upper bound, initialized to user-defined upper bound
    double f2 = function(x2);                               //f2, function evaluated at current lower bound
    double x1 = x2 - ((f2*(x2-x0))/(f2-f0));                //x1, false position point calculated based on x0 and x2
    double f1 = function(x1);                               //f1, function evaluated at false position point
    
    double errorCurrent = 10;
    int iterCurrent = 0;
    
    while(abs(errorCurrent) > errorTolerance && iterCurrent <= iterMax)
    {
        f0 = function(x0);                                  //re-evaluate f0 based on new x0
        f2 = function(x2);                                  //re-evaluate f2 based on new x2
        x1 = x2 - ((f2*(x2-x0))/(f2-f0));                   //re-calculate x1 based on new x0, f0, x2, f2
        f1 = function(x1);                                  //re-evaluate f1 based on new x1
        
        if(f0 * f1 > 0){                                    //Is the sign of f0 and f1 the same?
            x0 = x1;                                        //If so, the lower bound can be moved up to x1, so set x0 = x1
        } else{
            x2 = x1;                                        //Otherwise, f1 and f2 must be of the same size, the upper bound can be moved down to x1, set x1 = x2
        }
        
        errorCurrent = f1;
        iterCurrent++;
        
    }
    
    // Error handling for if max number of iterations  
    try{
        if(iterCurrent >= iterMax){
            throw 10000;
        }
        
    }
    catch(int err10000){
        std::cout << "Error: " << "Maximum iteration limit exceeded" << endl;
        std::cout << "Disregard return value." << endl;
        std::cout << "Possible sources of error: " << endl;
        std::cout << "Root is not bracketed. Root must be bracketed for regula-falsi." << endl;
        std::cout << "Not enough iterations were alloted." << endl;
        return 0;
    }
    
    std::cout << "iter: " << iterCurrent << endl;
    return x1;   
    
}

double findRoot_NewtonRaphson(double (*function)(double), double (*functionDerivative)(double), double guessInitial, double errorTolerance, int iterMax){
    /*
     INPUTS:
     *  function: The function that the user would like to find the root (0) of
     *  functionDerivative: The 1st derivative of the function that the user would like to find the root (0) of
     *  guessInitial: The initial guess for where the root might be
     *  errorTolerance: The amount by which the function evaluated at the best estimate of the root is allowed to differ from 0 by; f(x1); how far from 0 can f(x1) be
     *  iterMax: Maximum number of iterations for the Newton-Raphson algorithm before it stops
     
     OUTPUTS:
     *  x0: The approximated root of the function
     
     SOURCES:
     *  https://en.wikipedia.org/wiki/Newton%27s_method
     *  TAMU AERO 220 Summer 2016 Lecture Notes
     *  Numerical Methods in Engineering with MATLAB by Jaan Kiusalaas (3rd ed.)
     */
    
    double x0 = guessInitial;                                               //x0, current guess of the function root, initialized with user-inputed initial guess of root
    double f0 = function(x0);                                               //f0, function evaluated at current guess of root
    double fp0 = functionDerivative(x0);                                    //fp0, function 1st derivative evaluated at current guess of root
    
    // Error handling for entering an initial guess that causes fp0 to equal 0
    try{
        if(fp0 == 0){
            throw 10000;
        }
        
    }
    catch(int err10000){
        std::cout << "Error: " << "Function first derivative evaluated at initial guess equal to 0. Invalid input. Disregard return value." << endl;
        return 0;
    }
    
    double errorCurrent = f0;
    
    int iterCurrent = 0;
    
    while(abs(errorCurrent) > errorTolerance && iterCurrent <= iterMax){
        f0 = function(x0);                                                  //re-evaluate f0 based on new x0
        fp0 = functionDerivative(x0);                                       //re-evaluate fp0 based on new x0
        x0 = x0 - f0/fp0;                                                   //re-calculate x0 based on current x0, f0, fp0
        errorCurrent = f0;
        iterCurrent++;
    }
    
    // Error handling for if max number of iterations  
    try{
        if(iterCurrent >= iterMax){
            throw 10001;
        }
        
    }
    catch(int err10001){
        std::cout << "Error: " << "Maximum iteration limit exceeded" << endl;
        std::cout << "Disregard return value." << endl;
        std::cout << "Possible sources of error: " << endl;
        std::cout << "Inappropriate initial guess; Newton-Raphson algorithm can diverge." << endl;
        std::cout << "Function derivative may be incorrect." << endl;
        std::cout << "Function derivative may have been equal to 0 (division by 0) during an iteration." << endl;
        std::cout << "Not enough iterations were alloted." << endl;
        return 0;
    }
    
    std::cout << "iter: " << iterCurrent << endl;
    return x0;
    
}

double findRoot_Halley(double (*function)(double), double (*functionDerivative1)(double), double (*functionDerivative2)(double), double guessInitial, double errorTolerance, int iterMax){
     /*
     INPUTS:
     *  function: The function that the user would like to find the root (0) of
     *  functionDerivative1: The 1st derivative of the function that the user would like to find the root (0) of
     *  functionDerivative2: The 2nd derivative of the function that the user would like to find the root (0) of
     *  guessInitial: The initial guess for where the root might be
     *  errorTolerance: The amount by which the function evaluated at the best estimate of the root is allowed to differ from 0 by; f(x1); how far from 0 can f(x1) be
     *  iterMax: Maximum number of iterations for the Newton-Raphson algorithm before it stops
     
     OUTPUTS:
     *  x0: The approximated root of the function
     
     SOURCES:
     *  https://en.wikipedia.org/wiki/Newton%27s_method
     *  TAMU AERO 220 Summer 2016 Lecture Notes
     *  Numerical Methods in Engineering with MATLAB by Jaan Kiusalaas (3rd ed.)
     */
    
    double x0 = guessInitial;                                               //x0, current guess of the function root, initialized with user-inputed initial guess of root
    double f0 = function(x0);                                               //f0, function evaluated at current guess of root
    double fp0 = functionDerivative1(x0);                                   //fp0, function 1st derivative evaluated at current guess of root
    double fpp0 = functionDerivative2(x0);                                  //fpp0, function 2nd derivative evaluated at current guess of root
    
    double errorCurrent = f0;
    
    int iterCurrent = 0;
    
    while(abs(errorCurrent) > errorTolerance && iterCurrent <= iterMax){
        f0 = function(x0);                                                  //re-evaluate f0 based on new x0
        fp0 = functionDerivative1(x0);                                      //re-evaluate fp0 based on new x0
        fpp0 = functionDerivative2(x0);                                     //re-evaluate fpp0 based on new x0
        x0 = x0 - (2*f0*fp0)/(2*fp0*fp0-f0*fpp0);                           //re-calculate x0 based on current x0, f0, fp0
        errorCurrent = f0;
        iterCurrent++;
    }
    
    // Error handling for if max number of iterations  
    try{
        if(iterCurrent >= iterMax){
            throw 10001;
        }
        
    }
    catch(int err10001){
        std::cout << "Error: " << "Maximum iteration limit exceeded" << endl;
        std::cout << "Disregard return value." << endl;
        std::cout << "Possible sources of error: " << endl;
        std::cout << "Inappropriate initial guess; Halley algorithm can diverge." << endl;
        std::cout << "Function 1st or 2nd derivative may be incorrect." << endl;
        std::cout << "Denominator may have been equal to 0 (division by 0) during an iteration." << endl;
        std::cout << "Not enough iterations were alloted." << endl;
        return 0;
    }
    
    std::cout << "iter: " << iterCurrent << endl;
    return x0;
    
}