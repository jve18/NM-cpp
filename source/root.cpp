
#include "root.h"
#include <iostream>
#include <complex>
#include <cmath>

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
        std::cout << "(root.cpp, fzero_bisect) Error: " << "Input lower bound equal to input upper bound. Invalid input. Disregard return value." << std::endl;
        return 0;
    }
    
    // Error handling for entering a lower bound greater than the upper bound.
    try{
        if(bound_lower > bound_upper){
            throw 10002;
        }
        
    }
    catch(int err10002){
        std::cout << "(root.cpp, fzero_bisect) Warning: " << "Input lower bound greater than input upper bound. Evaluate reasonability of return value." << std::endl;
        
    }
    
    // Initialize lower bound, function at lower bound.
    double x0 = bound_lower;
    double f0 = func(x0);

    // Initialize upper bound, function at upper bound.
    double x2 = bound_upper;
    double f2 = func(x2);

    // Initialize mid-point, function at mid-point.
    double x1 = (x0 + x2)/2;
    double f1 = func(x1);
    
    // Initialize error and iteration count
    double err = f1;
    int iter = 0;

    std::cout << std::abs(err) << std::endl;
    
    while(std::abs(err) > err_tol && iter <= iter_max)
    {
        // Re-evaluate function at upper and lower bounds.
        f0 = func(x0);
        f2 = func(x2);

        // Re-evaluate midpoint and function at mid-point.
        x1 = (x0 + x2)/2; 
        f1 = func(x1);

        // Reset lower/upper bounds based on sign of f0, f2.
        if(f0 * f1 > 0)
        {
            x0 = x1; 
        } 
        else
        {
            x2 = x1;
        }
        
        // Compute error, increment iteration.
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
        std::cout << "(root.cpp, fzero_bisect) Error: " << "Maximum iteration limit exceeded. A high number of iterations may be necessary ~ log2((bound_upper-bound_lower)/err_tol)" << std::endl;
        return 0;
    }
   
    std::cout << "(root.cpp, fzero_bisect) iter: " << iter << std::endl;
    return x1;
    
}

double fzero_secant(double (*function)(double), double guess_init_1, double guess_init_2, double err_tol, int iter_max){
     /*
     INPUTS:
     *  function: The function that the user would like to find the root (0) of
     *  guess_init_1: First initial guess near the root. DOES NOT have to bracket the root. DOES NOT have to be less than guess_init_2
     *  guess_init_2: Second initial guess near the root. DOES NOT have to bracket the root. DOES NOT have to be greater than guess_init_1
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
        if(guess_init_1 == guess_init_2){
            throw 10001;
        }
        
    }
    catch(int err10001){
        std::cout << "(root.cpp, fzero_secant) Error: " << "Input first initial guess equal to second initial guess. Invalid input. Disregard return value." << std::endl;
        return 0;
    }
    
    // Initialize first guess, function evaluated at first guess.
    double x0 = guess_init_1;
    double f0 = function(x0);

    // Initialize second guess, function evaluated at second guess.
    double x1 = guess_init_2;
    double f1 = function(x1);

    // Initialize third guess, function evualuated at third guess.
    double x2 = (x0*f1 - x1*f0)/(f1 - f0); 
    double f2 = function(x2);

    // Initialize error and iteration count    
    double err = 10.0;
    int iter = 0;
    
    while(std::abs(err) > err_tol && iter <= iter_max){

        // Evaluate function at new first/second guesses.
        f0 = function(x0);
        f1 = function(x1);

        // Compute new third guess and evaluate function.
        x2 = (x0*f1 - x1*f0)/(f1 - f0);
        f2 = function(x2);

        // Set new first/second guesses.
        x0 = x1;
        x1 = x2;

        // Evaluate error, increment iteration.
        err = f2;
        iter++;
        
    }
    
    // Error handling for if max number of iterations  
    try{
        if(iter >= iter_max){
            throw 10000;
        }
        
    }
    catch(int err10000){
        std::cout << "(root.cpp, fzero_secant) Error: " << "Maximum iteration limit exceeded." << std::endl;
        return 0;
    }
    
    std::cout << "(root.cpp, fzero_secant) iter: " << iter << std::endl;
    return x2;
    
}

double fzero_false_pos(double (*function)(double), double bound_lower, double bound_upper, double err_tol, int iter_max){
    /*
     INPUTS:
     *  function: The function that the user would like to find the root (0) of
     *  bound_lower: The lower bound of the domain the root is known to exist within
     *  bound_upper: The upper bound of the domain the root is known to exist within
     *  err_tol: The amount by which the function evaluated at the best estimate of the root is allowed to differ from 0 by; f(x1); how far from 0 can f(x1) be
     *  iter_max: Maximum number of iterations for the regula falsi algorithm before it stops
     
     OUTPUTS:
     *  x1: The approximated root of the function
     
     SOURCES:
     *  https://en.wikipedia.org/wiki/False_position_method
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
        std::cout << "(root.cpp, fzero_false_pos) Error: " << "Input lower bound equal to input upper bound. Invalid input. Disregard return value." << std::endl;
        return 0;
    }
    
    // Error handling for entering a lower bound greater than the upper bound.
    try{
        if(bound_lower > bound_upper){
            throw 10002;
        }
        
    }
    catch(int err10002){
        std::cout << "(root.cpp, fzero_false_pos) Warning: " << "Input lower bound greater than input upper bound. Evaluate reasonability of return value." << std::endl;
        
    }
    
    // Initialize lower bound, function at lower bound.
    double x0 = bound_lower;
    double f0 = function(x0);

    // Initialize upper bound, function at upper bound.
    double x2 = bound_upper;
    double f2 = function(x2);

    // Initialize false position point, function at false position point.
    double x1 = x2 - ((f2*(x2-x0))/(f2-f0));
    double f1 = function(x1); 
    
    // Initialize error and iteration count.
    double err = 10;
    int iter = 0;
    
    while(std::abs(err) > err_tol && iter <= iter_max)
    {
        
        // Re-evaluate function at new lower and upper bounds.
        f0 = function(x0);
        f2 = function(x2);

        // Compute new false position point, and evaluate function.
        x1 = x2 - ((f2*(x2-x0))/(f2-f0));
        f1 = function(x1);
        
        // Set lower/upper bound for next iteration.
        if(f0 * f1 > 0)
        {
            x0 = x1;
        } 
        else
        {
            x2 = x1;
        }
        
        // Evaluate error, increment iteration.
        err = f1;
        iter++;
        
    }
    
    // Error handling for if max number of iterations.
    try{
        if(iter >= iter_max){
            throw 10000;
        }
        
    }
    catch(int err10000){
        std::cout << "(root.cpp, fzero_false_pos) Error: " << "Maximum iteration limit exceeded" << std::endl;
        return 0;
    }
    
    std::cout << "(root.cpp, fzero_false_pos) iter: " << iter << std::endl;
    return x1;   
    
}

double fzero_NR(double (*func)(double), double (*func_deriv)(double), double guess_init, double err_tol, int iter_max){
    /*
     INPUTS:
     *  func: The function that the user would like to find the root (0) of
     *  functionDerivative: The 1st derivative of the function that the user would like to find the root (0) of
     *  guess_init: The initial guess for where the root might be
     *  err_tol: The amount by which the function evaluated at the best estimate of the root is allowed to differ from 0 by; f(x1); how far from 0 can f(x1) be
     *  iter_max: Maximum number of iterations for the Newton-Raphson algorithm before it stops
     
     OUTPUTS:
     *  x0: The approximated root of the function
     
     SOURCES:
     *  https://en.wikipedia.org/wiki/Newton%27s_method
     *  TAMU AERO 220 Summer 2016 Lecture Notes
     *  Numerical Methods in Engineering with MATLAB by Jaan Kiusalaas (3rd ed.)
     */
    
    // Initialize guess, function at guess, function derivative at guess.
    double x0 = guess_init;
    double f0 = func(x0);
    double fp0 = func_deriv(x0);
    
    // Error handling for entering an initial guess that causes fp0 to equal 0.
    try{
        if(fp0 == 0){
            throw 10000;
        }
        
    }
    catch(int err10000){
        std::cout << "(root.cpp, fzero_NR) Error: " << "Function first derivative evaluated at initial guess equal to 0. Invalid input. Disregard return value." << std::endl;
        return 0;
    }
    
    // Initialize error and iteration count.
    double err = f0;
    int iter = 0;
    
    while(std::abs(err) > err_tol && iter <= iter_max)
    {
        // Evaluate function and function derivative at guess.
        f0 = func(x0);
        fp0 = func_deriv(x0);

        // New guess.
        x0 = x0 - f0/fp0;

        // Evaluate error, increment iteration.
        err = f0;
        iter++;
    }
    
    // Error handling for if max number of iterations.
    try{
        if(iter >= iter_max){
            throw 10001;
        }
        
    }
    catch(int err10001){
        std::cout << "(root.cpp, fzero_NR) Error: " << "Maximum iteration limit exceeded" << std::endl;
        return 0;
    }
    
    std::cout << "(root.cpp, fzero_NR) iter: " << iter << std::endl;
    return x0;
    
}

double fzero_Halley(double (*func)(double), double (*func_deriv)(double), double (*func_deriv2)(double), double guess_init, double err_tol, int iter_max){
     /*
     INPUTS:
     *  func: The function that the user would like to find the root (0) of
     *  func_deriv: The 1st derivative of the function that the user would like to find the root (0) of
     *  func_deriv2: The 2nd derivative of the function that the user would like to find the root (0) of
     *  guess_init: The initial guess for where the root might be
     *  err_tol: The amount by which the function evaluated at the best estimate of the root is allowed to differ from 0 by; f(x1); how far from 0 can f(x1) be
     *  iter_max: Maximum number of iterations for the Newton-Raphson algorithm before it stops
     
     OUTPUTS:
     *  x0: The approximated root of the function
     
     SOURCES:
     *  https://en.wikipedia.org/wiki/Newton%27s_method
     *  TAMU AERO 220 Summer 2016 Lecture Notes
     *  Numerical Methods in Engineering with MATLAB by Jaan Kiusalaas (3rd ed.)
     */
    
     // Initialize guess, function at guess, function derivative at guess, function second derivative at guess.
    double x0 = guess_init;
    double f0 = func(x0);
    double fp0 = func_deriv(x0);
    double fpp0 = func_deriv2(x0);

    // Initialize error and iteration count.
    double err = f0;
    int iter = 0;
    
    while(abs(err) > err_tol && iter <= iter_max)
    {
        // Evaluate function, function derivative, function second derivative at guess.
        f0 = func(x0);
        fp0 = func_deriv(x0);
        fpp0 = func_deriv2(x0);

        // New guess.
        x0 = x0 - (2*f0*fp0)/(2*fp0*fp0-f0*fpp0);

        // Evaluate error, increment iteration.
        err = f0;
        iter++;
    }
    
    // Error handling for if max number of iterations  
    try{
        if(iter >= iter_max){
            throw 10001;
        }
        
    }
    catch(int err10001){
        std::cout << "(root.cpp, fzero_Halley) Error: " << "Maximum iteration limit exceeded" << std::endl;
        return 0;
    }
    
    std::cout << "(root.cpp, fzero_Halley) iter: " << iter << std::endl;
    return x0;
    
}