double fzero_bisect(double (*func)(double), double bound_lower, double bound_upper, double err_tol, int iter_max);
double fzero_secant(double (*function)(double), double guess_init_1, double guess_init_2, double err_tol, int iter_max);
double fzero_false_pos(double (*function)(double), double bound_lower, double bound_upper, double err_tol, int iter_max);
double fzero_NR(double (*func)(double), double (*func_deriv)(double), double guess_init, double error_tol, int iter_max);
double fzero_Halley(double (*func)(double), double (*func_deriv)(double), double (*func_deriv_2)(double), double guess_init, double err_tol, int iter_max);