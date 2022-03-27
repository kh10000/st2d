/*
 Code for calculating symmetry function value
 For that, this code recieve the atomic coordinates and calculate nearest neighbor.
 */

#include <math.h>
//#include "mpi.h"
#include "symmetry_functions.h"

extern "C" int calculate_chebyshev(double **, double **, double **,
                                   int *, double *, int, int*, int,
                                   int*, double *, int, 
                                   double**);
