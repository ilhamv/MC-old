#ifndef _SOLVER_HEADER_
#define _SOLVER_HEADER_

#include <vector>

#include "Const.h" // MAX
#include "Random.h"
#include "Point.h"


// Quad solver for geometry-point evaluation
// return smallest positive real root if it exists; if it does not, return very big number
double solve_quad( const double a, const double b, const double c );


// Binary search a double location in bin grids
// Note:
// 	value < lowest  grid --> -1
// 	value > highest grid --> vector.size - 1 (or number of bins)
// 	value = grid points  --> location of bin whose upper bound is the value
// 	                         (-1 if value = lowest grid)
int Binary_Search( const double x, const std::vector<double>& vec );


// Scatter direction
// return final direction dir_f after scatter initial direction dir_i with scattering polar cosine mu
Point_t scatter_direction( const Point_t dir_i, const double mu0 );


// Lagrenge interpolation
double Linterpolate( const double x, const double x1, const double x2, const double y1, const double y2 );


#endif
