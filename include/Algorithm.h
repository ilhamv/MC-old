#ifndef _ALGORITHM_H
#define _ALGORITHM_H

#include <vector>

#include "Point.h"
#include "Geometry.h"
#include "Nuclide.h"


//=============================================================================
// Geometry 
//=============================================================================

double geometry_quad( const double a, const double b, const double c );


//=============================================================================
// Miscellany 
//=============================================================================

int binary_search( const double x, const std::vector<double>& vec );
Point scatter_direction( const Point dir_i, const double mu0 );
double interpolate( const double x, const double x1, const double x2,
                    const double y1, const double y2 );
void split_roulette( Particle& P, std::stack<Particle>& Pbank );
void normalize_point( Point& p );
bool point_equal( const Point p1, const Point p2 );

double exponential_sample( const double param );


#endif // ALGORITHM_H
