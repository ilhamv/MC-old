#ifndef _ALGORITHM_H
#define _ALGORITHM_H

#include <vector>

#include "Point.h"
#include "Geometry.h"


double geometry_quad( const double a, const double b, const double c );
int binary_search( const double x, const std::vector<double>& vec );
Point scatter_direction( const Point dir_i, const double mu0 );
double interpolate( const double x, const double x1, const double x2,
                    const double y1, const double y2 );
void split_roulette( Particle& P, std::stack<Particle>& Pbank );

std::shared_ptr<Cell> search_cell( const Point& p,
                            const std::vector<std::shared_ptr<Cell>>& Cell );


#endif // ALGORITHM_H
