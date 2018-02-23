#include <cmath>
#include <vector>
#include <iostream>

#include "Distribution.h"
#include "Random.h"
#include "Point.h"
#include "Const.h"       // PI2
#include "Solver.h"
#include "Solver.h" // Linterpolate

double Uniform_Distribution::sample( const double param /*= 0.0*/ ) 
{
    return a + Urand() * range;
}

double Watt_Distribution::sample( const double E /*= 0.0*/ )
{
    double a;
    double b;
    double g;
    double xi;   // xi_1 in formula
    double C;    // Acceptance parameter
    double Eout;

    // Binary search is not employed as there are only three grid points
    if ( E <= 1.0 )
    {
        // E <= 1 eV (thermal)
	a = vec_a[0];
	b = vec_b[0];
	g = vec_g[0];
    }
    else if ( E <= 1.0e6 )
    {
	// 1 eV < E <= 1 MeV
	a = Linterpolate( E, 1.0 , 1.0e6, vec_a[0], vec_a[1] );
	b = Linterpolate( E, 1.0 , 1.0e6, vec_b[0], vec_b[1] );
	g = Linterpolate( E, 1.0 , 1.0e6, vec_g[0], vec_g[1] );
    }
    else
    {
	// E >= 1 MeV, note for E > 14 MeV the values are extrapolated
	a = Linterpolate( E, 1.0e6, 14.0e6, vec_a[1], vec_a[2] );
	b = Linterpolate( E, 1.0e6, 14.0e6, vec_b[1], vec_b[2] );
	g = Linterpolate( E, 1.0e6, 14.0e6, vec_g[1], vec_g[2] );
    }
    
    do
    {
        xi = Urand();
	Eout  = -a*g * std::log( xi ); //MeV
        C = ( 1.0 - g ) * ( 1.0 - std::log( xi ) ) - std::log( Urand() );
    }
    while ( C*C > b*Eout );
    
    return ( Eout*1.0e6 ); //eV
}


Point IsotropicDirection_Distribution::sample( const double param /*= 0.0*/ )
{
    // Sample polar cosine and azimuthal angle uniformly
    const double mu  = 2.0 * Urand() - 1.0;
    const double azi = PI2 * Urand();

	
    // Convert to Cartesian coordinates
    double c = std::sqrt( 1.0 - mu * mu );
    Point p;
    p.y = std::cos( azi ) * c;
    p.z = std::sin( azi ) * c;
    p.x = mu;

    return p;
}

Point IndependentXYZ_Distribution::sample( const double param /*= 0.0*/ ) 
{
    return Point( dist_x->sample(), dist_y->sample(), dist_z->sample() );
}
