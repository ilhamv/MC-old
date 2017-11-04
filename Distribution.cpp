#include <cmath>
#include <vector>
#include <iostream>

#include "Distribution.h"
#include "Random.h"
#include "Point.h"
#include "Const.h"       // PI2
#include "Solver.h"
#include "Particle.h"
#include "Solver.h" // Linterpolate

double Uniform_Distribution::sample( const double param /*= 0.0*/ ) 
{
	return a + Urand() * range;
}


double Linear_Distribution::sample( const double param /*= 0.0*/ ) 
{
  	double r1 = Urand(), r2 = Urand();
  	double p  = 2.0 * std::fmin( fa, fb ) / ( fa + fb );
  	if ( r1 < p ) { return a + r2 * ( b - a ); }
  	else 
	{
    		if ( fb > fa ) { return a + ( b - a ) * std::sqrt( r2 ); }
    		else           { return a + ( b - a ) * ( 1.0 - std::sqrt( r2 )); }
  	}
}


double Cubic_Distribution::sample( const double param /*= 0.0*/ ) 
{
	double x;
	double y;
	
	do
	{
		x = a + Urand() * ( b - a );
		y = Urand() * fmax;
	}
	while ( y > f(x) );

	return x;
}


double Normal_Distribution::sample( const double param /*= 0.0*/ ) 
{
  	return mean + sigma * std::sqrt( -2.0 * std::log( Urand() ) ) * std::cos( twopi * Urand() ); 
}


double RayleighScatter_Distribution::sample( const double param /*= 0.0*/ ) 
{
  	// perform rejection sampling of Rayleigh scattering distribution bounded by uniform distribution
  	// empirical tests showed that this method was the fastest with optimized -O3 setting with g++ (on MacOS)
  	double x, y, f;  
  	do 
	{
    		x = 2.0  * Urand() - 1.0;
    		y = 0.75 * Urand();
    		f = 0.375 * ( 1.0 + x*x );
  	}
	while ( y > f );

  	return x;
}


double HGScatter_Distribution::sample( const double param /*= 0.0*/ ) 
{
  	// direct inversion of Henhey-Green scattering distribution
  	const double x = Urand() + E;

  	return A - D / ( x*x );
}


double LinearScatter_Distribution::sample( const double param /*= 0.0*/ )
{
	// Linear decomposition
	if ( Urand() < prob ) { return std::sqrt( 4.0 * Urand() ) - 1.0; }
	else { return 2.0 * Urand() - 1.0; }
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


Point_t IsotropicDirection_Distribution::sample( const double param /*= 0.0*/ )
{
	// Sample polar cosine and azimuthal angle uniformly
	const double mu  = 2.0 * Urand() - 1.0;
	const double azi = PI2 * Urand();

	// Convert to Cartesian coordinates
	double c = std::sqrt( 1.0 - mu * mu );
	Point_t p;
	p.y = std::cos( azi ) * c;
	p.z = std::sin( azi ) * c;
	p.x = mu;

  	return p;
}


Point_t AnisotropicDirection_Distribution::sample( const double param /*= 0.0*/ ) 
{
  	const double mu  = dist_mu->sample(); 
  	const double azi = PI2 * Urand();
  	const double cos_azi = std::cos(azi);
  	const double sin_azi = std::sin(azi);

  	// rotate the local particle coordinate system aligned along the incident direction
  	// to the global problem (x,y,z) coordinate system 
  	double sin_t0 = std::sqrt( 1.0 - mu * mu );
  	double c = sin_t0 / sin_t;

  	Point_t p;
  	p.x = axis.x * mu + ( axis.x * axis.z * cos_azi - axis.y * sin_azi ) * c;
  	p.y = axis.y * mu + ( axis.y * axis.z * cos_azi + axis.x * sin_azi ) * c;
  	p.z = axis.z * mu - cos_azi * sin_t0 * sin_t;
  	return p;
}


Point_t IndependentXYZ_Distribution::sample( const double param /*= 0.0*/ ) 
{
  	return Point_t( dist_x->sample(), dist_y->sample(), dist_z->sample() );
}
