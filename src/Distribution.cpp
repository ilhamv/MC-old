#include <cmath>

#include "Distribution.h"
#include "Random.h"
#include "Point.h"
#include "Constants.h"
#include "Algorithm.h"


//=============================================================================
// Constructor
//=============================================================================

DistributionWatt::DistributionWatt( const std::vector<double> p1, 
                                    const std::vector<double> p2, 
                                    const std::string label /*= ""*/ ):
    Distribution(label), vec_a(p1), vec_b(p2)
{
    for( int i = 0 ; i < vec_a.size() ; i++ ){
        const double C    = ( 1.0 + vec_a[i]*vec_b[i]/8.0 );
	vec_g.push_back( std::sqrt( C*C - 1.0 ) + C );
    }
}


//=============================================================================
// Sample
//=============================================================================

double DistributionUniform::sample( const double param /*= 0.0*/ ) 
{
    return a + Urand() * range;
}
double DistributionWatt::sample( const double E /*= 0.0*/ )
{
    double a;
    double b;
    double g;
    double xi;   // xi_1 in formula
    double C;    // Acceptance parameter
    double Eout;

    // Binary search is not employed as there are only three grid points
    // E <= 1 eV (thermal)
    if ( E <= 1.0 ){
	a = vec_a[0];
	b = vec_b[0];
	g = vec_g[0];
    }
    // 1 eV < E <= 1 MeV
    else if ( E <= 1.0e6 )
    {
	a = interpolate( E, 1.0 , 1.0e6, vec_a[0], vec_a[1] );
	b = interpolate( E, 1.0 , 1.0e6, vec_b[0], vec_b[1] );
	g = interpolate( E, 1.0 , 1.0e6, vec_g[0], vec_g[1] );
    }
    // E >= 1 MeV, note for E > 14 MeV the values are extrapolated
    else
    {
	a = interpolate( E, 1.0e6, 14.0e6, vec_a[1], vec_a[2] );
	b = interpolate( E, 1.0e6, 14.0e6, vec_b[1], vec_b[2] );
	g = interpolate( E, 1.0e6, 14.0e6, vec_g[1], vec_g[2] );
    }
    
    do{
        xi = Urand();
	Eout  = -a*g * std::log( xi ); //MeV
        C = ( 1.0 - g ) * ( 1.0 - std::log( xi ) ) - std::log( Urand() );
    }
    while ( C*C > b*Eout );
    
    return ( Eout*1.0e6 ); //eV
}
double DistributionIsotropicScatter::sample( const double param /*= 0.0*/ )
{ 
    return 2.0 * Urand() - 1.0; 
}
Point DistributionIsotropicDirection::sample( const double param /*= 0.0*/ )
{
    // Sample polar cosine and azimuthal angle uniformly
    const double mu  = 2.0 * Urand() - 1.0;
    const double azi = PI_2 * Urand();
	
    // Convert to Cartesian coordinates
    double c = std::sqrt( 1.0 - mu * mu );
    Point p;
    p.y = std::cos( azi ) * c;
    p.z = std::sin( azi ) * c;
    p.x = mu;

    return p;
}
Point DistributionIndepndentXYZ::sample( const double param /*= 0.0*/ ) 
{
    return Point( dist_x->sample(), dist_y->sample(), dist_z->sample() );
}
