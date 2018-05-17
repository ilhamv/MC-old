#include <vector>
#include <iostream>

#include "Random.h"
#include "Algorithm.h"
#include "Constants.h"
#include "Point.h"


//=============================================================================
// Geometry 
//=============================================================================

// Return smallest positive real root if it exists; 
//   if it does not, return very big number
double geometry_quad( const double a, const double b, const double c ) 
{
    // Determinant
    const double D = b*b - 4.0 * a * c;
    	
    // roots are complex, no intersection, return huge number
    // or identical roots, tangent, return huge number
    if ( D <= 0.0 ) { return MAX_float; }
    else{
        const double sqrtD = std::sqrt(D);
        const double ai = 0.5 / a;

        // roots
        double r1 = ai * ( -1.0 * b - sqrtD );
        double r2 = ai * ( -1.0 * b + sqrtD );

        // Negative roots return huge number (moving away from surface)
        if ( r1 < 0 ) { r1 = MAX_float; }
        if ( r2 < 0 ) { r2 = MAX_float; }

        return std::fmin( r1, r2 );
    }
}


//=============================================================================
// Miscellany 
//=============================================================================

// Binary search a double location in a bin grid
int binary_search( const double x, const std::vector<double>& vec )
{
    int left  = 0;
    int right = vec.size() - 1;
    int mid;

    while ( left <= right ){
        mid = ( left +  right ) / 2;        
        if ( vec[mid] < x ) { left  = mid + 1; }
        else                { right = mid - 1; }
    }
	
    return right; 
    // Note:
    // 	value < lowest  grid --> -1
    // 	value > highest grid --> vector.size - 1 (or number of bins)
    // 	value = grid points  --> location of bin whose upper bound is the value
    // 	                         (-1 if value = lowest grid)
}

// Return final direction with scattering cosine mu
Point scatter_direction( const Point dir_i, const double mu0 )
{
    // Sample azimuthal direction
    const double     azi = PI_2 * Urand();
    const double cos_azi = std::cos(azi);
    const double sin_azi = std::sin(azi);
    const double      Ac = std::sqrt( 1.0 - mu0 * mu0 );
    Point      dir_f; // Final direction

    if( dir_i.z != 1.0 ){
        const double       B = std::sqrt( 1.0 - dir_i.z * dir_i.z );
        const double       C = Ac / B;
		
        dir_f.x = dir_i.x * mu0 
                  + ( dir_i.x * dir_i.z * cos_azi - dir_i.y * sin_azi ) * C;
        dir_f.y = dir_i.y * mu0 
                  + ( dir_i.y * dir_i.z * cos_azi + dir_i.x * sin_azi ) * C;
        dir_f.z = dir_i.z * mu0 - cos_azi * Ac * B;
    }
	
    // If dir_i = 0i + 0j + k, interchange z and y in the scattering formula
    else{
        const double       B = std::sqrt( 1.0 - dir_i.y * dir_i.y );
        const double       C = Ac / B;
		
        Point            q; // to store new direction point
        
        dir_f.x = dir_i.x * mu0 
                  + ( dir_i.x * dir_i.y * cos_azi - dir_i.z * sin_azi ) * C;
        dir_f.z = dir_i.z * mu0 
                  + ( dir_i.z * dir_i.y * cos_azi + dir_i.x * sin_azi ) * C;
        dir_f.y = dir_i.y * mu0 - cos_azi * Ac * B;
    }
    return dir_f;
}

double interpolate( const double x, const double x1, const double x2,
                    const double y1, const double y2 )
{ return ( x - x2 ) / ( x1 - x2 ) * y1 + ( x - x1 ) / ( x2 - x1 ) * y2; }

void normalize_point( Point& p )
{
    const double denom = p.x*p.x + p.y*p.y + p.z*p.z;
    p.x /= denom;
    p.y /= denom;
    p.z /= denom;
}

bool point_equal( const Point p1, const Point p2 )
{
    if( p1.x != p2.x ){ return false; }
    if( p1.y != p2.y ){ return false; }
    if( p1.z != p2.z ){ return false; }
    return true;
}

double exponential_sample( const double param )
{
    return -std::log(Urand()) / param;
}

