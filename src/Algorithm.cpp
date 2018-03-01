#include <vector>
#include <iostream>

#include "Constants.h" // MAX
#include "Algorithm.h"
#include "Point.h"

// Quad solver for geometry-point evaluation
// return smallest positive real root if it exists; if it does not, return very big number
double solve_quad( const double a, const double b, const double c ) 
{
    // Determinant
    const double D = b*b - 4.0 * a * c;
    	
    // roots are complex, no intersection, return huge number
    // or identical roots, tangent, return huge number
    if ( D <= 0.0 ) { return MAX_float; }
  	
    else 
    {
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


// Binary search a double location in a bin grid
int Binary_Search( const double x, const std::vector<double>& vec )
{
    int left  = 0;
    int right = vec.size() - 1;
    int mid;

    while ( left <= right )
    {
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


// Scatter direction
// Return final direction dir_f after scatter initial direction dir_i with scattering cosine mu
Point scatter_direction( const Point dir_i, const double mu0 )
{
    // Sample azimuthal direction
    const double     azi = PI_2 * Urand();
    const double cos_azi = std::cos(azi);
    const double sin_azi = std::sin(azi);
    const double      Ac = std::sqrt( 1.0 - mu0 * mu0 );
    Point      dir_f; // Final direction

    if( dir_i.z != 1.0 )
    {
        const double       B = std::sqrt( 1.0 - dir_i.z * dir_i.z );
        const double       C = Ac / B;
		
        dir_f.x = dir_i.x * mu0 + ( dir_i.x * dir_i.z * cos_azi - dir_i.y * sin_azi ) * C;
        dir_f.y = dir_i.y * mu0 + ( dir_i.y * dir_i.z * cos_azi + dir_i.x * sin_azi ) * C;
        dir_f.z = dir_i.z * mu0 - cos_azi * Ac * B;
    }
	
    // If dir_i = 0i + 0j + k, interchange z and y in the scattering formula
    else
    {
        const double       B = std::sqrt( 1.0 - dir_i.y * dir_i.y );
        const double       C = Ac / B;
		
        Point            q; // to store new direction point
        
        dir_f.x = dir_i.x * mu0 + ( dir_i.x * dir_i.y * cos_azi - dir_i.z * sin_azi ) * C;
        dir_f.z = dir_i.z * mu0 + ( dir_i.z * dir_i.y * cos_azi + dir_i.x * sin_azi ) * C;
        dir_f.y = dir_i.y * mu0 - cos_azi * Ac * B;
    }
    return dir_f;
}

// Lagrenge interpolation
double Linterpolate( const double x, const double x1, const double x2, const double y1, const double y2 )
{ return ( x - x2 ) / ( x1 - x2 ) * y1 + ( x - x1 ) / ( x2 - x1 ) * y2; }


// Shannon entropy
void Shannon_Entropy_Mesh::clear()
{
    total = 0;
    for ( int i = 0; i < x_nmesh; i++ ){
        for ( int j = 0; j < y_nmesh; j++ ){
	    for ( int k = 0; k < z_nmesh; k++ ){
		mesh[i][j][k] = 0.0;
	    }
	}
    }
    return;
}

void Shannon_Entropy_Mesh::update( const Point& p, const double N )
{
    int x_index = floor( ( ( p.x - xmin ) * x_nmesh ) / (xmax - xmin) );
    int y_index = floor( ( ( p.y - ymin ) * y_nmesh ) / (ymax - ymin) );
    int z_index = floor( ( ( p.z - zmin ) * z_nmesh ) / (zmax - zmin) );

    mesh[x_index][y_index][z_index] += N;
    total += N;
    return;
}

double Shannon_Entropy_Mesh::entropy()
{
    double e = 0.0;
    for ( int i = 0; i < x_nmesh; i++ ){
	for ( int j = 0; j < y_nmesh; j++ ){
	    for ( int k = 0; k < z_nmesh; k++ ){
		const double p = mesh[i][j][k] / total;
		e -= p * std::log2(p);
	    }
	}
    }
    return e;
}
std::shared_ptr<Cell> search_cell( const Point& p, const std::vector<std::shared_ptr<Cell>>& Cell )
{
    for( const auto& C : Cell ){
        if ( C->testPoint( p ) ){ return C; }
    }
    std::cout<< "[WARNING] A particle is lost:\n( x, y, z )  (" << p.x << ", " << p.y << ", " << p.z << " )\n";
    std::exit(EXIT_FAILURE);
}
