#include <cmath>

#include "Entropy.h"
#include "Algorithm.h"


//=============================================================================
// Construction: initialize p
//=============================================================================

ShannonEntropy::ShannonEntropy( const std::vector<double>& vx,
                                const std::vector<double>& vy,
                                const std::vector<double>& vz ): 
                                x(vx),y(vy),z(vz),Ix(vx.size() - 1 ),
                                Iy(vy.size() - 1 ),Iz(vz.size() - 1 )
{
    I = Ix*Iy*Iz;
    Izy = Iz*Iy;
    p.resize(I, 0.0);
}


//=============================================================================
// Score p
//=============================================================================

void ShannonEntropy::score( const Point& po, const int N )
{
    const int ix = binary_search( po.x, x );
    const int iy = binary_search( po.y, y );
    const int iz = binary_search( po.z, z );

    int idx = ix * Izy + iy * Iz + iz;
    p[idx] += N;
}


//=============================================================================
// Calculate entropy
//=============================================================================


double ShannonEntropy::calculate_H()
{
    double sum = 0.0;

    // Normalize p
    for( int i = 0; i < p.size(); i++ ) { sum += p[i]; }
    if( sum == 0.0 ) { return sum; }
    for( int i = 0; i < p.size(); i++ ) { p[i] /= sum; }

    // Calculate H
    sum = 0;
    for( int i = 0; i < p.size(); i++ ){
        if( p[i] != 0 ){ sum -= p[i] * std::log2( p[i] ); }
    }

    // Reset p
    std::fill( p.begin(), p.end(), 0.0 ); 

    return sum;
}
