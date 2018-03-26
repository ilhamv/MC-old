#include <vector> 
#include <memory> 
#include <cassert>
#include <iostream>

#include "Random.h"
#include "Nuclide.h"
#include "Algorithm.h"


//==============================================================================
// XS Binary search tool
//==============================================================================

// Check energy at cross secton call
//   if it's another different energy, 
//   search the location on the table --> idx
void Nuclide::checkE( const double E )
{
    if ( E != E_current ){ 
	idx  = binary_search( E, n_E );
	E_current = E;
    }
}


//==============================================================================
// Getters
//==============================================================================

std::string Nuclide::name() { return n_name; }
std::shared_ptr<ReactionScatter> Nuclide::scatter() { return n_scatter; }
std::shared_ptr<ReactionFission> Nuclide::fission() { return n_fission; }

double Nuclide::sigmaS( const double E ){
    checkE( E );
    return n_scatter->xs(idx, E, n_E);
}
double Nuclide::sigmaC( const double E ){
    checkE( E );
    return n_capture->xs(idx, E, n_E);
}
double Nuclide::sigmaF( const double E ){
    checkE( E );
    return n_fission->xs(idx, E, n_E);
}
double Nuclide::sigmaA( const double E ) 
{ 
    checkE( E );
    return n_absorb->xs(idx, E, n_E);
}
double Nuclide::sigmaT( const double E ) 
{ 
    checkE( E );
    return n_total->xs(idx, E, n_E);
}
double Nuclide::nusigmaF( const double E ) 
{ 
    checkE( E );
    return n_fission->xs(idx, E, n_E) * n_fission->nu(idx, E, n_E);
}
double Nuclide::nusigmaF_prompt( const double E ) 
{ 
    checkE( E );
    return ( 1.0 - beta(E) ) 
           * n_fission->xs(idx, E, n_E) * n_fission->nu(idx, E, n_E);
}
double Nuclide::nusigmaF_delayed( const double E, const int i ) 
{ 
    checkE( E );
    return beta(E) * n_fission->fraction(i)
           * n_fission->xs(idx, E, n_E) * n_fission->nu(idx, E, n_E);
}
double Nuclide::beta( const double E ){
    checkE( E );
    return n_fission->beta(idx, E, n_E);
}
