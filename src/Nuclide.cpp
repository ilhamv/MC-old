#include <vector> 
#include <memory> 
#include <cassert>
#include <iostream>

#include "Random.h"
#include "Nuclide.h"
#include "Algorithm.h"

//==============================================================================
// Getters
//==============================================================================

std::string Nuclide::name() { return n_name; }
double Nuclide::sigmaS( const double E ) 
{ 
    checkE( E );
    for( auto& r : n_reactions )
    {
	if( r->type() == 1 ) { return r->xs( E, idx_help ); }
    }
    return 0.0;
}
double Nuclide::sigmaC( const double E ) 
{ 
    checkE( E );
    for( auto& r : n_reactions )
    {
        if( r->type() == 0 ) { return r->xs( E, idx_help ); }
    }
    return 0.0;
}
double Nuclide::sigmaF( const double E ) 
{ 
    checkE( E );
    for( auto& r : n_reactions )
    {
	if( r->type() == 2 ) { return r->xs( E, idx_help ); }
    }
    return 0.0;
}
double Nuclide::sigmaA( const double E ) 
{ 
    checkE( E );
    double sum = 0.0;
    for( auto& r : n_reactions )
    {
	if( r->type() == 0 || r->type() == 2 ){ 
            sum += r->xs( E, idx_help ); 
        }
    }
    return sum;
}
double Nuclide::sigmaT( const double E )
{ 
    checkE( E );
    double sum = 0.0;

    for( auto& r : n_reactions ){ sum += r->xs( E, idx_help ); }
    return sum; 
}
double Nuclide::nusigmaF( const double E ) 
{ 
    checkE( E );
    for ( auto& r : n_reactions )
    {
	if( r->type() == 2 ){ 
            return r->xs( E, idx_help ) * r->nu( E, idx_help ); 
        }
    }
    return 0.0;
}


// Check energy at cross secton call
//   if it's another different energy, 
//   search the location on the table --> idx_help
void Nuclide::checkE( const double E )
{
    if ( !E_table->empty() ){
	if ( E != E_current ){ 
	    idx_help  = binary_search( E, *E_table );
	    E_current = E;
	}
    }
}


// Sample Chi spectrum
double Nuclide::Chi( const double E )
{ 
    for ( auto& r : n_reactions )
    {
	if ( r->type() == 2 ) { return r->Chi(E); }
    }
    return 0.0;
};


// Set energy grids for table look-up XS
void Nuclide::setTable( const std::shared_ptr< std::vector<double> >& Evec )
{
    E_table = Evec;
}


// Add reaction
void Nuclide::addReaction( const std::shared_ptr< Reaction_t >& R ) 
{ 
    if ( R->type() == 1 ) { scatter = R; }
    n_reactions.push_back( R ); 
}


// Randomly sample a reaction type from the nuclide
std::shared_ptr< Reaction_t > Nuclide::reaction_sample( const double E, 
                                                          const bool ksearch )
{
    //Note: Implicit Capture/Absorption is implemented
    double         implicit  = sigmaC(E);
    if (ksearch) { implicit += sigmaF(E); }

    const double u = (sigmaT( E ) - implicit ) * Urand();

    double s = 0.0;
    for( auto& r : n_reactions ){
	if ( !r->type() == 0 && !( ksearch && r->type() == 2 ) ){
            s += r->xs( E );
	    if ( s > u ) { return r; }
        }
    }
    return nullptr;
}
