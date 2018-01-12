#include <vector> // vector
#include <memory> // shared_ptr
#include <cassert>
#include <iostream>

#include "Random.h"
#include "Nuclide.h"


// Getters
std::string Nuclide_t::name()   { return n_name; } // Name

// microXs
double Nuclide_t::sigmaS( const double E ) 
{ 
    checkE( E );
    for ( auto& r : reactions )
    {
	if ( r->type() == 1 ) { return r->xs( E, idx_help ); }
    }
    // Nuclide doesn't have the reaction
    return 0.0;
}
double Nuclide_t::sigmaC( const double E ) 
{ 
    checkE( E );
    for ( auto& r : reactions )
    {
        if ( r->type() == 0 ) { return r->xs( E, idx_help ); }
    }
    // Nuclide doesn't have the reaction
    return 0.0;
}
double Nuclide_t::sigmaF( const double E ) 
{ 
    checkE( E );
    for ( auto& r : reactions )
    {
	if ( r->type() == 2 ) { return r->xs( E, idx_help ); }
    }
    // Nuclide doesn't have the reaction
    return 0.0;
}
double Nuclide_t::sigmaA( const double E ) 
{ 
    checkE( E );
    double sum = 0.0;

    for ( auto& r : reactions )
    {
	if ( r->type() == 0 || r->type() == 2 ) { sum += r->xs( E, idx_help ); }
    }
    // Nuclide doesn't have the reaction
    return sum;
}
double Nuclide_t::sigmaT( const double E )
{ 
    checkE( E );
    double sum = 0.0;

    for ( auto& r : reactions )
    { 
	sum += r->xs( E, idx_help );
    }

    return sum; 
}
double Nuclide_t::nusigmaF( const double E ) 
{ 
    checkE( E );
    for ( auto& r : reactions )
    {
	if ( r->type() == 2 ) { return r->xs( E, idx_help ) * r->nu( E, idx_help ); }
    }
    // Nuclide doesn't have the reaction
    return 0.0;
}


// Check energy at cross secton call
//   if it's another different energy, search the location on the table --> idx_help
void Nuclide_t::checkE( const double E )
{
    if ( !E_table->empty() )
    {
	if ( E != E_current ) 
	{ 
	    idx_help  = Binary_Search( E, *E_table );
	    E_current = E;
	}
    }
}


// Sample Chi spectrum
double Nuclide_t::Chi( const double E )
{ 
    for ( auto& r : reactions )
    {
	if ( r->type() == 2 ) { return r->Chi(E); }
    }
    return 0.0;
};


// Set energy grids for table look-up XS
void Nuclide_t::setTable( const std::shared_ptr< std::vector<double> >& Evec )
{
    E_table = Evec;
}


// Add reaction
void Nuclide_t::addReaction( const std::shared_ptr< Reaction_t >& R ) 
{ 
    if ( R->type() == 1 ) { scatter = R; } // Attach pointer on scattering reaction
    reactions.push_back( R ); 
}


// Randomly sample a reaction type from the nuclide
std::shared_ptr< Reaction_t > Nuclide_t::reaction_sample( const double E ) 
{
    //Note: Implicit Capture/Absorption is implemented

    const double u = (sigmaT( E ) - sigmaC( E ) ) * Urand();

    double s = 0.0;
    for ( auto& r : reactions ) 
    {
	if ( !r->type() == 0 ) 
        {
            s += r->xs( E );
	    if ( s > u ) { return r; }
        }
    }
    return nullptr; // Purely capture, skip collision sample
}


// Simulate scattering for scattering matrix MGXS
void Nuclide_t::simulate_scatter( Particle_t& P )
{
    std::stack<Particle_t> null;
    scatter->sample(P,null);
}
