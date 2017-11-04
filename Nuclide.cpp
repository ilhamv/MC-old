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
		if ( r->type("scatter") ) { return r->xs( E, idx_help ); }
	}
	// Nuclide doesn't have the reaction
	return 0.0;
}
double Nuclide_t::sigmaC( const double E ) 
{ 
	checkE( E );
	for ( auto& r : reactions )
	{
		if ( r->type("capture") ) { return r->xs( E, idx_help ); }
	}
	// Nuclide doesn't have the reaction
	return 0.0;
}
double Nuclide_t::sigmaF( const double E ) 
{ 
	checkE( E );
	for ( auto& r : reactions )
	{
		if ( r->type("fission") ) { return r->xs( E, idx_help ); }
	}
	// Nuclide doesn't have the reaction
	return 0.0;
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
		if ( r->type("fission") ) { return r->xs( E, idx_help ) * r->nu( E, idx_help ); }
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
		if ( r->type( "fission" ) ) { return r->Chi(E); }
	}
	return 0.0;
};


// Set energy grids for table look-up XS
void Nuclide_t::setTable( const std::shared_ptr< std::vector<double> >& Evec )
{
	E_table = Evec;
}


// Add reaction
void Nuclide_t::addReaction( const std::shared_ptr< Reaction_t >& C ) 
{ 
	if ( C->type("scatter") ) { scatter = C; } // Attach pointer on scattering reaction
	reactions.push_back( C ); 
}


// Randomly sample a reaction type from the nuclide
std::shared_ptr< Reaction_t > Nuclide_t::reaction_sample( const double E ) 
{
	double u = sigmaT( E ) * Urand();
	double s = 0.0;
	for ( auto& r : reactions ) 
	{
		s += r->xs( E );
		if ( s > u ) { return r; }
	}
    assert( false ); // should never reach here
    return nullptr;
}


// Simulate scattering for scattering matrix MGXS
void Nuclide_t::simulate_scatter( Particle_t& P )
{
	std::stack<Particle_t> null;
	scatter->sample(P,null);
}
