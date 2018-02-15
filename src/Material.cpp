#include <vector>  // vector
#include <memory>  // shared_ptr
#include <cassert>

#include "Random.h"
#include "Particle.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"


// Getters
std::string Material_t::name()   { return m_name; }   // Name
// macroXsec
double Material_t::SigmaS( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaS( E ) * n.second;
	}	
	return sum;
}
double Material_t::SigmaC( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaC( E ) * n.second;
	}	
	return sum;
}
double Material_t::SigmaA( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaA( E ) * n.second;
	}	
	return sum;
}
double Material_t::SigmaF( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaF( E ) * n.second;
	}	
	return sum;
}
double Material_t::SigmaT( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaT( E ) * n.second;
	}	
	return sum;
}
double Material_t::nuSigmaF( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->nusigmaF( E ) * n.second;
	}	
	return sum;
}


// Add a nuclide
// the supplied variable are the nuclide and its nuclide density
void Material_t::addNuclide( const std::shared_ptr< Nuclide_t >& Nuclide, double N ) 
{ nuclides.push_back( std::make_pair( Nuclide, N ) ); }


// Sample collision distance
double Material_t::collision_distance_sample( const double E )
{ return - std::log( Urand() ) / SigmaT( E ); }


// Sample collided nuclide
std::shared_ptr< Nuclide_t > Material_t::nuclide_sample( const double E )
{
    double u = SigmaT( E ) * Urand();
    double s = 0.0;
    for ( auto& n : nuclides ){
	// first is pointer to nuclide, second is nuclide density
	s += n.first->sigmaT( E ) * n.second;
	if ( s > u ) { return n.first; }
    }
    //assert( false ); // should never reach here
    return nullptr;
}


//==============================================================================
// Collision
//==============================================================================

void Material_t::collision_sample( Particle_t& P, std::stack<Particle_t>& Pbank,
                                   const bool ksearch, Source_Bank& Fbank, 
                                   const double k )
{
    // Note that we implement Implicit Capture or Absorption (if ksearch)
    double implicit  = SigmaC(P.energy());

    // The implicit fission
    if (ksearch){ 
        implicit += SigmaF(P.energy()); 

        // Bank Fbank
        const double bank_nu = std::floor( P.weight() / k 
                                           * nuSigmaF(P.energy()) 
                                           / SigmaT(P.energy()) + Urand() );
        for ( int i = 0 ; i < bank_nu ; i++ )
        {
            // Determine the emitting nuclide 
            const double r = Urand();
            double       s = 0.0;
            for( auto& n : nuclides )
            {
                s += n.first->nusigmaF( P.energy() ) / nuSigmaF( P.energy() );
                if ( r < s )
                {
                    Fbank.addSource( 
                                    std::make_shared<Delta_Source>
                                    ( P.pos(), isotropic.sample(),
                                      n.first->Chi( P.energy() ), 1.0, 
                                      P.time() ) );
                    break;
                }
            }
        }
    }
    
    // The implicit capture/absorption
    P.setWeight( P.weight() * ( SigmaT(P.energy()) - implicit ) / SigmaT(P.energy()) );

    // First sample nuclide
    std::shared_ptr< Nuclide_t >  N = nuclide_sample( P.energy() );

    // Now get the reaction
    std::shared_ptr< Reaction_t > R = N->reaction_sample( P.energy(), ksearch );
	
    // Finally process the reaction on the Particle
    if( R ) { R->sample( P, Pbank ); }
}
		

// Simulate scattering for scattering matrix MGXS
void Material_t::simulate_scatter( Particle_t& P )
{
	// Sample the scattering nuclide
	double u = SigmaS( P.energy() ) * Urand();
	double s = 0.0;
	for ( auto& n : nuclides ) 
	{
		// first is pointer to nuclide, second is nuclide density
		s += n.first->sigmaS( P.energy() ) * n.second;
		if ( s > u ) 
		{ 
			// Simulate the scatter
			return n.first->simulate_scatter( P );
	       	}
	}
	// There is no scattering nuclide
}
