#include <stack> // stack

#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"


// Capture kills the working particle
void Capture_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{ P.kill(); }


// Scatter the working particle
void Scatter_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{ P.scatter( scatter_dist->sample(), A ); }


// Fission reaction sample
void Fission_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{
	// create random number of secondaries from multiplicity distributon nu_dist and
	// push all but one of them into the Particle bank, and reset the top particle 
	// if no secondaries, kill the particle

	int n = nu_dist->sample( r_nu->xs( P.energy() ) ); // sampled multiplicity
    
	if ( n != 0 ) 
    	{
        	// bank all but last particle (skips if n = 1)
	        for ( int i = 0 ; i < n - 1 ; i++ )
        	{
	            	Particle_t p( P.pos(), isotropic.sample(), Chi_dist->sample( P.energy() ), P.time(), P.weight() );
        	    	p.setCell( P.cell() );
            		Pbank.push( p );
        	}

		// reset the top particle
		P.setDirection( isotropic.sample() );
	}
	else
	{ P.kill(); }
}
   

// Check reaction type
bool Capture_Reaction::type( const std::string s ) { if( s.compare("capture") == 0 ) { return true; } return false; }
bool Scatter_Reaction::type( const std::string s ) { if( s.compare("scatter") == 0 ) { return true; } return false; }
bool Fission_Reaction::type( const std::string s ) { if( s.compare("fission") == 0 ) { return true; } return false; }
