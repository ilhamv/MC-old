#include <cmath>      // floor
#include <memory>     // shared_ptr
#include <stack>      // stack
#include <iostream>

#include "Random.h"
#include "Particle.h"
#include "Geometry.h"


// Perform Splitting & Rouletting variance reduction technique
// arguments: Old working cell, working particle, particle bank
void Split_Roulette( std::shared_ptr<Cell_t>& C, Particle_t& P, std::stack<Particle_t>& Pbank )
{
	const double Iold = C->importance();          // Previous importance
	const double Inew = P.cell()->importance(); // Current importance
	
	// Update working cell
	C = P.cell();

	// Same importance, do nothing
	if ( Inew == Iold ) { return; }
	
	const double rat = Inew / Iold; // Ratio of importances
	
	// Less important, Russian Roulette
	if ( rat < 1.0 )
	{
		if ( Urand() < rat ) { P.setWeight( P.weight() / rat ); }
		else                 { P.kill(); }
	}

	// More important, Splitting
	else
	{
		// Sample # of splitting
		const int n = std::floor( rat + Urand() );
		
		// Update working particle weight
		P.setWeight( P.weight() / double(n) );

		// Push n-1 identical particles into particle bank
		for ( int i = 0 ; i < n-1 ; i++ )
		{
			Pbank.push( P );
		}
	}
}

