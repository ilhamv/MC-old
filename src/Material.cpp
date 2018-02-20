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
std::shared_ptr<Nuclide_t> Material_t::nuclide_sample( const double E )
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

// Sample scattered nuclide
std::shared_ptr<Nuclide_t> Material_t::nuclide_scatter( const double E )
{
    double u = SigmaS( E ) * Urand();
    double s = 0.0;
    for ( auto& n : nuclides ){
	// first is pointer to nuclide, second is nuclide density
	s += n.first->sigmaS( E ) * n.second;
	if ( s > u ) { return n.first; }
    }
    //assert( false ); // should never reach here
    return nullptr;
}

// Sample nu-fission nuclide
std::shared_ptr<Nuclide_t> Material_t::nuclide_nufission( const double E )
{
    double u = nuSigmaF( E ) * Urand();
    double s = 0.0;
    for ( auto& n : nuclides ){
	// first is pointer to nuclide, second is nuclide density
	s += n.first->nusigmaF( E ) * n.second;
	if ( s > u ) { return n.first; }
    }
    //assert( false ); // should never reach here
    return nullptr;
}

//=============================================================================
// Collision
//=============================================================================

void Material_t::collision_sample( Particle_t& P,std::stack<Particle_t>& Pbank,
                                   const bool ksearch, Source_Bank& Fbank, 
                                   const double k )
{
    // Implicit Fission (k is always 1 in non ksearch modes)
    const double bank_nu = std::floor( P.weight() / k * nuSigmaF(P.energy()) 
                                       / SigmaT(P.energy()) + Urand() );

    std::shared_ptr<Nuclide_t> N_fission = nuclide_nufission( P.energy() );
    if (ksearch){ 
        for ( int i = 0 ; i < bank_nu ; i++ ){
            Fbank.addSource( std::make_shared<Delta_Source>
                    ( P.pos(), isotropic.sample(), N_fission->Chi(P.energy()),
                      1.0, P.time() ) );
        }
    } else{
        for ( int i = 0 ; i < bank_nu ; i++ ){
            Particle_t P_new ( P.pos(), isotropic.sample(),
                                    N_fission->Chi( P.energy() ), P.time(), 
                                    1.0, P.tdmc() );
            P_new.setCell( P.cell() );
            Pbank.push(P_new);
        }
    }
    
    // Implicit Absorption
    const double implicit = SigmaC(P.energy()) + SigmaF(P.energy());
    P.setWeight( P.weight() * ( SigmaT(P.energy()) - implicit ) 
                 / SigmaT(P.energy()) );

    std::shared_ptr<Nuclide_t> N_scatter = nuclide_scatter( P.energy() );
    if(!N_scatter){ return; }
    
    std::shared_ptr<Reaction_t > R = N_scatter->scatter;
    R->sample( P, Pbank );
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
