#include <vector>  // vector
#include <memory>  // shared_ptr
#include <cassert>

#include "Random.h"
#include "Particle.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"


// Getters
std::string Material::name()   { return m_name; }   // Name
// macroXsec
double Material::SigmaS( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaS( E ) * n.second;
	}	
	return sum;
}
double Material::SigmaC( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaC( E ) * n.second;
	}	
	return sum;
}
double Material::SigmaA( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaA( E ) * n.second;
	}	
	return sum;
}
double Material::SigmaF( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaF( E ) * n.second;
	}	
	return sum;
}
double Material::SigmaT( const double E ) 
{ 
	double sum = 0.0;
	for ( auto& n : nuclides )
	{
		sum += n.first->sigmaT( E ) * n.second;
	}	
	return sum;
}
double Material::nuSigmaF( const double E ) 
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
void Material::addNuclide( const std::shared_ptr< Nuclide_t >& Nuclide, double N ) 
{ nuclides.push_back( std::make_pair( Nuclide, N ) ); }


// Sample collision distance
double Material::collision_distance_sample( const double E )
{ return - std::log( Urand() ) / SigmaT( E ); }


// Sample collided nuclide
std::shared_ptr<Nuclide_t> Material::nuclide_sample( const double E )
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
std::shared_ptr<Nuclide_t> Material::nuclide_scatter( const double E )
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
std::shared_ptr<Nuclide_t> Material::nuclide_nufission( const double E )
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

void Material::collision_sample( Particle& P,std::stack<Particle>& Pbank,
                                   const bool ksearch, SourceBank& Fbank, 
                                   const double k )
{
    // Implicit Fission (k is always 1 in non ksearch modes)
    const double bank_nu = std::floor( P.weight() / k * nuSigmaF(P.energy()) 
                                       / SigmaT(P.energy()) + Urand() );

    std::shared_ptr<Nuclide_t> N_fission = nuclide_nufission( P.energy() );
    if (ksearch){ 
        for ( int i = 0 ; i < bank_nu ; i++ ){
            Particle P_new( P.pos(), isotropic.sample(),
                            N_fission->Chi(P.energy()), P.time(), 1.0, 0,
                            P.cell() );
            Fbank.add_source( std::make_shared<SourceDelta>(P_new), 1.0 );
        }
    } else{
        for ( int i = 0 ; i < bank_nu ; i++ ){
            Particle P_new( P.pos(), isotropic.sample(),
                            N_fission->Chi(P.energy()), P.time(), 1.0, 
                            P.tdmc(), P.cell() );
            Pbank.push(P_new);
        }
    }
    
    // Implicit Absorption
    const double implicit = SigmaC(P.energy()) + SigmaF(P.energy());
    P.set_weight( P.weight() * ( SigmaT(P.energy()) - implicit ) 
                 / SigmaT(P.energy()) );

    std::shared_ptr<Nuclide_t> N_scatter = nuclide_scatter( P.energy() );
    if(!N_scatter){ return; }
    
    std::shared_ptr<Reaction_t > R = N_scatter->scatter;
    R->sample( P, Pbank );
}
		

// Simulate scattering for scattering matrix MGXS
void Material::simulate_scatter( Particle& P )
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
