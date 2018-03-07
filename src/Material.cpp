#include <vector>  // vector
#include <memory>  // shared_ptr
#include <cassert>

#include "Random.h"
#include "Particle.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"


//==============================================================================
// Getters
//==============================================================================

std::string Material::name() { return m_name; }
double Material::SigmaS( const double E ) 
{ 
    double sum = 0.0;
    for ( auto& n : nuclides ){
	sum += n.first->sigmaS( E ) * n.second;
    }	
    return sum;
}
double Material::SigmaC( const double E ) 
{ 
    double sum = 0.0;
    for ( auto& n : nuclides ){
	sum += n.first->sigmaC( E ) * n.second;
    }	
    return sum;
}
double Material::SigmaA( const double E ) 
{ 
    double sum = 0.0;
    for ( auto& n : nuclides ){
	sum += n.first->sigmaA( E ) * n.second;
    }	
    return sum;
}
double Material::SigmaF( const double E ) 
{ 
    double sum = 0.0;
    for ( auto& n : nuclides ){
        sum += n.first->sigmaF( E ) * n.second;
    }	
    return sum;
}
double Material::SigmaT( const double E ) 
{ 
    double sum = 0.0;
    for ( auto& n : nuclides ){
	sum += n.first->sigmaT( E ) * n.second;
    }	
    return sum;
}
double Material::nuSigmaF( const double E ) 
{ 
    double sum = 0.0;
    for ( auto& n : nuclides ){
	sum += n.first->nusigmaF( E ) * n.second;
    }	
    return sum;
}

//==============================================================================
// Sample collided nuclide
//==============================================================================

std::shared_ptr<Nuclide> Material::nuclide_sample( const double E )
{
    double u = SigmaT( E ) * Urand();
    double s = 0.0;
    for ( auto& n : nuclides ){
	s += n.first->sigmaT( E ) * n.second;
	if ( s > u ) { return n.first; }
    }
    return nullptr;
}
std::shared_ptr<Nuclide> Material::nuclide_scatter( const double E )
{
    double u = SigmaS( E ) * Urand();
    double s = 0.0;
    for ( auto& n : nuclides ){
	s += n.first->sigmaS( E ) * n.second;
	if ( s > u ) { return n.first; }
    }
    return nullptr;
}
std::shared_ptr<Nuclide> Material::nuclide_nufission( const double E )
{
    double u = nuSigmaF( E ) * Urand();
    double s = 0.0;
    for ( auto& n : nuclides ){
	s += n.first->nusigmaF( E ) * n.second;
	if ( s > u ) { return n.first; }
    }
    return nullptr;
}
