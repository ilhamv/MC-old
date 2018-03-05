#ifndef _MATERIAL_H
#define _MATERIAL_H

#include <vector>
#include <memory>
#include <cmath> 
#include <stack> 
#include <cstring>

#include "Particle.h"
#include "Nuclide.h"
#include "Source.h"


class Material
{
    private:
	std::string m_name;
        DistributionIsotropicDirection isotropic;
	
        // Nuclides contained and its corresponding density
	std::vector< std::pair< std::shared_ptr<Nuclide_t>, double > > nuclides;        

    public:
     	// Constructor: pass the material name 
	Material( const std::string n ) : m_name(n) {};
	~Material() {};

	// Getters
	std::string name();// Name
	// macroXsec
	double      SigmaT( const double E );
	double      SigmaS( const double E );
	double      SigmaA( const double E );
	double      SigmaC( const double E );
	double      SigmaF( const double E );
	double    nuSigmaF( const double E );
	
	// Add a pair of nuclide and its total macroXs
	// the supplied variable are the nuclide and its nuclide density
	void addNuclide( const std::shared_ptr< Nuclide_t >& Nuclide, double N );

	// Sample collision distance
	double collision_distance_sample( const double E );
	
	// Sample collided nuclide
	std::shared_ptr<Nuclide_t> nuclide_sample( const double E );
	std::shared_ptr<Nuclide_t> nuclide_scatter( const double E );
	std::shared_ptr<Nuclide_t> nuclide_nufission( const double E );

	// Handle collision event
	// Sample entire collision (nuclide, then nuclide reaction)
	// Then, process the reaction on the Particle
	void collision_sample( Particle& P, std::stack< Particle >& Pbank, const bool ksearch, SourceBank& Fbank, const double k );
	
	// Simulate scattering for scattering matrix MGXS
	void simulate_scatter( Particle& P );
};


#endif // MATERIAL_H
