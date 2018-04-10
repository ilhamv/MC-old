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
	const std::string m_name;
	const std::vector<std::pair<std::shared_ptr<Nuclide>,double>>nuclides;

    public:
	Material( const std::string n, 
          const std::vector<std::pair<std::shared_ptr<Nuclide>,double>> nd ):
            m_name(n), nuclides(nd) {};
	~Material() {};

	// Getters
	std::string name();
	double   SigmaT( const double E );
	double   SigmaS( const double E );
	double   SigmaA( const double E );
	double   SigmaC( const double E );
	double   SigmaF( const double E );
	double nuSigmaF( const double E );
	double nuSigmaF_prompt( const double E );
	double nuSigmaF_delayed( const double E, const int i);
	double nuSigmaF_delayed_decay( const double E, const int i);
	
	// Sample collided nuclide
	std::shared_ptr<Nuclide> nuclide_sample( const double E );
	std::shared_ptr<Nuclide> nuclide_scatter( const double E );
	std::shared_ptr<Nuclide> nuclide_nufission( const double E );
	std::shared_ptr<Nuclide> nuclide_nufission_prompt( const double E );
	std::shared_ptr<Nuclide> nuclide_nufission_delayed( const double E,
                                                            const int i );
};


#endif // MATERIAL_H
