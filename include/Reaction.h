#ifndef _REACTION_HEADER_
#define _REACTION_HEADER_

#include <vector>  // vector
#include <memory>  // shared_ptr
#include <stack>   // stack
#include <cstring> // string

#include "Particle.h"
#include "Distribution.h"
#include "Point.h"
#include "XSec.h"

class Particle_t;

//=============================================================================
// Base
//=============================================================================

class Reaction_t
{
    private:
	std::shared_ptr<XSec_t> r_xs;
        const int r_type;
	
    protected:
	std::shared_ptr<XSec_t> r_nu;
	std::shared_ptr< Distribution_t<double> > Chi_dist;

    public:
	// Constructor: Pass the microXs
       	 Reaction_t( std::shared_ptr<XSec_t> x, const int t,
                     std::shared_ptr<XSec_t> n = std::make_shared<Constant_XSec>(0.0) ) : r_xs(x), r_type(t), r_nu(n) {}; // Pass the microXs
	~Reaction_t() {};

	// Get the microXs
	virtual double xs( const double E, const unsigned long long idx = 0 ) final { return r_xs->xs( E, idx ); };

	// Get the microXs
	virtual double nu( const double E, const unsigned long long idx = 0 ) final { return r_nu->xs( E, idx ); };

	// Sample the Chi spectrum
	virtual double Chi( const double E ) final { return Chi_dist->sample(E); }
	
	// Sample the reaction process on the working particle and the particle bank
	virtual void sample( Particle_t& P, std::stack< Particle_t >& Pbank ) = 0;
	
	// Check type
        //   0 -> Capture
        //   1 -> Scatter
        //   2 -> Fission
	virtual int type() final { return r_type; }
};


//=============================================================================
// Capture
//=============================================================================

class Capture_Reaction : public Reaction_t 
{
    public:
	// Constructor: Pass the microXs
	 Capture_Reaction( std::shared_ptr<XSec_t> x ) : Reaction_t(x,0) {}; // Pass the microXs
	~Capture_Reaction() {};

	// Kill the working particle upon reaction sample
	void  sample( Particle_t& P, std::stack< Particle_t >& Pbank );
};


//=============================================================================
// Scatter reaction, room temperature
//=============================================================================
class Scatter_Reaction : public Reaction_t 
{
    private:
	std::shared_ptr< Distribution_t<double> > scatter_dist; // Scattering angle distribution
	const double                              A;            // Nuclide mass
    public:
	// Constructor: Pass the microXs
	 Scatter_Reaction( std::shared_ptr<XSec_t> x, const std::shared_ptr< Distribution_t<double> >& D, const double a ) :
		Reaction_t(x,1), scatter_dist(D), A(a) {}; // Pass the microXs and scattering angle distribution
	~Scatter_Reaction() {};

	// Scatter the working particle
	void  sample( Particle_t& P, std::stack< Particle_t >& Pbank );
};


// Scatter reaction, zero K
class Scatter_Zero_Reaction : public Reaction_t 
{
    private:
	std::shared_ptr< Distribution_t<double> > scatter_dist; // Scattering angle distribution
	const double                              A;            // Nuclide mass
    public:
	// Constructor: Pass the microXs
	 Scatter_Zero_Reaction( std::shared_ptr<XSec_t> x, const std::shared_ptr< Distribution_t<double> >& D, const double a ) :
		Reaction_t(x,1), scatter_dist(D), A(a) {}; // Pass the microXs and scattering angle distribution
	~Scatter_Zero_Reaction() {};

	// Scatter the working particle
	void  sample( Particle_t& P, std::stack< Particle_t >& Pbank );
};


//=============================================================================
// Fission reaction
//=============================================================================
class Fission_Reaction : public Reaction_t 
{
    private:
	std::shared_ptr< Distribution_t<int> >    nu_dist;   // Fission multiplicity distribution
	IsotropicDirection_Distribution           isotropic; // Isotropic distribution for emerging neutron

    public:
	// Constructor: Pass the microXs and distributions
	 Fission_Reaction( std::shared_ptr<XSec_t> x, std::shared_ptr<XSec_t> n, const std::shared_ptr< Distribution_t<double> >& W  ) :
	 	Reaction_t(x,2,n)
	{
		Chi_dist = W;
	}

	~Fission_Reaction() {};

	// Sample fission multiplicity, then appropriately pushing new fission particles to the bank
	// --> Reaction.cpp
	void sample( Particle_t& P, std::stack< Particle_t >& Pbank );
};


#endif
