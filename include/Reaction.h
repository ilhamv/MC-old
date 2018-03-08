#ifndef _REACTION_H
#define _REACTION_H

#include <vector>  
#include <memory> 
#include <stack>   
#include <cstring> 

#include "Particle.h"
#include "Distribution.h"
#include "Point.h"
#include "XSec.h"


//=============================================================================
// Basic Reaction
//=============================================================================

class Reaction
{
    private:
	std::shared_ptr<XS> r_xs;

    public:
       	 Reaction( const std::shared_ptr<XS> x ) : r_xs(x) {};
	~Reaction() {};

	double xs( const unsigned long long idx, const double E,
                   const std::vector<double>& E_vec );
};


//=============================================================================
// Scatter reaction
//=============================================================================

class ReactionScatter : public Reaction 
{
    private:
	const std::shared_ptr< Distribution<double> > scatter_dist; 
	const double A;
    public:
	 ReactionScatter( const std::shared_ptr<XS> x, 
                          const std::shared_ptr<Distribution<double>>& D, 
                          const double a ): 
             Reaction(x), scatter_dist(D), A(a) {}; 
	~ReactionScatter() {};

	void  sample( Particle& P );
};


//=============================================================================
// Fission reaction
//=============================================================================
class ReactionFission : public Reaction 
{
    private:
	const std::shared_ptr<XS> r_nu;
        const std::shared_ptr<Distribution<double>> r_Chi;

    public:
	// Constructor: Pass the microXs and distributions
	 ReactionFission( const std::shared_ptr<XS> x, 
                          const std::shared_ptr<XS> n, 
                          const std::shared_ptr<Distribution<double>>& W  ) :
	 	Reaction(x), r_nu(n), r_Chi(W) {};

	~ReactionFission() {};


	double nu( const unsigned long long idx, const double E,
                   const std::vector<double>& E_vec );
        double Chi( const double E );
};


#endif // REACTION_H
