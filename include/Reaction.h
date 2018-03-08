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
        const std::vector<std::shared_ptr<Distribution<double>>> r_ChiD;
        const std::shared_ptr<XS> r_beta;
        const std::vector<double> r_lambda;
        const std::vector<double> r_fraction;

    public:
	// Constructor: Pass the microXs and distributions
	 ReactionFission( const std::shared_ptr<XS> x, 
                          const std::shared_ptr<XS> n, 
                          const std::shared_ptr<Distribution<double>> W,
                    const std::vector<std::shared_ptr<Distribution<double>>> D,
                          const std::shared_ptr<XS> b,
                          const std::vector<double> l,
                          const std::vector<double> f) :Reaction(x), 
         r_nu(n), r_Chi(W), r_ChiD(D), r_beta(b), r_lambda(l), r_fraction(f) {};

	~ReactionFission() {};


	double nu( const unsigned long long idx, const double E,
                   const std::vector<double>& E_vec );
        double Chi( const double E );
        double ChiD( const int g, const double E );
	double beta( const unsigned long long idx, const double E,
                     const std::vector<double>& E_vec );
        double lambda( const int g );
        double fraction( const int g );
};


#endif // REACTION_H
