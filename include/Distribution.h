#ifndef _DISTRIBUTION_HEADER_
#define _DISTRIBUTION_HEADER_

#include <vector>   // vector
#include <cmath>    // erf, floor
#include <string>
#include <memory>
#include <cassert>
#include <fstream>

#include "Random.h" // Urand
#include "Point.h"
#include "XSec.h"

// Sampling Distributions base class
template< class T >
class Distribution_t
{
    private:
	const std::string d_name;

    public:
     	 Distribution_t( const std::string label = "" ) : d_name(label) {};
    	~Distribution_t() {};
		
	virtual std::string name() final { return d_name; };
		
	// Get sample value
	virtual T sample( const double param = 0.0 ) = 0;
};


// Delta distribution
template <class T>
class Delta_Distribution : public Distribution_t<T> 
{
    private:
    	T result;

    public:
     	 Delta_Distribution( T val, const std::string label = "" ) : Distribution_t<T>(label), result(val) {};
    	~Delta_Distribution() {};

    	T sample( const double param = 0.0 ) { return result; }
};


// Uniform distribution in [a,b]
class Uniform_Distribution : public Distribution_t<double>
{
    private:
    	const double a, b, range;

    public:
     	 Uniform_Distribution( const double p1, const double p2, const std::string label = "" ) : a( p1 ), b( p2 ), range( p2 - p1 ), Distribution_t(label) {};
    	~Uniform_Distribution() {};
    	
        double sample( const double param = 0.0 );
};


// Watt distribution - Hyperbolic sin model
class Watt_Distribution : public Distribution_t <double>
{
    private:
	// Vector of parameters corresponding to thermal(< 1eV), 1 MeV, and 14 MeV
	std::vector<double> vec_a;    
	std::vector<double> vec_b;
	std::vector<double> vec_g;

    public:
	// Constructor: Pass the parameters a and b and also compute g
	Watt_Distribution( const std::vector<double> p1, const std::vector<double> p2, const std::string label = "" ): 
	    Distribution_t(label), vec_a(p1), vec_b(p2)
	{
	    for ( int i = 0 ; i < vec_a.size() ; i++ )
	    {
		const double C    = ( 1.0 + vec_a[i]*vec_b[i]/8.0 );
		const double C_sq = C*C;
		vec_g.push_back( std::sqrt( C_sq - 1.0 ) + C );
	    }
    	}
        ~Watt_Distribution() {};
    		
	double sample( const double E = 0.0 );
};
 

// Isotropic scattering distribution
class IsotropicScatter_Distribution : public Distribution_t<double>
{
    public:
         IsotropicScatter_Distribution( const std::string label = "" ) : Distribution_t(label) {};
	~IsotropicScatter_Distribution() {};

    	double sample( const double param = 0.0 ){ return 2.0 * Urand() - 1.0; }
};


// Isotropic direction distribution
class IsotropicDirection_Distribution : public Distribution_t<Point>
{
    public:
	 IsotropicDirection_Distribution( const std::string label = "" ) : Distribution_t(label) {};
    	~IsotropicDirection_Distribution() {};
    		
	Point sample( const double param = 0.0 );
};


// Independent 3point distribution
class IndependentXYZ_Distribution : public Distribution_t<Point>
{
    private:
        std::shared_ptr<Distribution_t<double>> dist_x, dist_y, dist_z;
    
    public:
        IndependentXYZ_Distribution( std::shared_ptr<Distribution_t<double>> dx,
            std::shared_ptr<Distribution_t<double>> dy, std::shared_ptr<Distribution_t<double>> dz, const std::string label = "" )
                : Distribution_t(label), dist_x(dx), dist_y(dy), dist_z(dz) {};
        ~IndependentXYZ_Distribution() {};
    
        Point sample( const double param = 0.0 );
};


#endif
