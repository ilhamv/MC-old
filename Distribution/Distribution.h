#ifndef _DISTRIBUTION_HEADER_
#define _DISTRIBUTION_HEADER_

#include <vector>   // vector
#include <cmath>    // erf, floor
#include <string>
#include <memory>
#include <cassert>
#include <fstream>

#include "Const.h" // PI
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


// Discrete distribution base class
template <class T>
class Discrete_Distribution : public Distribution_t<T> 
{
	protected:
		std::vector< std::pair< T, double > > cdf;

	public:
		 Discrete_Distribution( const std::vector< std::pair< T, double > >& data, const std::string label = "" );
		~Discrete_Distribution() {};
		T sample( const double param = 0.0 );
};


template < class T >
Discrete_Distribution<T>::Discrete_Distribution( const std::vector< std::pair< T, double > >& data, const std::string label /*= ""*/ ):
       Distribution_t<T>(label)	
{
	// convert pmf to cdf
	double c = 0.0;
	for ( auto& d : data ) 
	{
		// first is pointer to data type T and second is pmf input
		cdf.push_back( std::make_pair( d.first, d.second + c ) );
		c += d.second;
	}
    
}


template < class T >
T Discrete_Distribution<T>::sample( const double param /*= 0.0*/ ) 
{
	double   r = Urand() * cdf.back().second;
	for ( auto& c : cdf ) 
	{
		// first is pointer to data type T and second is cdf
		if ( r < c.second ) { return c.first; };
	}
    assert( false ); // should never reach here
    
}


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


// Linear distribution
class Linear_Distribution : public Distribution_t<double> 
{
  	private:
    		double a, b, fa, fb;
  	
	public:
    		Linear_Distribution( double x1, double x2, double y1, double y2, const std::string label = "" )  
      			: Distribution_t(label), a(x1), b(x2), fa(y1), fb(y2) {};
	    	~Linear_Distribution() {};
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
 

// Cubic distribution [rejection sampling]
class Cubic_Distribution : public Distribution_t<double>
{
  	private:
    		const double a, b, c3, c2, c1, c0, fmax;
  		double f( const double x )
		{
			const double x_sq = x*x;
			return c3 * x_sq*x + c2 * x_sq + c1 * x + c0;
		}	

	public:
    		Cubic_Distribution( const double p1, const double p2, const double p3, double const p4,
			       const double p5, const double p6, const double p7, const std::string label = "" )  
      			: Distribution_t(label), a(p1), b(p2), c3(p3), c2(p4), c1(p5), c0(p6), fmax(p7) {};
   		~Cubic_Distribution() {};
   		double sample( const double param = 0.0 );
};


// Rayleigh scattering distribution
class RayleighScatter_Distribution : public Distribution_t<double>
{
  	public:
     		 RayleighScatter_Distribution( const std::string label = "" ) : Distribution_t(label) {};
    		~RayleighScatter_Distribution() {};
    		double sample( const double param = 0.0 );
};


// Normal distribution with mean and standard deviation sigma
class Normal_Distribution : public Distribution_t<double>
{
  	private:
    		const double twopi = 2.0 * PI;
    		const double       mean, sigma;
  	
	public:
     		 Normal_Distribution( const double p1, const double p2, const std::string label = "" ) : Distribution_t(label), mean(p1), sigma(p2) {};
    		~Normal_Distribution() {};
	    	double sample( const double param = 0.0 );
};


// Henyey-Green scattering distribution with parameter g
class HGScatter_Distribution : public Distribution_t<double>
{
  	private:
    		const double g;
 		double A;
    		double B;
    		double C;
    		double D;
    		double E;

  	public:
     		HGScatter_Distribution( const double p, const std::string label = "" ) : g(p), Distribution_t(label)
		{
 			A = ( 1.0 + g*g ) / ( 2.0*g );
    			B = 2.0 * g / ( 1.0 - g*g );
    			C = 1.0 / sqrt( 2.0 * g + g*g + 1.0 );
			D = 1.0 / ( 2.0*g*B*B );
			E = C/B;
		};
    		~HGScatter_Distribution() {};
    		double sample( const double param = 0.0 );
};


// Isotropic scattering distribution
class IsotropicScatter_Distribution : public Distribution_t<double>
{
	public:
		 IsotropicScatter_Distribution( const std::string label = "" ) : Distribution_t(label) {};
		~IsotropicScatter_Distribution() {};

    		double sample( const double param = 0.0 ){ return 2.0 * Urand() - 1.0; }
};


// Linear scattering distribution
// using linear decomposition
class LinearScatter_Distribution : public Distribution_t<double>
{
	private:
		const double prob; // Probability for the first pdf

	public:
		 LinearScatter_Distribution( const double mubar, const std::string label = "" ) : prob( 3.0 * mubar ), Distribution_t(label) {};
		~LinearScatter_Distribution() {};

		double sample( const double param = 0.0 );
};


// Isotropic direction distribution
class IsotropicDirection_Distribution : public Distribution_t<Point_t>
{
	public:
		IsotropicDirection_Distribution( const std::string label = "" ) : Distribution_t(label) {};
    		~IsotropicDirection_Distribution() {};
    		
		Point_t sample( const double param = 0.0 );
};


// Average multiplicity distribution
class Average_Multiplicity_Distribution : public Distribution_t<int>
{
	public:
		 Average_Multiplicity_Distribution( const std::string label = "" ) : Distribution_t(label) {};
		~Average_Multiplicity_Distribution() {};
		
		int sample( const double nu = 0.0 ) { return std::floor( nu + Urand() ); }
};

// Terrel multiplicity distribution
class Terrel_Multiplicity_Distribution : public Discrete_Distribution<int>
{
	private:
		const double nubar, gamma, b, n;

	public:
		Terrel_Multiplicity_Distribution( const double p1, const double p2, const double p3, const int p4, const std::vector< std::pair< int, double > >& data, const std::string label = "" ) :
			nubar(p1), gamma(p2), b(p3), n(p4+1), Discrete_Distribution(data,label)
		{
			for ( int nu = 0; nu < n; nu++ )
			{
				double prob = 0.5 * ( erf( ( nu - nubar + 0.5 + b ) / ( gamma * sqrt(2.0) ) ) + 1.0 );
				cdf.push_back( std::make_pair( nu, prob ) );
			}
		}		
		~Terrel_Multiplicity_Distribution() {};
};


// Independent 3point distribution
class IndependentXYZ_Distribution : public Distribution_t<Point_t>
{
    	private:
        	std::shared_ptr<Distribution_t<double>> dist_x, dist_y, dist_z;
    
	public:
        	IndependentXYZ_Distribution( std::shared_ptr<Distribution_t<double>> dx,
                                std::shared_ptr<Distribution_t<double>> dy, std::shared_ptr<Distribution_t<double>> dz, const std::string label = "" )
        	: Distribution_t(label), dist_x(dx), dist_y(dy), dist_z(dz) {};
        	~IndependentXYZ_Distribution() {};
    
        	Point_t sample( const double param = 0.0 );
};


// Anisotropic direction distribution
class AnisotropicDirection_Distribution : public Distribution_t<Point_t>
{
    	private:
        	double sin_t;
    
        	Point_t axis;
        	std::shared_ptr<Distribution_t<double>> dist_mu;
    
    	public:
        	AnisotropicDirection_Distribution( Point_t p, std::shared_ptr<Distribution_t<double>> dmu, const std::string label = "" )
        	: Distribution_t(label), axis(p), dist_mu(dmu)
        	{ axis.normalize(); sin_t = std::sqrt( 1.0 - axis.z * axis.z ); };
        	~AnisotropicDirection_Distribution() {};
    
        	Point_t sample( const double param = 0.0 );
};


#endif
