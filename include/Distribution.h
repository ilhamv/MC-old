#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector> 
#include <string>
#include <memory>

#include "Point.h"


template<class T>
class Distribution
{
    private:
	const std::string d_name;
    public:
     	 Distribution( const std::string label = "" ) : d_name(label) {};
    	~Distribution() {};
		
	virtual std::string name() final { return d_name; };
	virtual T sample( const double param = 0.0 ) = 0;
};

template<class T>
class DistributionDelta : public Distribution<T>
{
    private:
    	const T result;
    public:
     	 DistributionDelta( T val, const std::string label = "" ): 
             Distribution<T>(label), result(val) {};
    	~DistributionDelta() {};
    	T sample( const double param = 0.0 ) { return result; }
};

class DistributionUniform : public Distribution<double>
{
    private:
    	const double a, b, range;
    public:
     	 DistributionUniform( const double p1, const double p2,
                              const std::string label = "" ): 
             a( p1 ), b( p2 ), range( p2 - p1 ), Distribution(label) {};
    	~DistributionUniform() {};
        double sample( const double param = 0.0 );
};

class DistributionWatt : public Distribution <double>
{
    private:
	// Vector of parameters corresponding to 
        //   thermal(< 1eV), 1 MeV, and 14 MeV
	std::vector<double> vec_a;    
	std::vector<double> vec_b;
	std::vector<double> vec_g;
    public:
	DistributionWatt( const std::vector<double> p1, 
                          const std::vector<double> p2, 
                          const std::string label = "" );
        ~DistributionWatt() {};
	double sample( const double E = 0.0 );
};
 
class DistributionIsotropicScatter : public Distribution<double>
{
    public:
         DistributionIsotropicScatter( const std::string label = "" ): 
             Distribution(label) {};
	~DistributionIsotropicScatter() {};
    	double sample( const double param = 0.0 );
};

class DistributionIsotropicDirection : public Distribution<Point>
{
    public:
	 DistributionIsotropicDirection( const std::string label = "" ): 
             Distribution(label) {};
    	~DistributionIsotropicDirection() {};
	Point sample( const double param = 0.0 );
};

class DistributionIndepndentXYZ : public Distribution<Point>
{
    private:
        std::shared_ptr<Distribution<double>> dist_x, dist_y, dist_z;
    public:
        DistributionIndepndentXYZ( std::shared_ptr<Distribution<double>> dx,
                                   std::shared_ptr<Distribution<double>> dy, 
                                   std::shared_ptr<Distribution<double>> dz, 
                                   const std::string label = "" ): 
            Distribution(label), dist_x(dx), dist_y(dy), dist_z(dz) {};
        ~DistributionIndepndentXYZ() {};
        Point sample( const double param = 0.0 );
};


class DistributionDelayedNeutron : public Distribution<double>
{
    private:
        const std::vector<double> CDF;
        const std::vector<double> v;
    public:
     	 DistributionDelayedNeutron( const std::vector<double> cdf, 
                              const std::vector<double> val,
                              const std::string label = "" ):
             CDF(cdf), v(val), Distribution(label) {};
    	~DistributionDelayedNeutron() {};
        double sample( const double param = 0.0 );
};

#endif // DISTRIBUTION_H
