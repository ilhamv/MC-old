#ifndef _SOURCE_HEADER_
#define _SOURCE_HEADER_

#include <vector>     // vector
#include <memory>     // shared_ptr

#include "Particle.h"
#include "Const.h"    // EPSILON
#include "Geometry.h"


// Particle source base class
class Source_t
{
  	public:
 		 Source_t() {};
		~Source_t() {};

		// Get the particle source
		virtual Particle_t getSource() = 0;
};


// Point Source
class Point_Source : public Source_t
{
	private:
		Point_t                         pos;       // Position
    		// Direction, energy and time distribution
    		const std::shared_ptr< Distribution_t<Point_t> > dist_dir;
    		const std::shared_ptr< Distribution_t<double>  > dist_enrg;
    		const std::shared_ptr< Distribution_t<double>  > dist_time;

	public:
		Point_Source( const double p1, const double p2, const double p3, const std::shared_ptr< Distribution_t<Point_t> > dir
				,const std::shared_ptr< Distribution_t<double> > enrg, const std::shared_ptr< Distribution_t<double> > time  ) :
		dist_dir(dir), dist_enrg(enrg), dist_time(time) { pos.x = p1; pos.y = p2; pos.z = p3; }
		~Point_Source() {};

		// Get the particle source
		Particle_t getSource();
};


// DiskX Source
class DiskX_Source : public Source_t
{
	private:
		// Center position radius and direction sense of the source
		const double x0, y0, z0, r;
    		// Direction, energy and time distribution
    		const std::shared_ptr< Distribution_t<Point_t> > dist_dir;
    		const std::shared_ptr< Distribution_t<double>  > dist_enrg;
    		const std::shared_ptr< Distribution_t<double>  > dist_time;

	public:
		 DiskX_Source( const double p1, const double p2, const double p3, const double p4, const std::shared_ptr< Distribution_t<Point_t> > dir
			,const std::shared_ptr< Distribution_t<double> > enrg, const std::shared_ptr< Distribution_t<double> > time  ) :
			x0(p1), y0(p2), z0(p3), r(p4), dist_dir(dir), dist_enrg(enrg), dist_time(time) {};
		~DiskX_Source() {};

		// Get the particle source with rejection sampling
		// Direct methods are way too costly (cos, sin, sqrt)
		// Acceptance probability is pretty good ~3.14/4
		Particle_t getSource();
};


// DiskZ Source
class DiskZ_Source : public Source_t
{
	private:
		// Center position radius and direction sense of the source
		const double x0, y0, z0, r;
    		// Direction, energy and time distribution
    		const std::shared_ptr< Distribution_t<Point_t> > dist_dir;
    		const std::shared_ptr< Distribution_t<double>  > dist_enrg;
    		const std::shared_ptr< Distribution_t<double>  > dist_time;

	public:
		 DiskZ_Source( const double p1, const double p2, const double p3, const double p4, const std::shared_ptr< Distribution_t<Point_t> > dir
			,const std::shared_ptr< Distribution_t<double> > enrg, const std::shared_ptr< Distribution_t<double> > time  ) :
			x0(p1), y0(p2), z0(p3), r(p4), dist_dir(dir), dist_enrg(enrg), dist_time(time) {};
		~DiskZ_Source() {};

		// Get the particle source with rejection sampling
		// Direct methods are way too costly (cos, sin, sqrt)
		// Acceptance probability is pretty good ~3.14/4
		Particle_t getSource();
};


// Spherical Shell Source
class Sphere_Shell_Source : public Source_t
{
	private:
		// Center position and outer radius
		const double x0, y0, z0, ro;
		// Inner radius normalized by the outer radius, then squared
		const double risq;
    		// Direction, energy and time distribution
    		const std::shared_ptr< Distribution_t<Point_t> > dist_dir;
    		const std::shared_ptr< Distribution_t<double>  > dist_enrg;
    		const std::shared_ptr< Distribution_t<double>  > dist_time;

	public:
		 Sphere_Shell_Source( const double p1, const double p2, const double p3, const double p4, const double p5, const std::shared_ptr< Distribution_t<Point_t> > dir
			,const std::shared_ptr< Distribution_t<double> > enrg, const std::shared_ptr< Distribution_t<double> > time  ) :
			x0(p1), y0(p2), z0(p3), risq( p4*p4/p5/p5 ), ro(p5), dist_dir(dir), dist_enrg(enrg), dist_time(time) {};
		~Sphere_Shell_Source() {};

		// Get the particle source with rejection sampling
		// Direct methods are way too costly (acos, sin, sqrt)
		// Acceptance probability is almost half, ~3.14/6, yet it is still faster
		Particle_t getSource();
};


// Generic source
class Generic_Source : public Source_t
{
  	private:
    		// Position, direction, energy and time distribution
		const std::shared_ptr< Distribution_t<Point_t> > dist_pos;
    		const std::shared_ptr< Distribution_t<Point_t> > dist_dir;
    		const std::shared_ptr< Distribution_t<double>  > dist_enrg;
    		const std::shared_ptr< Distribution_t<double>  > dist_time;
  	public:
     		Generic_Source( const std::shared_ptr< Distribution_t<Point_t> > pos, const std::shared_ptr< Distribution_t<Point_t> > dir
				,const std::shared_ptr< Distribution_t<double> > enrg, const std::shared_ptr< Distribution_t<double> > time )
			: dist_pos(pos), dist_dir(dir), dist_enrg(enrg), dist_time(time) {};
    		~Generic_Source() {};
    		Particle_t getSource();
};


// Source Bank
// A collection of sources and its probability (or ratio) 
// An interface for every individual sources to the simulation
class Source_Bank
{
    private:
        std::vector< std::pair< std::shared_ptr<Source_t>, double > >  sources;
        double                                                         total = 0.0; // total probability (or ratio)
	
    public:
        Source_Bank() {};
        ~Source_Bank() {};
		
        void addSource( const std::shared_ptr<Source_t>& S, const double prob = 1.0 )
        { 
            sources.push_back( std::make_pair( S, prob ) );
            total += prob;
        }
		
        // Get source
        // sources are sampled wrt to their probability
        // then, particle cell is searched and set
        Particle_t getSource( const std::vector<std::shared_ptr<Cell_t>>& Cell )
        {
            const double xi = total * Urand();
            double s  = 0.0;
            for ( auto& So : sources ) 
            {
                // first is source, second is ratio
                s += So.second;
                if ( s > xi ) 
                { 
                    Particle_t P = So.first->getSource();
                    P.searchCell( Cell );
                    return P;
                }
            }
            //this is added because there is a possibility that this class does not return anything.
            std::cout<< "[ERROR] Source weights are not normalized to one\n";
            std::exit(EXIT_FAILURE);
        }
};


#endif
