#ifndef _GEOMETRY_HEADER_
#define _GEOMETRY_HEADER_

#include <memory>      // shared_ptr
#include <vector>      // vector
#include <cstring>     // string

#include "Const.h"     // EPSILON
#include "Particle.h"
#include "Estimator.h"
#include "Point.h"
#include "Material.h"


// Forward declaration
class Estimator_t;


// Geometry base class
class Geometry_t
{
	private:
		const std::string g_name;                                 // Geometry name
	
	protected:
		std::vector< std::shared_ptr< Estimator_t > > estimators; // Estimators attached
	
	public:
     		 Geometry_t( const std::string n ) : g_name(n) {};        // Pass the name
    		~Geometry_t() {};

    		// Get name
		virtual std::string name()  final { return g_name; };
		
		// Add an estimator
		virtual void addEstimator( const std::shared_ptr< Estimator_t >& E ) final 
		{ estimators.push_back( E ); }
};



////////////////
/// Surfaces ///
///////////////

// Forward declaration
class Cell_t;

// Surface base class
class Surface_t : public Geometry_t
{
	public:
		const std::string bc; // Boundary condition
		                      //   "transmission"
				      //   "reflective"
	
	protected:
		// Crossing the surface --> an epsilon kick to the working particle
		void cross ( Particle_t& P );

		// Reflect and cross the working particle
		virtual void reflect ( Particle_t& P ) = 0;
		
	public:
     		 Surface_t( const std::string n, const std::string b ) : Geometry_t(n), bc(b) {}; // Pass the name
    		~Surface_t() {};

		// Hit implementation
		virtual void hit( Particle_t& P, const std::vector<std::shared_ptr<Cell_t>>& Cell );

		// Evaluate point location via the "S" equation
		virtual double eval( const Point_t& p ) = 0;

		// Get moving particle distance to surface
    		virtual double distance( const Particle_t& P ) = 0;
};


// Plane-X
class PlaneX_Surface : public Surface_t
{
  	private:
    		const double x; // Parameters for S equation

	protected:
		void reflect ( Particle_t& P );

  	public:
     		 PlaneX_Surface( const std::string n, const std::string bc, const double loc ) : 
			 Surface_t(n,bc), x(loc) {};
    		~PlaneX_Surface() {};

		double eval    ( const Point_t& p );
     		double distance( const Particle_t& P );
};


// Plane-Y
class PlaneY_Surface : public Surface_t
{
  	private:
    		const double y; // Parameters for S equation

	protected:
		void reflect ( Particle_t& P );

  	public:
     		 PlaneY_Surface( const std::string n, const std::string bc, const double loc ) : 
			 Surface_t(n,bc), y(loc) {};
    		~PlaneY_Surface() {};

		double eval    ( const Point_t& p );
     		double distance( const Particle_t& P );
};


// Plane-Z
class PlaneZ_Surface : public Surface_t
{
  	private:
    		const double z; // Parameters for S equation

	protected:
		void reflect ( Particle_t& P );

  	public:
     		 PlaneZ_Surface( const std::string n, const std::string bc, const double loc ) : 
			 Surface_t(n,bc), z(loc) {};
    		~PlaneZ_Surface() {};

		double eval    ( const Point_t& p );
     		double distance( const Particle_t& P );
};


// Generic Plane
class Plane_Surface : public Surface_t
{
  	private:
    		const double a, b, c, d; // Parameters for S equation
		double modx, mody, modz; // Parameters for reflection

	protected:
		void reflect ( Particle_t& P );

  	public:
     		 Plane_Surface( const std::string n, const std::string bc, const double pa, const double pb, const double pc, const double pd ) : 
			 Surface_t(n,bc), a(pa), b(pb), c(pc), d(pd)
		{
			const double L = 2.0 / ( a*a + b*b + c*c );
			modx = L * a;
			mody = L * b;
			modz = L * c;
		};
    		~Plane_Surface() {};

		double eval    ( const Point_t& p );
     		double distance( const Particle_t& P );
};


// Sphere
class Sphere_Surface : public Surface_t 
{
	private:
    		const double x0, y0, z0, rad, rad_sq;
  	
	protected:
		void reflect ( Particle_t& P );

	public:
     		 Sphere_Surface( const std::string n, const std::string b, const double p1, const double p2, const double p3, const double p4 ) : 
       			Surface_t(n,b), x0(p1), y0(p2), z0(p3), rad(p4), rad_sq(p4*p4) {};
    		~Sphere_Surface() {};

     		double eval    ( const Point_t& p );
     		double distance( const Particle_t& P );
};


// Infinite Cylinder-X
class CylinderX_Surface : public Surface_t
{
	private:
		const double y0, z0, rad, rad_sq;

	protected:
		void reflect ( Particle_t& P );

	public:
		 CylinderX_Surface( const std::string n, const std::string b, const double p1, const double p2, const double p3 ) :
			Surface_t(n,b), y0(p1), z0(p2), rad(p3), rad_sq(p3*p3) {};
		~CylinderX_Surface() {};

		double eval    ( const Point_t& p );
		double distance( const Particle_t& P );
};


// Infinite Cylinder-Y
class CylinderY_Surface : public Surface_t
{
	private:
		const double x0, z0, rad, rad_sq;

	protected:
		void reflect ( Particle_t& P );

	public:
		 CylinderY_Surface( const std::string n, const std::string b, const double p1, const double p2, const double p3 ) :
			Surface_t(n,b), x0(p1), z0(p2), rad(p3), rad_sq(p3*p3) {};
		~CylinderY_Surface() {};

		double eval    ( const Point_t& p );
		double distance( const Particle_t& P );
};


// Infinite Cylinder-Z
class CylinderZ_Surface : public Surface_t
{
	private:
		const double x0, y0, rad, rad_sq;

	protected:
		void reflect ( Particle_t& P );

	public:
		 CylinderZ_Surface( const std::string n, const std::string b, const double p1, const double p2, const double p3 ) :
			Surface_t(n,b), x0(p1), y0(p2), rad(p3), rad_sq(p3*p3) {};
		~CylinderZ_Surface() {};

		double eval    ( const Point_t& p );
		double distance( const Particle_t& P );
};


// Infinite Cone-X 
class ConeX_Surface : public Surface_t
{
	private:
		const double x0, y0, z0, rad, rad_sq;

	protected:
		void reflect ( Particle_t& P );

	public:
		 ConeX_Surface( const std::string n, const std::string b, const double p1, const double p2, const double p3, const double p4 ) :
       			Surface_t(n,b), x0(p1), y0(p2), z0(p3), rad(p4), rad_sq(p4*p4) {};
		~ConeX_Surface() {};

		double eval    ( const Point_t& p );
		double distance( const Particle_t& P );
};


// Infinite Cone-Y
class ConeY_Surface : public Surface_t
{
	private:
		const double x0, y0, z0, rad, rad_sq;

	protected:
		void reflect ( Particle_t& P );

	public:
		 ConeY_Surface( const std::string n, const std::string b, const double p1, const double p2, const double p3, const double p4 ) :
       			Surface_t(n,b), x0(p1), y0(p2), z0(p3), rad(p4), rad_sq(p4*p4) {};
		~ConeY_Surface() {};

		double eval    ( const Point_t& p );
		double distance( const Particle_t& P );
};


// Infinite Cone-Z
class ConeZ_Surface : public Surface_t
{
	private:
		const double x0, y0, z0, rad, rad_sq;

	protected:
		void reflect ( Particle_t& P );

	public:
		 ConeZ_Surface( const std::string n, const std::string b, const double p1, const double p2, const double p3, const double p4 ) :
       			Surface_t(n,b), x0(p1), y0(p2), z0(p3), rad(p4), rad_sq(p4*p4) {};
		~ConeZ_Surface() {};

		double eval    ( const Point_t& p );
		double distance( const Particle_t& P );
};



//////////////
/// Cell ///
//////////////


// Cell base class
class Cell_t : public Geometry_t
{
	private:
		const double                                                  r_importance;    // Cell importance (for variance reduction)
		std::vector< std::pair< std::shared_ptr< Surface_t >, int > > surfaces;        // Pairs of bounding surface and its sense
		std::shared_ptr< Material_t >                                 material = NULL; // Contained material, default NULL means vacuum

	public:
     		 Cell_t( const std::string n, const double imp ) : // Pass name and importance
			 Geometry_t(n), r_importance(imp) {};
    		~Cell_t() {};

    		// Getters
		double importance(); // importance
		// MacroXsec of the contained material
		double SigmaT  ( const double E );
		double SigmaS  ( const double E );
		double SigmaC  ( const double E );
		double SigmaF  ( const double E );
		double SigmaA  ( const double E );
		double nuSigmaF( const double E );
		
		// Set the material
    		void setMaterial( const std::shared_ptr< Material_t >& M );
		
		// Add a bounding surface
		void addSurface ( const std::shared_ptr< Surface_t  >& S, const int sense );
    
        	// Return pointers to surfaces that belong to certain cell
        	std::vector< std::pair< std::shared_ptr< Surface_t >, int > > listSurfaces () { return surfaces; };

    		// Test if particle is in the cell
		bool testPoint( const Point_t& p );
		
		// Move particle and score any estimators
		void moveParticle( Particle_t& P, const double dmove );
		
		// Return the closest bounding surface and the corresponding particle hit distance 
		std::pair< std::shared_ptr< Surface_t >, double > surface_intersect( const Particle_t& P );
		
		// Return particle collision distance
		double collision_distance( const double E );

		// Let the Material take care of the collision sample and reaction process
		void collision( Particle_t& P, std::stack< Particle_t >& Pbank );

		// Simulate scattering for scattering matrix MGXS
		void simulate_scatter( Particle_t& P );
};


#endif
