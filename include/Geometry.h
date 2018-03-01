#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <memory> 
#include <vector> 
#include <stack>
#include <cstring>

#include "Constants.h"
#include "Particle.h"
#include "Point.h"
#include "Estimator.h"

class Estimator;
class SourceBank;
class Material_t;


//=============================================================================
// Geometry
//=============================================================================

class Geometry
{
    private:
	const std::string g_name;
        const int         g_ID;
        
    public:
     	Geometry( const std::string n, const int i ) : g_name(n), g_ID(i) {};
    	~Geometry() {};

	// Getters
	virtual std::string name() final { return g_name; }
	virtual int ID() final { return g_ID; }	
};


//=============================================================================
// Surface
//=============================================================================

class Surface : public Geometry
{
    private:
	const int s_bc; // 0:"transmission"
		        // 1:"reflective"
		
    public:
     	Surface( const std::string n, const int i, const int b ): 
            Geometry(n,i), s_bc(b) {}; // Pass the name
    	~Surface() {};

        virtual int bc() final { return s_bc; }
        
        std::vector<std::shared_ptr<Estimator>> estimators;
	virtual void attach_estimator( const std::shared_ptr<Estimator>& E ) 
            final{ estimators.push_back( E ); }

	// Type specific
	virtual double eval( const Point& p ) = 0;
    	virtual double distance( const Particle& P ) = 0;
        virtual void reflect( Particle& P ) = 0;
};

// Plane-X
class SurfacePlaneX : public Surface
{
    private:
    	const double x;
    public:
     	SurfacePlaneX( const std::string n, const int i, const int bc, 
                        const double loc ) : 
            Surface(n,i,bc), x(loc) {};
    	~SurfacePlaneX() {};	
        double eval( const Point& p );
     	double distance( const Particle& P );
	void reflect ( Particle& P );
};

// Plane-Y
class SurfacePlaneY : public Surface
{
    private:
    	const double y;
    public:
     	SurfacePlaneY( const std::string n, const int i, const int bc,
                        const double loc ): 
            Surface(n,i,bc), y(loc) {};
    	~SurfacePlaneY() {};
	double eval( const Point& p );
     	double distance( const Particle& P );
        void reflect ( Particle& P );
};

// Plane-Z
class SurfacePlaneZ : public Surface
{
    private:
    	const double z;
    public:
     	SurfacePlaneZ( const std::string n, const int i, const int bc, 
                        const double loc ) : 
            Surface(n,i,bc), z(loc) {};
    	~SurfacePlaneZ() {};
	double eval( const Point& p );
     	double distance( const Particle& P );
	void reflect( Particle& P );
};

// Generic Plane
class SurfacePlane : public Surface
{
    private:
    	const double a, b, c, d; 
	double modx, mody, modz; 
    public:
     	SurfacePlane( const std::string n, const int i, const int bc, 
                       const double pa, const double pb, const double pc, 
                       const double pd );
    	~SurfacePlane() {};
	double eval( const Point& p );
     	double distance( const Particle& P );
	void reflect( Particle& P );
};

// Sphere
class SurfaceSphere : public Surface 
{
    private:
    	const double x0, y0, z0, rad, rad_sq;

    public:
     	SurfaceSphere( const std::string n, const int i, const int b, 
                        const double p1, const double p2, const double p3, 
                        const double p4 ): 
       			Surface(n,i,b), x0(p1), y0(p2), z0(p3), rad(p4), rad_sq(p4*p4) {};
    	~SurfaceSphere() {};

     	double eval( const Point& p );
     	double distance( const Particle& P );
	void reflect( Particle& P );
};

// Infinite Cylinder-X
class SurfaceCylinderX : public Surface
{
    private:
	const double y0, z0, rad, rad_sq;
    public:
	SurfaceCylinderX( const std::string n, const int i, const int b, 
                           const double p1, const double p2, const double p3 ):
            Surface(n,i,b), y0(p1), z0(p2), rad(p3), rad_sq(p3*p3) {};
	~SurfaceCylinderX() {};
	double eval( const Point& p );
	double distance( const Particle& P );
	void reflect( Particle& P );
};

// Infinite Cylinder-Y
class SurfaceCylinderY : public Surface
{
    private:
        const double x0, z0, rad, rad_sq;
    public:
	SurfaceCylinderY( const std::string n, const int i, const int b, 
                           const double p1, const double p2, const double p3 ):
            Surface(n,i,b), x0(p1), z0(p2), rad(p3), rad_sq(p3*p3) {};
	~SurfaceCylinderY() {};
	double eval( const Point& p );
	double distance( const Particle& P );
	void reflect( Particle& P );
};

// Infinite Cylinder-Z
class SurfaceCylinderZ : public Surface
{
    private:
	const double x0, y0, rad, rad_sq;
    public:
	SurfaceCylinderZ( const std::string n, const int i, const int b, 
                           const double p1, const double p2, const double p3 ):
            Surface(n,i,b), x0(p1), y0(p2), rad(p3), rad_sq(p3*p3) {};
	~SurfaceCylinderZ() {};
	double eval( const Point& p );
	double distance( const Particle& P );
	void reflect( Particle& P );
};


//=============================================================================
// Cell
//=============================================================================

class Cell : public Geometry
{
    private:
	const double r_importance;
	std::vector<std::pair<std::shared_ptr<Surface>, int>> surfaces;
        std::shared_ptr<Material_t> c_material = NULL;

    public:
     	Cell( const std::string n, const int i, const double imp ) : // Pass name and importance
            Geometry(n,i), r_importance(imp) {};
    	~Cell() {};

	// Getters
	double importance();
	double SigmaT  ( const double E );
	double SigmaS  ( const double E );
	double SigmaC  ( const double E );
	double SigmaF  ( const double E );
	double SigmaA  ( const double E );
	double nuSigmaF( const double E );
        std::shared_ptr<Material_t> material();
	
	// Set the material
	void setMaterial( const std::shared_ptr< Material_t >& M );
	
	// Add a bounding surface
	void addSurface ( const std::shared_ptr< Surface  >& S, const int sense );

	// Return pointers to surfaces that belong to certain cell
	std::vector< std::pair< std::shared_ptr< Surface >, int > > listSurfaces () { return surfaces; };

	// Test if particle is in the cell
	bool testPoint( const Point& p );
	
	// Move particle and score any estimators
	void moveParticle( Particle& P, const double dmove, const bool tally );
	
	// Return the closest bounding surface and the corresponding particle hit distance 
	std::pair< std::shared_ptr< Surface >, double > surface_intersect( const Particle& P );
	
	// Return particle collision distance
	double collision_distance( const double E );

	// Let the Material take care of the collision sample and reaction process
	void collision( Particle& P, std::stack< Particle >& Pbank, const bool ksearch, SourceBank& Fbank, const double k );

	// Simulate scattering for scattering matrix MGXS
	void simulate_scatter( Particle& P );	
        
        // Attached estimators
        std::vector<std::shared_ptr<Estimator>> estimators_C;
        std::vector<std::shared_ptr<Estimator>> estimators_TL;
        // Attaching estimators
	void attach_estimator_C( const std::shared_ptr<Estimator>& E ){ estimators_C.push_back( E ); }
	void attach_estimator_TL( const std::shared_ptr<Estimator>& E ){ estimators_TL.push_back( E ); }
};


#endif // GEOMETRY_H
