#include "Geometry.h"
#include "Particle.h"
#include "Point.h"
#include "Const.h"  // MAX, EPSILON
#include "Solver.h" // solve_quad



////////////////
/// Surfaces ///
///////////////

// Crossing the surface --> an epsilon kick to the working particle
void Surface_t::cross ( Particle_t& P ) { P.move( EPSILON ); }

// Hit implementation
void Surface_t::hit( Particle_t& P, const std::vector<std::shared_ptr<Cell_t>>& Cell, const bool tally )
{
	// Note: new particle cell search is only performed in transmission
	// Transmission
	if ( bc == "transmission" )
	{
		// Cross the surface (Epsilon kick)
		cross( P );
		// Search and set new cell
		P.searchCell( Cell );
	}
	// Reflective
	else
	{
		// Reflect angle
		reflect( P );
		// Cross the surface (Epsilon kick)
		cross ( P );
	}
	
	// Score estimators 
        if (tally)
        {
	    const double told = P.time(); // crossing surface happends instantly
	    for ( auto& e : estimators ) { e->score( P, told ); }
        }
}


// Plane-X
double PlaneX_Surface::eval( const Point_t& p ) { return p.x - x; }

double PlaneX_Surface::distance( const Particle_t& P )
{
	const double pos = P.pos().x;
	const double dir = P.dir().x;

	// Check if particle moves in a direction that is (or very close to) parallel to the surface
  	if ( std::fabs( dir ) > EPSILON ) 
	{
    		const double dist = ( x - pos ) / dir;
    		// Check if particle moves away from the surface
		if ( dist > 0.0 ) { return dist; }
    		else { return MAX; }
  	}    	
	// It does! (Parallel)
  	else 
	{ return MAX; }
}

void PlaneX_Surface::reflect( Particle_t& P )
{
	Point_t q( -P.dir().x, P.dir().y, P.dir().z );
	P.setDirection(q);
}


// Plane-Y
double PlaneY_Surface::eval( const Point_t& p ) { return p.y - y; }

double PlaneY_Surface::distance( const Particle_t& P )
{
	const double pos = P.pos().y;
	const double dir = P.dir().y;

	// Check if particle moves in a direction that is (or very close to) parallel to the surface
  	if ( std::fabs( dir ) > EPSILON ) 
	{
    		const double dist = ( y - pos ) / dir;
    		// Check if particle moves away from the surface
		if ( dist > 0.0 ) { return dist; }
    		else { return MAX; }
  	}    	
	// It does! (Parallel)
  	else 
	{ return MAX; }
}

void PlaneY_Surface::reflect( Particle_t& P )
{
	Point_t q( P.dir().x, -P.dir().y, P.dir().z );
	P.setDirection(q);
}


// Plane-Z
double PlaneZ_Surface::eval( const Point_t& p ) { return p.z - z; }

double PlaneZ_Surface::distance( const Particle_t& P )
{
	const double pos = P.pos().z;
	const double dir = P.dir().z;

	// Check if particle moves in a direction that is (or very close to) parallel to the surface
  	if ( std::fabs( dir ) > EPSILON ) 
	{
    		const double dist = ( z - pos ) / dir;
    		// Check if particle moves away from the surface
		if ( dist > 0.0 ) { return dist; }
    		else { return MAX; }
  	}    	
	// It does! (Parallel)
  	else 
	{ return MAX; }
}

void PlaneZ_Surface::reflect( Particle_t& P )
{
	Point_t q( P.dir().x, P.dir().y, -P.dir().z );
	P.setDirection(q);
}


// Generic Plane
double Plane_Surface::eval( const Point_t& p )
{
	return a * p.x  +  b * p.y  +  c * p.z  - d;
}

double Plane_Surface::distance( const Particle_t& P )
{
	Point_t pos = P.pos();
	Point_t dir = P.dir();

	const double denom = a * dir.x  +  b * dir.y  +  c * dir.z;

	// Check if particle moves in a direction that is (or very close to) parallel to the surface
  	if ( std::fabs( denom ) > EPSILON ) 
	{
    		const double dist = ( d - a * pos.x - b * pos.y - c * pos.z ) / denom;
    		// Check if particle moves away from the surface
		if ( dist > 0.0 ) { return dist; }
    		else { return MAX; }
  	}    	
	// It does! (Parallel)
  	else 
	{ return MAX; }
}

void Plane_Surface::reflect( Particle_t& P )
{
	const double K = ( a * P.dir().x + b * P.dir().y + c * P.dir().z );
	Point_t q;
	
	q.x = P.dir().x - modx * K;
	q.y = P.dir().y - mody * K;
	q.z = P.dir().z - modz * K;

	P.setDirection(q);
}


// Sphere
double Sphere_Surface::eval( const Point_t& p ) 
{
  	const double x_t = p.x - x0;
  	const double y_t = p.y - y0;
  	const double z_t = p.z - z0;
	return x_t*x_t + y_t*y_t + z_t*z_t - rad_sq;
}

double Sphere_Surface::distance( const Particle_t& P ) 
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	// put into quadratic equation form: a*s^2 + b*s + c = 0, where a = 1
  	double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.y - y0 ) * u.y + ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return solve_quad( 1.0, b, c );
}

void Sphere_Surface::reflect( Particle_t& P ) { return; }


// Cylinder-X
double CylinderX_Surface::eval( const Point_t& p )
{
	const double y_t = p.y - y0;
	const double z_t = p.z - z0;
	return y_t*y_t + z_t*z_t - rad_sq;
}

double CylinderX_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	double a = 1.0 - u.x*u.x;
	double b = 2.0 * ( ( p.y - y0 ) * u.y + ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}

void CylinderX_Surface::reflect( Particle_t& P ) { return; }


// Cylinder-Y
double CylinderY_Surface::eval( const Point_t& p )
{
	const double x_t = p.x - x0;
	const double z_t = p.y - z0;
	return x_t*x_t + z_t*z_t - rad_sq;
}

double CylinderY_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	double a = 1.0 - u.y*u.y;
	double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}

void CylinderY_Surface::reflect( Particle_t& P ) { return; }


// Cylinder-Z
double CylinderZ_Surface::eval( const Point_t& p )
{
	const double x_t = p.x - x0;
	const double y_t = p.y - y0;
	return x_t*x_t + y_t*y_t - rad_sq;
}

double CylinderZ_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	double a = 1.0 - u.z*u.z;
	double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.y - y0 ) * u.y );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}

void CylinderZ_Surface::reflect( Particle_t& P ) { return; }


// Cone-X
double ConeX_Surface::eval( const Point_t& p )
{
	const double x_t = p.x - x0;
	const double y_t = p.y - y0;
	const double z_t = p.z - z0;
	return - rad_sq * x_t*x_t + y_t*y_t + z_t*z_t;
}

double ConeX_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	double a = 1.0 - ( rad_sq + 1.0 ) * u.x*u.x;
  	double b = 2.0 * ( - rad_sq * ( p.x - x0 ) * u.x + ( p.y - y0 ) * u.y + ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}

void ConeX_Surface::reflect( Particle_t& P ) { return; }


// Cone-Y
double ConeY_Surface::eval( const Point_t& p )
{
	const double x_t = p.x - x0;
	const double y_t = p.y - y0;
	const double z_t = p.z - z0;
	return x_t*x_t - rad_sq * y_t*y_t + z_t*z_t;
}

double ConeY_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	double a = 1.0 - ( rad_sq + 1.0 ) * u.y*u.y;
  	double b = 2.0 * ( ( p.x - x0 ) * u.x - rad_sq * ( p.y - y0 ) * u.y + ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}

void ConeY_Surface::reflect( Particle_t& P ) { return; }


// Cone-Z
double ConeZ_Surface::eval( const Point_t& p )
{
	const double x_t = p.x - x0;
	const double y_t = p.y - y0;
	const double z_t = p.z - z0;
	return x_t*x_t + y_t*y_t - rad_sq * z_t*z_t;
}

double ConeZ_Surface::distance( const Particle_t& P )
{
  	Point_t p = P.pos();
  	Point_t u = P.dir();

  	double a = 1.0 - ( rad_sq + 1.0 ) * u.z*u.z;
  	double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.y - y0 ) * u.y - rad_sq * ( p.z - z0 ) * u.z );
  	double c = eval( p );

  	return solve_quad( a, b, c );
}

void ConeZ_Surface::reflect( Particle_t& P ) { return; }



//////////////
/// Cell ///
//////////////

// Getters
double      Cell_t::importance() { return r_importance; } // importance
// Get macroXsec of the contained material
double      Cell_t::SigmaT  ( const double E ) { return material->SigmaT( E ); }
double      Cell_t::SigmaA  ( const double E ) { return material->SigmaA( E ); }
double      Cell_t::SigmaS  ( const double E ) { return material->SigmaS( E ); }
double      Cell_t::SigmaC  ( const double E ) { return material->SigmaC( E ); }
double      Cell_t::SigmaF  ( const double E ) { return material->SigmaF( E ); }
double      Cell_t::nuSigmaF( const double E ) { return material->nuSigmaF( E ); }


// Take in a pair of surface pointer and integer describing sense
// and append to vector of surfaces
void Cell_t::addSurface( const std::shared_ptr< Surface_t >& S, const int sense )
{ surfaces.push_back( std::make_pair( S, sense ) ); }


// Add the material
void Cell_t::setMaterial( const std::shared_ptr< Material_t >& M ) 
{ material = M; }


// Test if point is inside the cell
bool Cell_t::testPoint( const Point_t& p )
{
  	// Loop over surfaces in cell, if not on correct side return false
  	// if on correct side of all surfaces, particle is in the cell and return true
  	for ( const auto& S : surfaces ) 
	{
    		// first = surface pointer, second = +/- 1 indicating sense
    		if ( S.first->eval( p ) * S.second < 0 ) { return false; }  
  	}
  	return true;
}


// Move particle and score any estimators
void Cell_t::moveParticle( Particle_t& P, const double dmove, const bool tally )
{
	const double told = P.time();
	P.move( dmove );

        // Score track length estimator
        if (tally)
        { for ( const auto& e : estimators ) { e->score( P, told, dmove ); } }
}


// Find the closest surface and travel distance for particle p to reach
std::pair< std::shared_ptr< Surface_t >, double > Cell_t::surface_intersect( const Particle_t& P ) 
{
  	double dist = MAX;
	std::shared_ptr< Surface_t > S = nullptr;
  	for ( const auto& s : surfaces ) 
	{
    		double d = s.first->distance( P );
    		if ( d < dist )
		{ 
			dist = d;
			S    = s.first;
		}
  	}
	return std::make_pair( S, dist ); 
}


// Return particle collision distance
double Cell_t::collision_distance( const double E )
{ 
	if ( material ) 
	{ return material->collision_distance_sample( E ); }
	// Vacuum --> return sligthly less than very large number for collision distance
	// to ensure collision if no surface intersection
	else { return MAX_less; } // MAX_less = 0.9 MAX
}


// Collision
// Let the Material take care of the collision sample and reaction process
void Cell_t::collision( Particle_t& P, std::stack< Particle_t >& Pbank, const bool ksearch, Source_Bank& Fbank, const double k )
{ 
	if ( material ) 
	{ material->collision_sample( P, Pbank, ksearch, Fbank, k ); }
	// Vacuum --> Kill particle at collision
	else { return P.kill(); }
}	


// Simulate scattering for scattering matrix MGXS
void Cell_t::simulate_scatter( Particle_t& P )
{ material->simulate_scatter( P ); }
