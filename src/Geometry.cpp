#include "Geometry.h"
#include "Particle.h"
#include "Point.h"
#include "Constants.h" 
#include "Algorithm.h" 
#include "Estimator.h"
#include "Material.h"


//=============================================================================
// Surface: Constructors
//=============================================================================

SurfacePlane::SurfacePlane( const std::string n, const int i, const int bc, 
                              const double pa, const double pb, const double pc, 
                              const double pd ):		 
    Surface(n,i,bc), a(pa), b(pb), c(pc), d(pd)
{
    const double L = 2.0 / ( a*a + b*b + c*c );
    modx = L * a;
    mody = L * b;
    modz = L * c;
}

//=============================================================================
// Surface: Evaluation
//=============================================================================

double SurfacePlaneX::eval( const Point& p ) 
{ 
    return p.x - x; 
}
double SurfacePlaneY::eval( const Point& p ) 
{ 
    return p.y - y; 
}
double SurfacePlaneZ::eval( const Point& p ) 
{ 
    return p.z - z; 
}
double SurfacePlane::eval( const Point& p )
{
    return a * p.x + b * p.y + c * p.z - d;
}
double SurfaceSphere::eval( const Point& p ) 
{
    const double x_t = p.x - x0;
    const double y_t = p.y - y0;
    const double z_t = p.z - z0;
    return x_t*x_t + y_t*y_t + z_t*z_t - rad_sq;
}
double SurfaceCylinderX::eval( const Point& p )
{
    const double y_t = p.y - y0;
    const double z_t = p.z - z0;
    return y_t*y_t + z_t*z_t - rad_sq;
}
double SurfaceCylinderY::eval( const Point& p )
{
    const double x_t = p.x - x0;
    const double z_t = p.y - z0;
    return x_t*x_t + z_t*z_t - rad_sq;
}
double SurfaceCylinderZ::eval( const Point& p )
{
    const double x_t = p.x - x0;
    const double y_t = p.y - y0;
    return x_t*x_t + y_t*y_t - rad_sq;
}


//=============================================================================
// Surface: Distance to Hit
//=============================================================================

double SurfacePlaneX::distance( const Particle& P )
{
    const double pos = P.pos().x;
    const double dir = P.dir().x;

    // Check if particle moves in a direction that is (or very close to) 
    //   parallel to the surface
    if ( std::fabs( dir ) > EPSILON_float ) {
    	const double dist = ( x - pos ) / dir;
    	// Check if particle moves away from the surface
	if ( dist > 0.0 ) { return dist; }
    	else { return MAX_float; }
    }    	
    // It does! (Parallel)
    else 
    { return MAX_float; }
}
double SurfacePlaneY::distance( const Particle& P )
{
    const double pos = P.pos().y;
    const double dir = P.dir().y;

    // Check if particle moves in a direction that is (or very close to) 
    //   parallel to the surface
    if ( std::fabs( dir ) > EPSILON_float ) {
    	const double dist = ( y - pos ) / dir;
    	// Check if particle moves away from the surface
	if ( dist > 0.0 ) { return dist; }
    	else { return MAX_float; }
    }    	
    // It does! (Parallel)
    else 
    { return MAX_float; }
}
double SurfacePlaneZ::distance( const Particle& P )
{
    const double pos = P.pos().z;
    const double dir = P.dir().z;

    // Check if particle moves in a direction that is (or very close to)
    //   parallel to the surface
    if ( std::fabs( dir ) > EPSILON_float ) {
    	const double dist = ( z - pos ) / dir;
    	// Check if particle moves away from the surface
	if ( dist > 0.0 ) { return dist; }
    	else { return MAX_float; }
    }    	
    // It does! (Parallel)
    else 
    { return MAX_float; }
}
double SurfacePlane::distance( const Particle& P )
{
    Point pos = P.pos();
    Point dir = P.dir();
    const double denom = a * dir.x  +  b * dir.y  +  c * dir.z;

    // Check if particle moves in a direction that is (or very close to)
    //   parallel to the surface
    if ( std::fabs( denom ) > EPSILON_float ){
    	const double dist = ( d - a * pos.x - b * pos.y - c * pos.z ) / denom;
    	// Check if particle moves away from the surface
	if ( dist > 0.0 ) { return dist; }
    	else { return MAX_float; }
    }    	
    // It does! (Parallel)
    else 
    { return MAX_float; }
}
double SurfaceSphere::distance( const Particle& P ) 
{
    Point p = P.pos();
    Point u = P.dir();

    // put into quadratic equation form: a*s^2 + b*s + c = 0, where a = 1
    double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.y - y0 ) * u.y + ( p.z - z0 ) * u.z );
    double c = eval( p );

    return geometry_quad( 1.0, b, c );
}
double SurfaceCylinderX::distance( const Particle& P )
{
    Point p = P.pos();
    Point u = P.dir();

    double a = 1.0 - u.x*u.x;
    double b = 2.0 * ( ( p.y - y0 ) * u.y + ( p.z - z0 ) * u.z );
    double c = eval( p );

    return geometry_quad( a, b, c );
}
double SurfaceCylinderY::distance( const Particle& P )
{
    Point p = P.pos();
    Point u = P.dir();

    double a = 1.0 - u.y*u.y;
    double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.z - z0 ) * u.z );
    double c = eval( p );

    return geometry_quad( a, b, c );
}
double SurfaceCylinderZ::distance( const Particle& P )
{
  	Point p = P.pos();
  	Point u = P.dir();

  	double a = 1.0 - u.z*u.z;
	double b = 2.0 * ( ( p.x - x0 ) * u.x + ( p.y - y0 ) * u.y );
  	double c = eval( p );

  	return geometry_quad( a, b, c );
}


//=============================================================================
// Surface: Reflect
//=============================================================================

void SurfacePlaneX::reflect( Particle& P )
{
    Point q( -P.dir().x, P.dir().y, P.dir().z );
    P.set_direction(q);
}
void SurfacePlaneY::reflect( Particle& P )
{
    Point q( P.dir().x, -P.dir().y, P.dir().z );
    P.set_direction(q);
}
void SurfacePlaneZ::reflect( Particle& P )
{
    Point q( P.dir().x, P.dir().y, -P.dir().z );
    P.set_direction(q);
}
void SurfacePlane::reflect( Particle& P )
{
    const double K = ( a * P.dir().x + b * P.dir().y + c * P.dir().z );
    Point q;
    q.x = P.dir().x - modx * K;
    q.y = P.dir().y - mody * K;
    q.z = P.dir().z - modz * K;
    P.set_direction(q);
}
void SurfaceSphere::reflect( Particle& P ) { return; }
void SurfaceCylinderX::reflect( Particle& P ) { return; }
void SurfaceCylinderY::reflect( Particle& P ) { return; }
void SurfaceCylinderZ::reflect( Particle& P ) { return; }


//=============================================================================
// Cell
//=============================================================================

double Cell::importance() { return c_importance; } 
std::shared_ptr<Material> Cell::material() { return c_material; }
std::vector<std::pair<std::shared_ptr<Surface>,int>>& Cell::surfaces()
{ return c_surfaces; }

void Cell::add_surface( const std::shared_ptr< Surface >& S, const int sense )
{ 
    c_surfaces.push_back( std::make_pair( S, sense ) ); 
}
void Cell::set_material( const std::shared_ptr< Material >& M ) 
{ 
    c_material = M; 
}
bool Cell::test_point( const Point& p )
{
    // Loop over surfaces in cell, if not on correct side return false
    for ( const auto& S : c_surfaces ) {
    	if ( S.first->eval( p ) * S.second < 0 ) { return false; }  
    }
    return true;
}
double Cell::collision_distance( const double E )
{ 
    if ( c_material ) { return c_material->collision_distance_sample( E ); }
    // Vacuum --> return sligthly less than very large number
    //            to ensure collision (kill) if no surface intersection
    else { return MAX_float_less; }
}
std::pair<std::shared_ptr<Surface>, double> 
Cell::surface_intersect( const Particle& P ) 
{
    double dist = MAX_float;
    std::shared_ptr< Surface > S = nullptr;
    for ( const auto& s : c_surfaces ) {
    	double d = s.first->distance( P );
    	if ( d < dist ){ 
	    dist = d;
	    S    = s.first;
	}
    }
    return std::make_pair( S, dist ); 
}
void Cell::collision( Particle& P, std::stack< Particle >& Pbank, 
                      const bool ksearch, SourceBank& Fbank, const double k )
{ 
    if ( c_material ){ 
        c_material->collision_sample( P, Pbank, ksearch, Fbank, k ); 
    }
    // Vacuum --> Kill particle at collision
    else { return P.kill(); }
}	
void Cell::simulate_scatter( Particle& P )
{ c_material->simulate_scatter( P ); }
