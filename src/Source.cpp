#include <cmath> // sqrt

#include "Distribution.h"
#include "Random.h"
#include "Particle.h"
#include "Source.h"
#include "Point.h"


Particle_t Point_Source::getSource()
{
    Point p = dist_dir->sample();
	
    Particle_t P( pos, p, dist_enrg->sample(), dist_time->sample(), 1.0, 0 );
    return P;
}


Particle_t Delta_Source::getSource()
{
    Particle_t P( pos, dir, E, t, w, 0 );
    return P;
}


Particle_t DiskX_Source::getSource()
{
    Point p = dist_dir->sample();

    double z,y;
    // Rejection sampling: square --> circle
    do
    {
	z = 2.0 * Urand() - 1.0;
	y = 2.0 * Urand() - 1.0;
    }
    while ( z*z + y*y > 1.0 );

    z = z * r;
    y = y * r;
    Point pos( x0, y0 + y, z0 + z );
	
    Particle_t P( pos, p, dist_enrg->sample(), dist_time->sample(), 1.0, 0 );
    return P;
}


Particle_t DiskZ_Source::getSource()
{
    Point p = dist_dir->sample();

    double x,y;
    // Rejection sampling: square --> circle
    do
    {
	x = 2.0 * Urand() - 1.0;
	y = 2.0 * Urand() - 1.0;
    }
    while ( x*x + y*y > 1.0 );

    x = x * r;
    y = y * r;
    Point pos( x0 + x, y0 + y, z0 );
	
    Particle_t P( pos, p, dist_enrg->sample(), dist_time->sample(), 1.0, 0 );
    return P;
}


Particle_t Sphere_Shell_Source::getSource()
{
    Point p = dist_dir->sample();

    double x,y,z;
    // Rejection sampling: Cube --> Spherical shell
    do
    {
	x = 2.0 * Urand() - 1.0;
    	y = 2.0 * Urand() - 1.0;
    	z = 2.0 * Urand() - 1.0;
    }
    while ( ( x*x + y*y + z*z > 1.0 ) || ( x*x + y*y + z*z < risq ) );

    x = x * ro;
    y = y * ro;
    z = z * ro;

    Point pos( x0 + x, y0 + y, z + z0 );

    Particle_t P( pos, p, dist_enrg->sample(), dist_time->sample(), 1.0, 0 );
    return P;
}


Particle_t Generic_Source::getSource() 
{
    Point p = dist_dir->sample();
	
    Particle_t P ( dist_pos->sample(), p, dist_enrg->sample(), dist_time->sample(), 1.0, 0 );
    return P;
}


// Get source
// sources are sampled wrt to their probability
// then, particle cell is searched and set
Particle_t Source_Bank::getSource( const std::vector<std::shared_ptr<Cell>>& Cell )
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
    std::cout<< "[ERROR] There is no sample found\n";
    std::exit(EXIT_FAILURE);
}
