#include <cmath>      // cos, sin, sqrt

#include "Geometry.h"
#include "Particle.h"
#include "Random.h"
#include "Point.h"
#include "Const.h"    // PI2, PI_sqrt, PI_half
#include "Solver.h"   // scatter_direction


// Getters
Point_t                 Particle_t::pos()    const { return p_pos; }    // Position
Point_t                 Particle_t::dir()    const { return p_dir; }    // Direction
bool                    Particle_t::alive()  const { return p_alive; }  // Alive/dead flag
double                  Particle_t::weight() const { return p_weight; } // Weight
double                  Particle_t::time()   const { return p_time; }   // Time
double                  Particle_t::energy() const { return p_energy; } // Energy
double                  Particle_t::speed()  const { return p_speed; }  // Speed
std::shared_ptr<Cell_t> Particle_t::cell()   const { return p_cell; }   // Particle cell


// Setters
void Particle_t::setDirection( const Point_t& p )            { p_dir    = p; } // Direction
void Particle_t::setWeight( const double w )                 { p_weight = w; } // Weight
void Particle_t::setCell( const std::shared_ptr<Cell_t>& C ) { p_cell = C; }   // Particle cell
void Particle_t::setTime( const double t )                   { p_time   = t; } // Elapsed time
void Particle_t::setEnergy( const double E )                                       // Energy, and speed
{ 
	p_energy = E; // eV
	p_speed  = std::sqrt( p_energy * 191312955.067 ) * 100.0; // cm/s
	// note, the constant above is 2.0 * ( 1.60217662e-19 J/eV ) / ( 1.674927471e-27 kg )
}
void Particle_t::setSpeed( const double v )                                        // Speed, and energy
{ 
	p_speed  = v; // cm/s
	p_energy = 5.2270376e-13 * v * v; // eV
	// note, the constant above is 0.5 / ( 1.60217662e-19 J/eV ) * ( 1.674927471e-27 kg ) / ( 10000 cm^2/m^2 )
}


// Move particle a distance dmove along its current trajectory
void Particle_t::move( const double dmove ) 
{
	p_pos.x += p_dir.x * dmove;
	p_pos.y += p_dir.y * dmove;
	p_pos.z += p_dir.z * dmove;
	
	// Advance particle time
	setTime( p_time + dmove / p_speed );
}


// Kill particle
void Particle_t::kill()
{ 
	p_alive  = false; 
	p_weight = 0.0;
}


// Search and set particle cell
void Particle_t::searchCell( const std::vector<std::shared_ptr<Cell_t>>& Cell )
{
    for( const auto& C : Cell )
	{
		// check if particle is in the current cell C
		if ( C->testPoint( p_pos ) )
		{
			p_cell = C;
			return;
		}
	}
	std::cout<< "[WARNING] A particle is lost:\n( x, y, z )  (" << p_pos.x << ", " << p_pos.y << ", " << p_pos.z << " )\n";
        std::exit(EXIT_FAILURE);
}
