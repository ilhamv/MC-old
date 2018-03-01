#include <cmath>      // cos, sin, sqrt

#include "Geometry.h"
#include "Particle.h"
#include "Random.h"
#include "Point.h"
#include "Constants.h"    // PI2, PI_sqrt, PI_half
#include "Algorithm.h"   // scatter_direction

//=============================================================================
// Getters
//=============================================================================

Point Particle::pos() const { return p_pos; }
Point Particle::dir() const { return p_dir; }
bool Particle::alive() const { return p_alive; }
double Particle::weight() const { return p_weight; }
double Particle::time() const { return p_time; }  
double Particle::time_old() const { return p_time_old; }   
double Particle::energy() const { return p_energy; }
double Particle::energy_old() const { return p_energy_old; }
double Particle::speed() const { return p_speed; } 
std::shared_ptr<Cell> Particle::cell() const { return p_cell; }
std::shared_ptr<Cell> Particle::cell_old() const { return p_cell_old; }  
std::shared_ptr<Surface_t> Particle::surface_old() const 
{ 
    return p_surface_old; 
}
int Particle::tdmc() const { return p_tdmc; }

//=============================================================================
// Setters
//=============================================================================

void Particle::set_direction( const Point& p ) { p_dir = p; }
void Particle::set_weight( const double w ) { p_weight = w; }
void Particle::set_cell( const std::shared_ptr<Cell> C )
{ 
    p_cell_old = p_cell;
    p_cell     = C;
}
void Particle::set_energy( const double E )
{ 
    p_energy_old = p_energy;
    p_energy     = E;
    p_speed  = std::sqrt( p_energy * 191312955.067 ) * 100.0;
    // constant: 2.0 * ( 1.60217662e-19 J/eV ) / ( 1.674927471e-27 kg )
}
void Particle::set_speed( const double v )
{ 
    p_speed  = v;
    p_energy_old = p_energy;
    p_energy = 5.2270376e-13 * v * v;
    // constant: 0.5 / ( 1.60217662e-19 J/eV ) * ( 1.674927471e-27 kg ) 
    //           / ( 10000 cm^2/m^2 )
}
void Particle::set_surface_old( const std::shared_ptr<Surface_t> S )
{
    p_surface_old = S;
}
                               
//=============================================================================
// Modifiers
//=============================================================================

void Particle::move( const double dmove ) 
{
    // Move particle a distance dmove
    p_pos.x += p_dir.x * dmove;
    p_pos.y += p_dir.y * dmove;
    p_pos.z += p_dir.z * dmove;
	
    // Advance particle time
    p_time_old = p_time;
    p_time += dmove / p_speed;
}
void Particle::kill()
{ 
    p_alive  = false; 
    p_weight = 0.0;
}
void Particle::increment_tdmc() { p_tdmc++; }
