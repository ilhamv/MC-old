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


// Scatter particle with scattering angle mu0, with nucleus having mass A
// Scattering angle mu0 is sampled in and passed by the Reaction (see Reaction.h)
// Scattering is trated in Center of mass (COM) frame
// Current model: Free gas scattering with constant cross section
void Particle_t::scatter( const double mu0, const double A )
{
	///////////////////////////////
	// Sampling nuclide velocity //
	
	double V_tilda;  // Nuclide speed candidate
	double mu_tilda; // Nuclide-neutron polar cosine candidate

	const double beta = std::sqrt( 2.0659834e-11 * A ); // Eq. 19
	// note, the constant above is (1.674927471e-27 kg) / (1.38064852e-19 cm^2 kg s^-2 K^-1) / (293.6 K) / 2
	// 	293.6 comes from JANIS data temperature
	
	const double y = beta * p_speed; // Eq. 32
	
	// Sample candidate V_tilda and mu_tilda?
	do
	{
		double x;
		if ( Urand() < 2.0 / ( 2.0 + PI_sqrt * y ) ) // Eq. 37
		{
			// w2 -> sample g2
			x = std::sqrt( -std::log( Urand()*Urand() ) ); // Eq. 39
		}
		else
		{
			// w1 --> sample g1
			const double cos_arg = PI_half * Urand();
			const double cos_val = std::cos( cos_arg );
		 	x = std::sqrt( -std::log( Urand() ) - std::log( Urand() ) * cos_val*cos_val ); // Eq. 38
		}

		V_tilda  = x / beta; // Eq. 32
		mu_tilda = 2.0*Urand() -1.0; // Eq. 40
	}
	
	// Accept candidate V_tilda and mu_tilda?
	while ( Urand() > std::sqrt( p_speed*p_speed + V_tilda*V_tilda - 2.0 * p_speed * V_tilda * mu_tilda ) / ( p_speed + V_tilda ) ); // Eq. 41

	
	Point_t nuclide_dir = scatter_direction( p_dir, mu_tilda );                            // Nuclide direction, Eq. 42
	Point_t V_lab ( nuclide_dir.x*V_tilda, nuclide_dir.y*V_tilda, nuclide_dir.z*V_tilda ); // Nuclide velocity - LAB, Eq. 43

	// Sampling nuclide velocity done //
	////////////////////////////////////


	// Particle velocity - LAB
	Point_t v_lab( p_speed * p_dir.x, p_speed * p_dir.y, p_speed * p_dir.z );
	
	// COM velocity
	const Point_t u ( ( v_lab.x + A*V_lab.x ) / ( 1.0 + A ), ( v_lab.y + A*V_lab.y ) / ( 1.0 + A ), ( v_lab.z + A*V_lab.z ) / ( 1.0 + A ) ); // Eq. 6
	
	// Particle velocity - COM
	Point_t v_c( v_lab.x - u.x, v_lab.y - u.y, v_lab.z - u.z );
	
	// Particle speed - COM
	const double speed_c = std::sqrt( v_c.x*v_c.x+ v_c.y*v_c.y+ v_c.z*v_c.z );

	// Particle initial direction - COM
	const Point_t dir_c( v_c.x / speed_c, v_c.y / speed_c, v_c.z / speed_c );
	
	// Scattering the direction in COM
	Point_t dir_cNew = scatter_direction( dir_c, mu0 ); // Final direction - COM

	// Final velocity - COM
	v_c.x = speed_c * dir_cNew.x;
	v_c.y = speed_c * dir_cNew.y;
	v_c.z = speed_c * dir_cNew.z;
	
	// Final velocity - LAB
	v_lab.x = v_c.x + u.x;
	v_lab.y = v_c.y + u.y;
	v_lab.z = v_c.z + u.z;

	// Final speed - LAB
	setSpeed( std::sqrt( v_lab.x*v_lab.x+ v_lab.y*v_lab.y+ v_lab.z*v_lab.z ) ); // Final energy is computed as well

	// Final direction - LAB
	p_dir.x = v_lab.x / p_speed;
	p_dir.y = v_lab.y / p_speed;
	p_dir.z = v_lab.z / p_speed;
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
