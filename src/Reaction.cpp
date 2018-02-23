#include <stack> // stack

#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"


// Capture kills the working particle
void Capture_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{ P.kill(); }


// Scatter the particle with sampled scattering angle mu0, with nucleus having mass A
// Scattering is trated in Center of mass (COM) frame
// Current model: Free gas scattering with constant cross section
void Scatter_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{
    const double mu0 = scatter_dist->sample();
	
    ///////////////////////////////
    // Sampling nuclide velocity //
	
    double V_tilda;  // Nuclide speed candidate
    double mu_tilda; // Nuclide-neutron polar cosine candidate

    const double beta = std::sqrt( 2.0659834e-11 * A ); // Eq. 19
    // note, the constant above is (1.674927471e-27 kg) / (1.38064852e-19 cm^2 kg s^-2 K^-1) / (293.6 K) / 2
    // 	 293.6 comes from room temperature JANIS data
	
    const double y = beta * P.speed(); // Eq. 32
	
    // Sample candidate V_tilda and mu_tilda
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
    while ( Urand() > std::sqrt( P.speed()*P.speed() + V_tilda*V_tilda - 2.0 * P.speed() * V_tilda * mu_tilda ) / ( P.speed() + V_tilda ) ); // Eq. 41

	
    Point nuclide_dir = scatter_direction( P.dir(), mu_tilda );                            // Nuclide direction, Eq. 42
    Point V_lab ( nuclide_dir.x*V_tilda, nuclide_dir.y*V_tilda, nuclide_dir.z*V_tilda ); // Nuclide velocity - LAB, Eq. 43

    // Sampling nuclide velocity done //
    ////////////////////////////////////


    // Particle velocity - LAB
    Point v_lab( P.speed() * P.dir().x, P.speed() * P.dir().y, P.speed() * P.dir().z );
		
    // COM velocity	
    const Point u ( ( v_lab.x + A*V_lab.x ) / ( 1.0 + A ), ( v_lab.y + A*V_lab.y ) / ( 1.0 + A ), ( v_lab.z + A*V_lab.z ) / ( 1.0 + A ) ); // Eq. 6
	
    // Particle velocity - COM
    Point v_c( v_lab.x - u.x, v_lab.y - u.y, v_lab.z - u.z );
	
    // Particle speed - COM
    const double speed_c = std::sqrt( v_c.x*v_c.x+ v_c.y*v_c.y+ v_c.z*v_c.z );

    // Particle initial direction - COM
    const Point dir_c( v_c.x / speed_c, v_c.y / speed_c, v_c.z / speed_c );
	
    // Scattering the direction in COM
    Point dir_cNew = scatter_direction( dir_c, mu0 ); // Final direction - COM

    // Final velocity - COM
    v_c.x = speed_c * dir_cNew.x;
    v_c.y = speed_c * dir_cNew.y;
    v_c.z = speed_c * dir_cNew.z;
	
    // Final velocity - LAB
    v_lab.x = v_c.x + u.x;
    v_lab.y = v_c.y + u.y;
    v_lab.z = v_c.z + u.z;

    // Final speed - LAB
    P.setSpeed( std::sqrt( v_lab.x*v_lab.x+ v_lab.y*v_lab.y+ v_lab.z*v_lab.z ) ); // Final energy is computed as well

    // Final direction - LAB
    P.setDirection( Point( v_lab.x / P.speed(), v_lab.y / P.speed(), v_lab.z / P.speed() ) );
}


// Scatter the particle with sampled scattering angle mu0, with nucleus having mass A
// Scattering is trated in Center of mass (COM) frame
// Current model: Free gas scattering with constant cross section
void Scatter_Zero_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{
    const double mu0 = scatter_dist->sample();

    // Particle velocity - LAB
    Point v_lab( P.speed() * P.dir().x, P.speed() * P.dir().y, P.speed() * P.dir().z );
		
    // COM velocity	
    const Point u ( v_lab.x / ( 1.0 + A ), v_lab.y / ( 1.0 + A ), v_lab.z / ( 1.0 + A ) ); // Eq. 6
	
    // Particle velocity - COM
    Point v_c( v_lab.x - u.x, v_lab.y - u.y, v_lab.z - u.z );
	
    // Particle speed - COM
    const double speed_c = std::sqrt( v_c.x*v_c.x+ v_c.y*v_c.y+ v_c.z*v_c.z );

    // Particle initial direction - COM
    const Point dir_c( v_c.x / speed_c, v_c.y / speed_c, v_c.z / speed_c );
	
    // Scattering the direction in COM
    Point dir_cNew = scatter_direction( dir_c, mu0 ); // Final direction - COM

    // Final velocity - COM
    v_c.x = speed_c * dir_cNew.x;
    v_c.y = speed_c * dir_cNew.y;
    v_c.z = speed_c * dir_cNew.z;
	
    // Final velocity - LAB
    v_lab.x = v_c.x + u.x;
    v_lab.y = v_c.y + u.y;
    v_lab.z = v_c.z + u.z;

    // Final speed - LAB
    P.setSpeed( std::sqrt( v_lab.x*v_lab.x+ v_lab.y*v_lab.y+ v_lab.z*v_lab.z ) ); // Final energy is computed as well

    // Final direction - LAB
    P.setDirection( Point( v_lab.x / P.speed(), v_lab.y / P.speed(), v_lab.z / P.speed() ) );
}


// Fission reaction sample
void Fission_Reaction::sample( Particle_t& P, std::stack< Particle_t >& Pbank )
{
	// create random number of secondaries from multiplicity distributon nu_dist and
	// push all but one of them into the Particle bank, and reset the top particle 
	// if no secondaries, kill the particle

	const int n = std::floor( r_nu->xs( P.energy() ) + Urand() ); // sampled multiplicity
    
	if ( n != 0 ) 
    	{
        	// bank all but last particle (skips if n = 1)
	        for ( int i = 0 ; i < n - 1 ; i++ )
        	{
	            	Particle_t p( P.pos(), isotropic.sample(), Chi_dist->sample( P.energy() ), P.time(), P.weight(), P.tdmc() );
        	    	p.setCell( P.cell() );
            		Pbank.push( p );
        	}

		// reset the top particle
		P.setDirection( isotropic.sample() );
		P.setEnergy( Chi_dist->sample( P.energy() ) );
	}
	else
	{ P.kill(); }
}
