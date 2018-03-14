#include <stack> 

#include "Random.h"
#include "Constants.h"
#include "Algorithm.h"
#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"


//=============================================================================
// Basic Reaction
//=============================================================================

double Reaction::xs( const unsigned long long idx, const double E, 
           const std::vector<double>& E_vec ) 
{ return r_xs->xs( idx, E, E_vec ); }


//=============================================================================
// Scatter Reaction
//=============================================================================

// Scatter with sampled scattering angle mu0, with nucleus mass A
// Scattering is trated in Center of mass (COM) frame
// Current model: Free gas scattering with constant cross section
void ReactionScatter::sample( Particle& P )
{
    const double mu0 = scatter_dist->sample();
	
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Sampling nuclide velocity
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
    double V_tilda;  // Nuclide speed candidate
    double mu_tilda; // Nuclide-neutron polar cosine candidate

    const double beta = std::sqrt( 2.0659834e-11 * A ); // Eq. 19
    // note, the constant above is 
    // (1.674927471e-27 kg) / (1.38064852e-19 cm^2 kg s^-2 K^-1) / (293.6 K)/2
    // 	 293.6 comes from room temperature
	
    const double y = beta * P.speed(); // Eq. 32
	
    // Sample candidate V_tilda and mu_tilda
    do{
	double x;
	if ( Urand() < 2.0 / ( 2.0 + PI_sqrt * y ) ){ // Eq. 37
	    // w2 -> sample g2
	    x = std::sqrt( -std::log( Urand()*Urand() ) ); // Eq. 39
	} else{
	    // w1 --> sample g1
	    const double cos_arg = PI_half * Urand();
	    const double cos_val = std::cos( cos_arg );
	    x = std::sqrt( -std::log( Urand() ) 
                           - std::log( Urand() ) * cos_val*cos_val ); // Eq. 38
	}
	V_tilda  = x / beta; // Eq. 32
	mu_tilda = 2.0*Urand() -1.0; // Eq. 40
    }	
    // Accept candidate V_tilda and mu_tilda?
    while ( Urand() > std::sqrt( P.speed()*P.speed() + V_tilda*V_tilda 
                                 - 2.0 * P.speed() * V_tilda * mu_tilda ) 
                      / ( P.speed() + V_tilda ) ); // Eq. 41
    // Nuclide direction, Eq. 42
    Point nuclide_dir = scatter_direction( P.dir(), mu_tilda ); 
    // Nuclide velocity - LAB, Eq. 43
    Point V_lab ( nuclide_dir.x*V_tilda, nuclide_dir.y*V_tilda, 
                  nuclide_dir.z*V_tilda ); 


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // COM Kinematics
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Particle velocity - LAB
    Point v_lab( P.speed() * P.dir().x, P.speed() * P.dir().y, 
                 P.speed() * P.dir().z );
		
    // COM velocity	
    const Point u ( ( v_lab.x + A*V_lab.x ) / ( 1.0 + A ), 
                    ( v_lab.y + A*V_lab.y ) / ( 1.0 + A ), 
                    ( v_lab.z + A*V_lab.z ) / ( 1.0 + A ) ); // Eq. 6
	
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
	
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Convert to LAB
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Final velocity - LAB
    v_lab.x = v_c.x + u.x;
    v_lab.y = v_c.y + u.y;
    v_lab.z = v_c.z + u.z;

    // Final speed - LAB
    P.set_speed( std::sqrt( v_lab.x*v_lab.x + v_lab.y*v_lab.y 
                            + v_lab.z*v_lab.z ) ); // Final energy is computed

    // Final direction - LAB
    P.set_direction( Point( v_lab.x / P.speed(), v_lab.y / P.speed(),
                            v_lab.z / P.speed() ) );
}


//=============================================================================
// Fission Reaction
//=============================================================================

double ReactionFission::nu( const unsigned long long idx, const double E, 
           const std::vector<double>& E_vec ) 
{ return r_nu->xs( idx, E, E_vec ); }

double ReactionFission::Chi( const double E ) 
{ return r_Chi->sample( E ); }
double ReactionFission::ChiD( const int g, const double E ) 
{ return r_ChiD[g]->sample( E ); }
double ReactionFission::beta( const unsigned long long idx, const double E, 
          const std::vector<double>& E_vec )
{ return r_beta->xs( idx, E, E_vec ); }
double ReactionFission::lambda( const int g ) 
{ return r_lambda[g]; }
double ReactionFission::fraction( const int g ) 
{ return r_fraction[g]; }
double ReactionFission::f_lambda( const int g ) 
{ return r_f_lambda[g]; }
