#include <vector>  // vector
#include <iostream>// cout
#include <cstring> // string, strcmp
#include <memory>  // shared_ptr, make_shared
#include <stack>   // stack
#include <cmath>   // exp
#include <sstream> // ostringstream
#include <fstream> // ifstream

#include "Simulator.h"
#include "XMLparser.h"

// Constructor: Set up the simulator with XML parser
Simulator_t::Simulator_t( const std::string input_file )
{
    XML_input( input_file, simName, nhist, Ecut_off, tcut_off, Sbank, Surface, Cell, Nuclide, Material, Estimator, double_distributions, int_distributions, point_distributions );
}

void Simulator_t::start()
{
    // Simulation loop
    for ( unsigned int isample = 0 ; isample < nhist ; isample++ )
    {
        Pbank.push( Sbank.getSource( Cell ) );
		
        // History loop
        while ( !Pbank.empty() )
        {
            Particle_t                P = Pbank.top(); // Working particle
            std::shared_ptr<Cell_t> C = P.cell();  // Working cell
            Pbank.pop();

            // Particle loop
            while ( P.alive() )
            {
                std::pair< std::shared_ptr< Surface_t >, double > SnD; // To hold nearest surface and its distance
				
                // Determine nearest surface and its distance
                SnD = C->surface_intersect( P );

                // Determine collision distance
                double dcol = C->collision_distance( P.energy() );
				
                // Hit surface?
                if ( dcol > SnD.second )
                {	
                    // Surface hit! Move particle to surface, tally if there is any Cell Tally
                    C->moveParticle( P, SnD.second );

                    // Implement surface hit:
                    // 	Reflect angle for reflective surface
                    // 	Cross the surface (move epsilon distance)
                    // 	Search new particle cell for transmission surface
                    // 	Tally if there is any Surface Tally
                    // 	Note: particle weight and working cell are not updated yet
                    SnD.first->hit( P, Cell );

                    // Splitting & Roulette Variance Reduction Technique
                    // 	More important : split
                    // 	Less important : roulette
                    // 	Same importance: do nothing
                    // 	Note: old working cell has the previous cell importance and will be updated
                    Split_Roulette( C, P, Pbank );
                }
				
                // Collide!!
                else
                {
                    // Move particle to collision site and sample the collision and tally if there is any cell tally
                    C->moveParticle( P, dcol );
                    C->collision( P, Pbank );			
                }
			
                // add # of tracks
                tracks++;

                // Cut-off working particle?
                if ( P.energy() < Ecut_off || P.time() > tcut_off || P.weight() == 0.0 ) { P.kill();}

            } // Particle is dead, end of particle loop		
            // Transport next Particle in the bank

        } // Particle bank is empty, end of history loop

        // Estimator history closeout
        for ( auto& E : Estimator ) { E->endHistory(); }
        // Start next history

    } // All histories are done, end of simulation loop
}

void Simulator_t::report()
{
    // Generate outputs
    std::ostringstream output;                       // Output text
    std::ofstream file( simName + " - output.txt" ); // .txt file
	
    // Header
    output << "\n";
    for ( int i = 0 ; i < simName.length()+6 ; i++ ) { output << "="; }
    output << "\n";
    output << "== " + simName + " ==\n";
    for ( int i = 0 ; i < simName.length()+6 ; i++ ) { output << "="; }
    output << "\n";
    output << "Number of histories: " << nhist << "\n";
    output << "Number of tracks: " << tracks << "\n";

    // Report tallies
    for ( auto& E : Estimator ) { E->report( output, tracks ); } // TrackTime is passed for F.O.M.
	
    // Print on monitor and file
    std::cout<< output.str();
    file<< output.str();
}
