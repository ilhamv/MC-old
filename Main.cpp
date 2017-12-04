#include <vector>  // vector
#include <iostream>// cout
#include <cstring> // string, strcmp
#include <memory>  // shared_ptr, make_shared
#include <stack>   // stack
#include <cmath>   // exp
#include <sstream> // ostringstream
#include <fstream> // ifstream

//#include "boost/filesystem.hpp" // exists
#include "VReduction.h"         // Split_Roulette
#include "Const.h"              // MAX
#include "pugixml.hpp"
#include "Geometry.h"
#include "Particle.h"
#include "Distribution.h"
#include "Source.h"
#include "Nuclide.h"
#include "Material.h"
#include "Reaction.h"
#include "Estimator.h"
#include "XSec.h"
#include "XMLparser.h"


int main( int argc, char* argv[] )
{
    // Check # of arguments passed
    if ( argc == 1 ) { std::cout<< "Please provide input file\n"; return 1; }

    // Variable declarations
    std::string                                  input_file = argv[1]; // Input file name 
    std::string                                  input_dir;            // Input file directory*/
    std::string                                  simName;              // Simulation name
    unsigned long long                           nhist;                // Number of particle samples
    double                                       Ecut_off  = 0.0;      // Energy cut-off
    double                                       tcut_off  = MAX;      // Time cut-off
    unsigned long long                           trackTime = 0;        // "Computation time" (particle track) for variance reduction
    Source_Bank                                  Sbank;                // Source bank
    std::stack  < Particle_t >                   Pbank;                // Particle bank
    std::vector < std::shared_ptr<Surface_t>   > Surface;              // Surfaces
    std::vector < std::shared_ptr<Cell_t>    > Cell;                 // Cells
    std::vector < std::shared_ptr<Nuclide_t>   > Nuclide;              // Nuclides
    std::vector < std::shared_ptr<Material_t>  > Material;             // Materials
    std::vector < std::shared_ptr<Estimator_t> > Estimator;            // Estimators  	
    // User-defined distributions
    std::vector < std::shared_ptr<Distribution_t<double>> > double_distributions;
    std::vector < std::shared_ptr<Distribution_t<int>>    > int_distributions;
    std::vector < std::shared_ptr<Distribution_t<Point_t>>> point_distributions;
    	
    // XML input parser	
    XML_input( input_file, simName, nhist, Ecut_off, tcut_off, Sbank, Surface, Cell, Nuclide, Material, Estimator, double_distributions, int_distributions, point_distributions );
	
    std::cout<<"\nSimulation setup done,\nNow running the simulation...\n\n";
    std::cout.flush();

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
					
                    // Accumulate "computation time"
                    trackTime++;
                }
				
                // Collide!!
                else
                {
                    // Move particle to collision site and sample the collision and tally if there is any cell tally
                    C->moveParticle( P, dcol );
                    C->collision( P, Pbank );
					
                    // Accumulate "computation time"
                    trackTime++;
                }
			
                // Cut-off working particle?
                if ( P.energy() < Ecut_off || P.time() > tcut_off ) { P.kill();}

            } // Particle is dead, end of particle loop		
            // Transport next Particle in the bank

        } // Particle bank is empty, end of history loop

        // Estimator history closeout
        for ( auto& E : Estimator ) { E->endHistory(); }
        // Start next history

    } // All histories are done, end of simulation loop


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
    output << "Track time: " << trackTime << "\n";

    // Report tallies
    for ( auto& E : Estimator ) { E->report( output, trackTime ); } // TrackTime is passed for F.O.M.
	
    // Print on monitor and file
    std::cout<< output.str();
    file<< output.str();

    return 0;
}
