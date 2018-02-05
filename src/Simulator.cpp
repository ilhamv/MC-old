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
    XML_input( input_file, simulation_name, Nsample, ksearch, Ncycle, Npassive, Ecut_off, tcut_off, Fbank, Surface, Cell, Nuclide, Material, estimator, Distribution_Double, Distribution_Point );
    
    k_cycle.resize(Ncycle, 0.0);
}

void Simulator_t::start()
{
    // Simulation loop
    for ( unsigned long long icycle = 0; icycle < Ncycle ; icycle++ ){
        if ( icycle == Npassive ) { tally = true; }
        Sbank = Fbank; Fbank.reset();

        // Cycle loop
        for ( unsigned long long isample = 0 ; isample < Nsample ; isample++ ){
            Pbank.push( Sbank.getSource( Cell ) );
                    
            // History loop
            while ( !Pbank.empty() ){
                Particle_t              P = Pbank.top(); // Working particle
                std::shared_ptr<Cell_t> C = P.cell();    // Working cell
                Pbank.pop();

                // Particle loop
                while ( P.alive() ){
                    std::pair< std::shared_ptr< Surface_t >, double > SnD; // To hold nearest surface and its distance
                                    
                    // Determine nearest surface and its distance
                    SnD = C->surface_intersect( P );

                    // Determine collision distance
                    const double dcol = C->collision_distance( P.energy() );
                                    
                    // Hit surface?
                    if ( dcol > SnD.second ){	
                        // Surface hit! Move particle to surface, tally if there is any Cell Tally
                        C->moveParticle( P, SnD.second, tally );

                        // Implement surface hit:
                        // 	Reflect angle for reflective surface
                        // 	Cross the surface (move epsilon distance)
                        // 	Search new particle cell for transmission surface
                        // 	Tally if there is any Surface Tally
                        // 	Note: particle weight and working cell are not updated yet
                        SnD.first->hit( P, Cell, tally );

                        // Splitting & Roulette Variance Reduction Technique
                        // 	More important : split
                        // 	Less important : roulette
                        // 	Same importance: do nothing
                        // 	Note: old working cell has the previous cell importance and will be updated
                        Split_Roulette( C, P, Pbank );
                    }
                                    
                    // Collide!!
                    else{
                        // Move particle to collision site and sample the collision and tally if there is any cell tally
                        C->moveParticle( P, dcol, tally );
                        
                        // New estimate k
                        if (ksearch)
                        { k_cycle[icycle] += C->nuSigmaF(P.energy()) * P.weight() / C->SigmaT(P.energy()); }
                        C->collision( P, Pbank, ksearch, Fbank, k );			
                    }
                            
                    // add # of tracks
                    tracks++;

                    // Cut-off or kill working particle?
                    if ( P.energy() < Ecut_off || P.time() > tcut_off || P.weight() == 0.0 ) { P.kill();}
                    else{
                        // Weight rouletting
                        if( P.weight() < wr ){
                            if( Urand() < P.weight() / ws ) { P.setWeight(ws); }
                            else { P.kill(); }
                        }
                    }

                } // Particle is dead, end of particle loop		
                // Transport next Particle in the bank

            } // Particle bank is empty, end of history loop

            // Estimator history closeout
            if (tally)
            { for ( auto& E : estimator ) { E->end_history(); } }

            // Start next history

        } // All histories are done, end of cycle loop

        // Estimator history closeout
        if (tally)
        { for ( auto& E : estimator ) { E->end_cycle(Nsample, tracks); } }

        if (ksearch){
            k_cycle[icycle] /= Nsample;
            k = k_cycle[icycle];
        }

        // Start next cycle
        tracks = 0.0;
    
    } // All cycles are done, end of simulation loop
    { for ( auto& E : estimator ) { E->end_simulation(Ncycle); } }
}

void Simulator_t::report()
{
    // Generate outputs
    std::ostringstream output;                       // Output text
    std::ofstream file( simulation_name + " - output.txt" ); // .txt file
	
    // Header
    output << "\n";
    for ( int i = 0 ; i < simulation_name.length()+6 ; i++ ) { output << "="; }
    output << "\n";
    output << "== " + simulation_name + " ==\n";
    for ( int i = 0 ; i < simulation_name.length()+6 ; i++ ) { output << "="; }
    output << "\n";
    output << "Number of passive cycles   : " << Npassive << "\n";
    output << "Number of active cycles    : " << Ncycle - Npassive << "\n";
    output << "Number of samples per cycle: " << Nsample << "\n\n";
    
    if (ksearch ) { 
        k = 0.0;
        for( int i = Npassive; i < Ncycle; i++ ) { k += k_cycle[i]; }
        k /= (Ncycle - Npassive);
        output << "k-eigenvalue: " << k << "\n"; 
    }

    // Report tallies
    for ( auto& E : estimator ) { E->report( output ); }
	
    // Print on monitor and file
    std::cout<< output.str();
    file<< output.str();
}
