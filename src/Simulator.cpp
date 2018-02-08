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
    XML_input( input_file, simulation_name, Nsample, ksearch, tdmc, Ncycle, Npassive,
               Ecut_off, tcut_off, Fbank, Surface, Cell, Nuclide, Material,
               estimator, Distribution_Double, Distribution_Point, tdmc_time );
    
    if(ksearch){
        k_cycle.resize(Ncycle, 0.0);
        k_C = std::make_shared<Estimator>("k-eigenvalue: collision");
        k_C->add_score(std::make_shared<ScoreNuFission>("k_eigenvalue: collision",
            std::make_shared<ScoreKernelCollision>()));
        k_C->initialize_tallies();
        estimator.push_back(k_C);
    }
}

// Move particle
void Simulator_t::move_particle( Particle_t& P, const double l )
{
    P.move( l );
    // Score track length estimators
    if(tally) { for( auto& e : P.cell()->estimators ) { e->score( P, l ); } }
}

void Simulator_t::collision( Particle_t& P )
{;
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
                Particle_t P = Pbank.top(); // Working particle
                Pbank.pop();

                // Particle loop
                while ( P.alive() ){
                     // To hold nearest surface and its distance
                    std::pair< std::shared_ptr< Surface_t >, double > SnD;
                                    
                    // Determine nearest surface and its distance
                    SnD = P.cell()->surface_intersect( P );

                    // Determine collision distance
                    double dcol = P.cell()->collision_distance( P.energy() );

                    // Exceeding TDMC time boundary?
                    if(tdmc){
                        double dbound = (tdmc_time[P.tdmc()] - P.time()) 
                                        * P.speed();
                        if( std::min(SnD.second,dcol) > dbound ){
                            move_particle( P, dbound );
                            Fbank.addSource( std::make_shared<Delta_Source>
                                             ( P.pos(), P.dir(), P.energy(), 
                                               P.weight(), P.time() ) );
                            move_particle( P, dbound );
                            goto tdmc_skip;
                        }
                    }
                                    
                    // Hit surface?
                    if ( dcol > SnD.second ){	
                        // Surface hit!
                        P.set_surface_old(SnD.first);

                        // Move particle to surface, 
                        //   tally if there is any cell tally
                        move_particle( P, SnD.second );

                        // Implement surface hit:
                        //   Reflect angle for reflective surface
                        //   Cross the surface (move epsilon distance)
                        //   Search new particle cell for transmission surface
                        //   Tally if there is any surface tally
                        SnD.first->hit( P, Cell, tally );

                        // Splitting & Roulette Variance Reduction Technique
                        //   More important : split
                        //   Less important : roulette
                        //   Same importance: do nothing
                        Split_Roulette( P, Pbank );
                    }
                                    
                    // Collide!!
                    else{
                        // Move particle to collision site
                        //   tally if there is any surface tally
                        move_particle( P, dcol );
                        
                        // New estimate k
                        if (ksearch){ 
                            k_cycle[icycle] += P.cell()->nuSigmaF(P.energy()) * P.weight() / P.cell()->SigmaT(P.energy());
                            if (tally) { k_C->score(P,dcol); }
                        }
                        
                        P.cell()->collision( P, Pbank, ksearch, Fbank, k );			
                    }
                            
                    tdmc_skip:
                    // add # of tracks
                    tracks++;

                    // Cut-off or kill working particle?
                    if ( P.energy() <= Ecut_off || P.time() >= tcut_off || P.weight() == 0.0 ) { P.kill();}
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
            std::cout<<icycle<<"  "<<k<<"\n";
        }

        // Start next cycle
        tracks = 0.0;
    
    } // All cycles are done, end of simulation loop
    { for ( auto& E : estimator ) { E->end_simulation(Ncycle-Npassive); } }
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
