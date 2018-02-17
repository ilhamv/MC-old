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
#include "H5Cpp.h"


//=============================================================================
// Set up
//=============================================================================

Simulator_t::Simulator_t( const std::string input_dir )
{
    io_dir = input_dir+"/";
    XML_input( io_dir, simulation_name, Nsample, ksearch, tdmc, Ncycle, 
               Npassive, Ecut_off, tcut_off, Fbank, Surfaces, Cells, Nuclides, 
               Materials, Estimators, Distribution_Double, Distribution_Point, 
               tdmc_time, tdmc_split );
    
    if(ksearch){
        k_estimator = std::make_shared<EstimatorK>(Ncycle, Ncycle-Npassive, 
                                                   Nsample);
    }
}


//=============================================================================
// Move particle
//=============================================================================

void Simulator_t::move_particle( Particle_t& P, const double l )
{
    P.move( l );
    if(ksearch) { k_estimator->estimate_TL(P,l); }
    if(tally){ 
        for( auto& e : P.cell()->estimators_TL ) { e->score( P, l ); }
    }
    Ntrack++;
}


//=============================================================================
// Collision
//=============================================================================
void Simulator_t::collision( Particle_t& P )
{
    if (ksearch) { k_estimator->estimate_C(P); }
    if(tally){ 
        for( auto& e : P.cell()->estimators_C ) { e->score( P, 0 ); }
    }
    
    P.cell()->collision( P, Pbank, ksearch, Fbank, k );			
}


//=============================================================================
// Cut-off and weight rouletting
//=============================================================================

void Simulator_t::cut_off( Particle_t& P )
{
    if ( P.energy() <= Ecut_off || P.time() >= tcut_off || P.weight() == 0.0 ){
        P.kill();
    }
    else{
        // Weight rouletting
        if( P.weight() < wr ){
            if( Urand() < P.weight() / ws ) { P.setWeight(ws); }
            else { P.kill(); }
        }
    }
}

//=============================================================================
// THE Simulation
//=============================================================================

void Simulator_t::start()
{
    // Simulation loop
    for ( icycle = 0; icycle < Ncycle ; icycle++ ){
        if ( icycle == Npassive ) { tally = true; }
        Sbank = Fbank; Fbank.reset();

        // Cycle loop
        for ( isample = 0 ; isample < Nsample ; isample++ ){
            Pbank.push( Sbank.getSource( Cells ) );
                    
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
                            P.set_tdmc(P.tdmc()+1);
                            cut_off( P );
                            // Time splitting
                            if(P.alive()){
                                P.setWeight(P.weight()/tdmc_split);
                                for( int i = 0; i < tdmc_split - 1; i++ ){
                                    Pbank.push(P);
                                }
                            }
                            continue;
                        }
                    }
                                    
                    // Hit surface?
                    if ( dcol > SnD.second ){	
                        P.set_surface_old(SnD.first);
                        move_particle( P, SnD.second );
                        SnD.first->hit( P, Cells, tally );
                        Split_Roulette( P, Pbank );
                    }
                    // Collide!!
                    else{
                        move_particle( P, dcol );
                        collision( P );
                    }        
                    cut_off( P );
                } 
            }

            // Estimator history closeout
            if (tally)
            { for ( auto& E : Estimators ) { E->end_history(); } }
            if (ksearch) { k_estimator->end_history(); }
        }

        // Estimator cycle closeout
        if (tally) { for ( auto& E : Estimators ) { E->end_cycle(Ntrack); } }
        if (ksearch){ 
            k_estimator->report_cycle(tally);
            k = k_estimator->k;
        }
        Ntrack = 0.0;
    
    } // All cycles are done, end of simulation loop
    for ( auto& E : Estimators ) { E->end_simulation(); }
}


//=============================================================================
// Report results
//=============================================================================

void Simulator_t::report()
{
    // Generate outputs
    std::ostringstream output;                       // Output text
    std::ofstream file( io_dir + "output.txt" ); // .txt file

    io_dir += "output.h5";
    H5std_string FILE_NAME(io_dir);
    H5::H5File output_H5(FILE_NAME, H5F_ACC_TRUNC);
	
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
    
    // Report tallies
    for ( auto& E : Estimators ) { E->report( output, output_H5 ); }
	
    // Print on monitor and file
    std::cout<< output.str();
    file<< output.str();
}
