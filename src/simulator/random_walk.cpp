#include <memory>  

#include "simulator.h"
#include "Algorithm.h"

//=============================================================================
// Random walk
//=============================================================================

void Simulator::start()
{
    // Simulation loop
    for ( icycle = 0; icycle < Ncycle ; icycle++ ){
        if ( icycle == Npassive ) { tally = true; }
        Sbank = Fbank; Sbank.set_up(); Fbank.reset();

        // Cycle loop
        for ( isample = 0 ; isample < Nsample ; isample++ ){
            Pbank.push( Sbank.get_source() );
                    
            // History loop
            while ( !Pbank.empty() ){
                Particle P = Pbank.top(); // Working particle
                Pbank.pop();

                // Particle loop
                while ( P.alive() ){
                    // Determine nearest surface and its distance
                    const std::pair< std::shared_ptr< Surface >, double > SnD =
                        surface_intersect( P );

                    // Determine collision distance
                    double dcol = collision_distance( P );

                    // Exceeding TDMC time boundary?
                    if(tdmc){
                        double dbound = (tdmc_time[P.tdmc()] - P.time()) 
                                        * P.speed();
                        if( std::min(SnD.second,dcol) > dbound ){
                            move_particle( P, dbound );
                            P.increment_tdmc();
                            if( P.time() == tdmc_time.back() ) { P.kill(); }
                            // Time splitting
                            if(P.alive()){
                                P.set_weight(P.weight()/tdmc_split);
                                for( int i = 0; i < tdmc_split - 1; i++ ){
                                    Particle P_new = P;
                                    weight_roulette( P_new );
                                    if( P_new.alive() ){ 
                                        push_particle_bank(P_new); 
                                    }
                                }
                                weight_roulette( P );
                            }
                            continue;
                        }
                    }
                                    
                    // Hit surface?
                    if ( dcol > SnD.second ){	
                        move_particle( P, SnD.second );
                        surface_hit( P, SnD.first );
                        split_roulette( P, Pbank );
                    }
                    // Collide!!
                    else{
                        move_particle( P, dcol );
                        collision( P );
                    }        
                    weight_roulette( P );
                } 
            }

            // Estimator history closeout
            if (tally)
            { for ( auto& E : Estimators ) { E->end_history(); } }
            if (ksearch) { k_estimator->end_history(); }
        }

        // Estimator cycle closeout
        if (tally) { for ( auto& E : Estimators ) { E->end_cycle(); } }
        if (ksearch){ 
            k_estimator->report_cycle(tally);
            k = k_estimator->k;
        }
    
    } // All cycles are done, end of simulation loop
    for ( auto& E : Estimators ) { E->end_simulation(); }

    // Save source?
    if( save_source ){

    }
}
