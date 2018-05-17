#include <memory>  

#include "simulator.h"
#include "Algorithm.h"


//=============================================================================
// Handler
//=============================================================================

void Simulator::start()
{
    // Cycle loop
    for( icycle = 0; icycle < Ncycle ; icycle++ ){
        if( icycle == Npassive ) { tally = true; }
        Sbank = Fbank; Sbank.set_up(); Fbank.reset();

        // Sample loop
        for( isample = 0 ; isample < Nsample ; isample++ ){
            Pbank.push_back( Sbank.get_source() );
            
            while( !Pbank.empty() ){
//std::cout<<isample<<"    "<<Pbank.size()<<"\n";
                Particle P = Pbank.back(); Pbank.pop_back();
                // Particle random walk
                random_walk( P );
                // Particle control: Combing
                if( comb ) { particle_comb( Pbank ); }
            }

            // Estimator history closeout
            if( tally ){ 
                for ( auto& E : Estimators ) { E->end_history(); } 
            }
            if( ksearch ) { k_estimator->end_history(); }
        }

        // Estimator cycle closeout
        if( tally ){ for ( auto& E : Estimators ) { E->end_cycle(); } }
        if( ksearch ){ 
            k_estimator->report_cycle(tally);
            k = k_estimator->k;
        }
    } 
    
    // End of simulation closeout
    for( auto& E : Estimators ) { E->end_simulation(); }
}
