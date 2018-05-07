#include <memory>

#include "simulator.h"
#include "Random.h"


//==============================================================================
// add fission source
//==============================================================================
void Simulator::add_fission_source( const std::shared_ptr<Source>& S,
                                    const double i )
{
    Fbank.add_source( S, i );
}

//=============================================================================
// ksearch collision
//=============================================================================

void Simulator::collision_ksearch( Particle& P, const double bank_nu, 
                                   const std::shared_ptr<Nuclide>& N_fission )
{
    if( Urand() > N_fission->beta(P.energy()) ){
        // Prompt energy spectrum
        for ( int i = 0 ; i < bank_nu ; i++ ){
            Particle P_new( P.pos(), isotropic_direction.sample(),
                            N_fission->fission()->Chi(P.energy()), P.time(),
                            1.0, P.tdmc(), P.cell() );
            add_fission_source( std::make_shared<SourceDelta>(P_new),1.0);
        }
    } else{
        // Delayed energy spectrum
        double p_energy;
        // Precursor group cg
        int cg;
        double sum = 0;
        const double xi = Urand();
        for( int i = 0; i < 6; i++ ){
            sum += N_fission->fission()->fraction(i);
            if( sum > xi ){ cg = i; break; }
        }
        for( int i = 0 ; i < bank_nu ; i++ ){
            p_energy = N_fission->fission()->ChiD(cg,P.energy());
            Particle P_new( P.pos(), isotropic_direction.sample(),
                            p_energy, P.time(), 1.0, P.tdmc(), P.cell() );
            add_fission_source( std::make_shared<SourceDelta>(P_new),1.0);
        }
    }
    
    // Tallies
    k_estimator->estimate_C(P);

    // Entropy
    if( bank_nu > 0 ){ k_estimator->entropy->score( P.pos(), bank_nu ); }
}
