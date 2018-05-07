#include <memory>  

#include "simulator.h"
#include "Algorithm.h"
#include "Random.h"


//=============================================================================
// fixed source collision
//=============================================================================

void Simulator::collision_fixed_source( Particle& P, const double bank_nu, 
                                    const std::shared_ptr<Nuclide>& N_fission )
{
    if( Urand() > N_fission->beta(P.energy()) ){
        // Prompt
        for ( int i = 0 ; i < bank_nu ; i++ ){
            Particle P_new( P.pos(), isotropic_direction.sample(),
                            N_fission->fission()->Chi(P.energy()), P.time(),
                            1.0, P.tdmc(), P.cell() );
            push_particle_bank(P_new);
        }
    } else{
        // Delayed
        if(tdmc){
            // Combined and forced decay
            for( int i = 0 ; i < bank_nu ; i++ ){
                Particle P_new = forced_decay( P, N_fission, P.time(), 
                                            tdmc_time[P.tdmc()] - P.time(), 
                                            P.tdmc() );
                weight_roulette( P_new );
                if( P_new.alive() ) { push_particle_bank( P_new ); }
                for( int j = P.tdmc(); j < tdmc_time.size()-1; j++ ){
                    P_new = forced_decay(P, N_fission, tdmc_time[j],
                                         tdmc_interval[j+1], j+1);
                    weight_roulette( P_new );
                    if( P_new.alive() ) {push_particle_bank( P_new ); }
                }
            }
        } else{
            int cg, p_tdmc;
            double p_energy, p_time;
            // Precursor group cg
            const double xi = Urand();
            double sum = 0;
            for( int i = 0; i < 6; i++ ){
                sum += N_fission->fission()->fraction(i);
                if( sum > xi ){ cg = i; break; }
            }
            for( int i = 0 ; i < bank_nu ; i++ ){
                p_energy = N_fission->fission()->ChiD(cg,P.energy());
                p_time = exponential_sample( 
                        N_fission->fission()->lambda(cg) );
                for( int j = P.tdmc(); j < tdmc_time.size(); j++ ){
                    if( p_time < tdmc_time[P.tdmc()] ){
                        p_tdmc = j;
                        Particle P_new( P.pos(), isotropic_direction.sample(), 
                                        p_energy, p_time, 1.0, p_tdmc, 
                                        P.cell() );
                        push_particle_bank(P_new);
                    }
                }
            }
        }
    }
}
