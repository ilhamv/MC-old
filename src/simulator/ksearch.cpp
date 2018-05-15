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
// ksearch implicit fission
//=============================================================================

void Simulator::implicit_fission_ksearch( const Particle& P, 
        const double bank_nu, const std::shared_ptr<Nuclide>& N_fission )
{
    // Get energy spectrum: Prompt or Delayed? If Delayed, which group?
    std::shared_ptr<Distribution<double>> Chi;
    if( Urand() > N_fission->beta(P.energy()) ){
        Chi = N_fission->fission()->r_Chi;
    } else{
        double p_energy;
        // Precursor group cg
        int cg;
        double sum = 0;
        const double xi = Urand();
        for( int i = 0; i < 6; i++ ){
            sum += N_fission->fission()->fraction(i);
            if( sum > xi ){ cg = i; break; }
        }
        Chi = N_fission->fission()->r_ChiD[cg];
    }

    // Generate fission neutrons
    for ( int i = 0 ; i < bank_nu ; i++ ){
        Particle P_new( P.pos(), isotropic_direction.sample(),
                        N_fission->fission()->Chi(P.energy()), P.time(),
                        1.0, P.tdmc(), P.cell() );
        add_fission_source( std::make_shared<SourceDelta>(P_new),1.0);
    }
}
