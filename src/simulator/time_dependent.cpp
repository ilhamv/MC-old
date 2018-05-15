#include <vector>  
#include <memory>  
#include <cmath>   
#include <numeric> // accumulate
#include <Eigen/Dense>

#include "simulator.h"
#include "Random.h"
#include "eigen3-hdf5.hpp"


//=============================================================================
// Forced delay neutron decay
//=============================================================================

Particle Simulator::forced_decay( const Particle& P, 
                                  const std::shared_ptr<Nuclide>& N, 
                                  const double initial, const double interval,
                                  const int p_tdmc )
{
    double p_time, p_energy;
    int cg;
    double p_weight = 0.0;
    std::vector<double> prob(6);
    p_time = initial + Urand() * interval;
    for( int k = 0; k < 6; k++ ){
        prob[k] += N->fission()->f_lambda(k) * std::exp(-N->fission()->lambda(k)
                                                        *(p_time - P.time()) );
        p_weight += prob[k];
    } 
    p_weight *= interval;
    
    const double xi = Urand() 
                      * std::accumulate( prob.begin(), prob.end(), 0.0 );
    double sum = 0;
    for( int k = 0; k < 6; k++ ){
        sum += prob[k];
        if( sum > xi ){ cg = k; break; }
    }

    p_energy = N->fission()->ChiD(cg,P.energy());
    return Particle( P.pos(), isotropic_direction.sample(), 
                     p_energy, p_time, p_weight, p_tdmc, 
                     P.cell() );
}

//=============================================================================
// Time boundary hit
//=============================================================================

void Simulator::time_hit( Particle& P )
{
    P.increment_tdmc();
    if( P.time() == tdmc_time.back() ) { P.kill(); }
}
