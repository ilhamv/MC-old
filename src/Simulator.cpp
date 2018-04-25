#include <vector>  
#include <iostream>
#include <cstring> 
#include <memory>  
#include <stack>   
#include <cmath>   
#include <numeric>

#include "Simulator.h"
#include "Algorithm.h"
#include "Random.h"
#include "eigen3-hdf5.hpp"


//=============================================================================
// Test Point
//=============================================================================

bool Simulator::test_point( const Point& p, const std::shared_ptr<Cell>& C )
{
    // Loop over surfaces in cell, if not on correct side return false
    for ( const auto& S : C->surfaces() ) {
    	if ( S.first->eval( p ) * S.second < 0 ) { return false; }  
    }
    return true;
}


//=============================================================================
// Push Particle Bank
//=============================================================================
void Simulator::push_particle_bank( const Particle& P )
{
    Pbank.push(P);
}
        

//==============================================================================
// Add Fission Source
//==============================================================================
void Simulator::add_fission_source( const std::shared_ptr<Source>& S,
                                    const double i )
{
    Fbank.add_source( S, i );
}

//=============================================================================
// Search Cell
//=============================================================================

std::shared_ptr<Cell> Simulator::search_cell( const Point& p )
{
    for( const auto& C : Cells ){
        if ( test_point( p, C ) ){ return C; }
    }
    std::cout<< "[WARNING] A particle is lost:\n( x, y, z )  (" << p.x << ", " 
             << p.y << ", " << p.z << " )\n";
    std::exit(EXIT_FAILURE);
}


//=============================================================================
// Collision Distance
//=============================================================================

double Simulator::collision_distance( const Particle& P )
{
    if ( P.cell()->material() ){ 
        return exponential_sample( P.cell()->material()->SigmaT( P.energy() ) );
    }
    // Vacuum --> return sligthly less than very large number
    //            to ensure collision (kill) if no surface intersection
    else { return MAX_float_less; }
}

 
//=============================================================================
// Surface Intersect
//=============================================================================

std::pair<std::shared_ptr<Surface>, double> Simulator::surface_intersect
    ( const Particle& P ) 
{
    double dist = MAX_float;
    std::shared_ptr< Surface > S = nullptr;
    for ( const auto& s : P.cell()->surfaces() ) {
    	double d = s.first->distance( P );
    	if ( d < dist ){ 
	    dist = d;
	    S    = s.first;
	}
    }
    return std::make_pair( S, dist ); 
}

//=============================================================================
// Move particle
//=============================================================================

void Simulator::move_particle( Particle& P, const double l )
{
    P.move( l );
    if(ksearch) { k_estimator->estimate_TL(P,l); }
    if(tally){ 
        for( auto& e : P.cell()->estimators_TL ) { e->score( P, l ); }
    }
    Ntrack++;
}


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
    return Particle( P.pos(), isotropic.sample(), 
                     p_energy, p_time, p_weight, p_tdmc, 
                     P.cell() );
}


//=============================================================================
// Collision
//=============================================================================

void Simulator::collision( Particle& P )
{
    // Vacuum --> Kill particle at collision
    if( !P.cell()->material() ){ return P.kill(); }
    
    // Collision tallies
    if(tally){ 
        for( auto& e : P.cell()->estimators_C ) { e->score( P, 0 ); }
    }
   
    // Get material
    std::shared_ptr<Material> M = P.cell()->material();
    
    // Implicit Fission
    const double bank_nu = std::floor( P.weight() / k * M->nuSigmaF(P.energy())
                                       / M->SigmaT(P.energy()) + Urand() );

    std::shared_ptr<Nuclide> N_fission = M->nuclide_nufission( P.energy() );

    if(ksearch){ 
        // k-eigenvalue
        if( Urand() > N_fission->beta(P.energy()) ){
            // Prompt energy spectrum
            for ( int i = 0 ; i < bank_nu ; i++ ){
                Particle P_new( P.pos(), isotropic.sample(),
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
                Particle P_new( P.pos(), isotropic.sample(),
                                p_energy, P.time(), 1.0, P.tdmc(), P.cell() );
                add_fission_source( std::make_shared<SourceDelta>(P_new),1.0);
            }
        }
    } else{
        // Fixed Source
        if( Urand() > N_fission->beta(P.energy()) ){
            // Prompt
            for ( int i = 0 ; i < bank_nu ; i++ ){
                Particle P_new( P.pos(), isotropic.sample(),
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
                            Particle P_new( P.pos(), isotropic.sample(), 
                                            p_energy, p_time, 1.0, p_tdmc, 
                                            P.cell() );
                            push_particle_bank(P_new);
                        }
                    }
                }
            }
        }
    }
    
    //=========================================================================
    // ksearch estimate
    //=========================================================================

    if(ksearch){
        k_estimator->estimate_C(P);
        if( bank_nu > 0 ){ k_estimator->entropy->score( P.pos(), bank_nu ); }
    }
    
    // Implicit Absorption
    const double implicit = M->SigmaC(P.energy()) + M->SigmaF(P.energy());
    P.set_weight( P.weight() * ( M->SigmaT(P.energy()) - implicit ) 
                  / M->SigmaT(P.energy()) );

    std::shared_ptr<Nuclide> N_scatter = M->nuclide_scatter( P.energy() );
    if(!N_scatter){ return; }
    
    // The only analog reaction
    N_scatter->scatter()->sample( P );
}


//=============================================================================
// Surface Hit
//=============================================================================
        
void Simulator::surface_hit( Particle& P, const std::shared_ptr<Surface>& S )
{
    P.set_surface_old(S);
    if ( S->bc() == 0 ){
	P.move( EPSILON_float );
	P.set_cell( search_cell( P.pos() ) );
    } else if ( S->bc() == -1 ){
	P.kill();
        P.set_cell( P.cell() );
    } else{
	S->reflect( P );
	P.move( EPSILON_float );
        P.set_cell( P.cell() );
    }
    if (tally){
	for ( auto& e : S->estimators ){ 
            e->score( P, 0.0 ); 
        }
    }
}


//=============================================================================
// Weight Roulette
//=============================================================================

void Simulator::weight_roulette( Particle& P )
{
    // Weight rouletting
    if( P.weight() < wr ){
        if( Urand() < P.weight() / ws ) { P.set_weight(ws); }
        else { P.kill(); }
    }
}


//=============================================================================
// THE Simulation
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
}


//=============================================================================
// Report results
//=============================================================================

void Simulator::report( H5::H5File& output )
{
    // H5 tools
    H5::DataSet dataset;
    H5::Group group;
    H5::DataSpace space_scalar(H5S_SCALAR);
    H5::DataType type_ull    = H5::PredType::NATIVE_ULLONG;
    H5::DataType type_double = H5::PredType::NATIVE_DOUBLE;
    H5::StrType type_str(0, H5T_VARIABLE);

    // Summary
    group = output.createGroup("/summary");
    dataset = group.createDataSet( "Ncycle", type_ull, space_scalar );
    dataset.write(&Ncycle, type_ull);
    dataset = group.createDataSet( "Nsample", type_ull, space_scalar );
    dataset.write(&Nsample, type_ull);
    dataset = group.createDataSet( "Npassive",type_ull, space_scalar );
    dataset.write(&Npassive, type_ull);
    dataset = group.createDataSet( "Ntrack",type_ull, space_scalar );
    dataset.write(&Ntrack, type_ull);
    dataset = group.createDataSet( "mode", type_str, space_scalar );
    dataset.write(mode, type_str);
    group = output.createGroup("/summary/survival_roulette");
    dataset = group.createDataSet( "wr", type_double, space_scalar);
    dataset.write(&wr, type_double);
    dataset = group.createDataSet( "ws", type_double, space_scalar);
    dataset.write(&ws, type_double);
    if(tdmc){
        group = output.createGroup("/summary/tdmc");
        dataset = group.createDataSet( "split", type_ull, space_scalar);
        dataset.write(&tdmc_split, type_ull);
        hsize_t dimsv[1]; dimsv[0] = tdmc_time.size();
        H5::DataSpace data_spacev(1,dimsv);
        dataset = group.createDataSet( "time", type_double, data_spacev);
        dataset.write(tdmc_time.data(), type_double);
    }

    // Report estimators
    for ( auto& E : Estimators ) { E->report( output ); }
    if(ksearch){k_estimator->report(output);}

    // TRMM report
    if( trmm ){
        const int score_N = trmm_estimator_simple->score_size();
        // TRM
        int G = trmm_estimator_simple->tally_size() / score_N;
        int J = 6;
        int trm_N = G+J;
        std::vector<double> TRM(trm_N*trm_N);

        // M
        int idx = 0;
        for( int f = 0; f < G; f++ ){
            for( int i = 0; i < G; i++ ){
                if( i == f ){
                    TRM[idx]  = -trmm_estimator_simple->tally(i).mean;
                    TRM[idx] +=  trmm_estimator_scatter->tally(i+i*G).mean;
                    TRM[idx] +=  trmm_estimator_fission_prompt->tally(i+i*G).mean;
                    TRM[idx] /=  trmm_estimator_simple->tally(i+G).mean;
                } else{
                    TRM[idx]  = trmm_estimator_scatter->tally(f+i*G).mean;
                    TRM[idx] += trmm_estimator_fission_prompt->tally(f+i*G).mean;
                    TRM[idx] /= trmm_estimator_simple->tally(i+G).mean;
                }
                idx++;
            }
            idx += J;
        }

        // D
        for( int j = 0; j < J; j++ ){
            for( int g = 0; g < G; g++ ){
                TRM[idx] = trmm_estimator_simple->tally(2*G+j*G+g).mean;
                TRM[idx] /=  trmm_estimator_simple->tally(G+g).mean;
                idx++;
            }
            idx += J;
        }

        // L
        std::vector<double> lambda(J);
        for( int j = 0; j < J; j++ ){
            double num = 0.0;
            double denom = 0.0;
            idx = trm_N*(G+j) + G + j;
            for( int g = 0; g < G; g++ ){
                num += trmm_estimator_simple->tally(2*G+j*G+g).mean;
                denom += trmm_estimator_simple->tally(2*G+J*G+j*G+g).mean;
            }
            lambda[j] = num / denom;
            TRM[idx] = -lambda[j];
        }

        // P
        for( int g = 0; g < G; g++ ){
            for( int j = 0; j < J; j++ ){
                double num = 0.0;
                double denom = 0.0;
                for( int gp = 0; gp < G; gp++ ){
                    num += trmm_estimator_fission_delayed[j]->tally(g+gp*G).mean;
                    denom += trmm_estimator_simple->tally(2*G+j*G+gp).mean;
                }
                idx = trm_N*g + G + j;
                TRM[idx] = num / denom * lambda[j];
            }
        }



        hsize_t dimsM[2]; dimsM[0] = trm_N; dimsM[1] = trm_N;
        H5::DataSpace data_spaceM(2,dimsM);
        dataset = output.createDataSet( "TRM", type_double, data_spaceM);
        dataset.write(TRM.data(), type_double);

        // Inverse speed
        std::vector<double> invspeed(G);
        for( int i = 0; i < G; i++ ){
            invspeed[i] = trmm_estimator_simple->
                                  tally( i + (score_N-1) * G).mean
                          / trmm_estimator_simple->tally(i+G).mean;
        }
        hsize_t dimsv[1]; dimsv[0] = G;
        H5::DataSpace data_spacev(1,dimsv);
        dataset = output.createDataSet("inverse_speed",type_double,data_spacev);
        dataset.write(invspeed.data(), type_double);
    }
}
