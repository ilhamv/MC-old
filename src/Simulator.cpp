#include <vector>  
#include <iostream>
#include <cstring> 
#include <memory>  
#include <stack>   
#include <cmath>   

#include "Simulator.h"
#include "Algorithm.h"
#include "Random.h"
#include "H5Cpp.h"


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
// Collision
//=============================================================================

void Simulator::collision( Particle& P )
{
    if(ksearch) { k_estimator->estimate_C(P); }
    
    if(tally){ 
        for( auto& e : P.cell()->estimators_C ) { e->score( P, 0 ); }
    }
    
    // Vacuum --> Kill particle at collision
    if( !P.cell()->material() ){ return P.kill(); }
    
    std::shared_ptr<Material> M = P.cell()->material();
    
    // Implicit Fission
    const double bank_nu = std::floor( P.weight() / k * M->nuSigmaF(P.energy())
                                       / M->SigmaT(P.energy()) + Urand() );

    std::shared_ptr<Nuclide> N_fission = M->nuclide_nufission( P.energy() );

    if(ksearch){ 
        // k-eigenvalue
        if( Urand() > N_fission->beta(P.energy()) ){
            // Prompt
            for ( int i = 0 ; i < bank_nu ; i++ ){
                Particle P_new( P.pos(), isotropic.sample(),
                                N_fission->fission()->Chi(P.energy()), P.time(),
                                1.0, P.tdmc(), P.cell() );
                add_fission_source( std::make_shared<SourceDelta>(P_new),1.0);
            }
        } else{
            // Delayed
            const double xi = Urand();
            double sum = 0;
            int cg, p_tdmc;
            double p_energy, p_time;
            // Precursor group cg
            for( int i = 0; i < 6; i++ ){
                sum += N_fission->fission()->fraction(i);
                if( sum > xi ){ cg = i; break; }
            }
            for( int i = 0 ; i < bank_nu ; i++ ){
                p_energy = N_fission->fission()->ChiD(cg,P.energy());
                p_time = exponential_sample( N_fission->fission()->lambda(cg) );
                for( int j = P.tdmc(); j < tdmc_time.size(); j++ ){
                    if( p_time < tdmc_time[P.tdmc()] ){
                        p_tdmc = j;
                        Particle P_new( P.pos(), isotropic.sample(),
                                     p_energy, p_time, 1.0, p_tdmc, P.cell() );
                        add_fission_source( std::make_shared<SourceDelta>(P_new)
                                                                          ,1.0);
                    }
                }
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
            const double xi = Urand();
            double sum = 0;
            int cg, p_tdmc;
            double p_energy, p_time;
            // Precursor group cg
            for( int i = 0; i < 6; i++ ){
                sum += N_fission->fission()->fraction(i);
                if( sum > xi ){ cg = i; break; }
            }
            for( int i = 0 ; i < bank_nu ; i++ ){
                p_energy = N_fission->fission()->ChiD(cg,P.energy());
                p_time = exponential_sample( N_fission->fission()->lambda(cg) );
                for( int j = P.tdmc(); j < tdmc_time.size(); j++ ){
                    if( p_time < tdmc_time[P.tdmc()] ){
                        p_tdmc = j;
                        Particle P_new( P.pos(), isotropic.sample(),
                                     p_energy, p_time, 1.0, p_tdmc, P.cell() );
                        push_particle_bank(P_new);
                    }
                }
            }
        }
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
    } else{
	S->reflect( P );
        P.set_cell( P.cell() );
	P.move( EPSILON_float );
    }
    if (tally){
	for ( auto& e : S->estimators ){ 
            e->score( P, 0.0 ); 
        }
    }
}


//=============================================================================
// Cut-off and weight rouletting
//=============================================================================

void Simulator::cut_off( Particle& P )
{
    if ( P.energy() <= Ecut_off || P.time() >= tcut_off || P.weight() == 0.0 ){
        P.kill();
    }
    else{
        // Weight rouletting
        if( P.weight() < wr ){
            if( Urand() < P.weight() / ws ) { P.set_weight(ws); }
            else { P.kill(); }
        }
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
        Sbank = Fbank; Fbank.reset();

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
                            cut_off( P );
                            // Time splitting
                            if(P.alive()){
                                P.set_weight(P.weight()/tdmc_split);
                                for( int i = 0; i < tdmc_split - 1; i++ ){
                                    push_particle_bank(P);
                                }
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
                    cut_off( P );
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

void Simulator::report()
{
    // H5 output
    io_dir += "output.h5";
    H5std_string FILE_NAME(io_dir);
    H5::H5File output(FILE_NAME, H5F_ACC_TRUNC);
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
    dataset = group.createDataSet( "cut_off-E", type_double, space_scalar);
    dataset.write(&Ecut_off, type_double);
    dataset = group.createDataSet( "cut_off-t", type_double, space_scalar);
    dataset.write(&tcut_off, type_double);
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

    if(!trmm){return;}

    // Set TRM
    unsigned long long trm_N = trmm_estimator_collision->tally_size()/2;
    Eigen::MatrixXd TRM,TRM_real;
    TRM = Eigen::MatrixXd(trm_N,trm_N);

    for( int f = 0; f < trm_N; f++ ){
        for( int i = 0; i < trm_N; i++ ){
            if( i == f ){
                TRM(i,i)  = -trmm_estimator_collision->tally(i).mean;
                TRM(i,i) +=  trmm_estimator_scatter->tally(i+i*trm_N).mean;
                TRM(i,i) +=  trmm_estimator_fission->tally(i+i*trm_N).mean;
                TRM(i,i) /=  trmm_estimator_collision->tally(i+trm_N).mean;
            } else{
                TRM(f,i)  = trmm_estimator_scatter->tally(f+i*trm_N).mean;
                TRM(f,i) += trmm_estimator_fission->tally(f+i*trm_N).mean;
                TRM(f,i) /= trmm_estimator_collision->tally(i+trm_N).mean;
            }
        }
    }
    TRM_real = TRM.real();

    // Solve eigenvalue of TRM
    Eigen::MatrixXcd phi_mode;
    Eigen::VectorXcd alpha;
    Eigen::EigenSolver<Eigen::MatrixXd> eSolve(TRM);
    phi_mode = eSolve.eigenvectors();
    alpha    = eSolve.eigenvalues();
    std::vector<double> alpha_real(trm_N);
    std::vector<double> alpha_imag(trm_N);
    for( int i = 0; i < trm_N; i++ ){
        alpha_real[i] = alpha[i].real();
        alpha_imag[i] = alpha[i].imag();
    }

    // Solve coefficients via initial condition
    Eigen::VectorXcd phi0;
    Eigen::VectorXcd A;
    phi0 = Eigen::VectorXcd::Zero(trm_N);
    phi0(trm_N-1) = 1.0 * std::sqrt( 14.1E6 * 191312955.067 ) * 100.0;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> dec(phi_mode);
    A = dec.solve(phi0);

    // Construct solution in time
    std::vector<double> t = {0.0, 3E-8, 15E-8, 4E-6, 1E-4};
    Eigen::MatrixXcd phi = Eigen::MatrixXcd::Zero(t.size(),trm_N);
    std::vector<double> phi_real(trm_N*t.size());

    unsigned long long idx = 0;
    for( int i = 0; i < t.size(); i++){
        for(int g = 0; g < trm_N; g++){
            for(int n = 0; n < trm_N; n++){
                phi(i,g) += A(n) * phi_mode(g,n) * std::exp( alpha[n] * t[i] );
            }
            phi_real[idx] = phi(i,g).real();
            idx++;
        }
    }
    
    // TRMM results
    group = output.createGroup("/TRMM");
    hsize_t dimsM[2]; dimsM[0] = trm_N; dimsM[1] = trm_N;
    H5::DataSpace data_spaceM(2,dimsM);
    dataset = group.createDataSet( "TRM", type_double, data_spaceM);
    dataset.write(TRM_real.data(), type_double);
    hsize_t dims[2]; dims[0] = t.size(); dims[1] = trm_N;
    H5::DataSpace data_spacev(2,dims);
    dataset = group.createDataSet( "flux", type_double, data_spacev);
    dataset.write(phi_real.data(), type_double);
    
    hsize_t dims_alpha[1]; dims_alpha[0] = trm_N;
    H5::DataSpace space_alpha(1,dims_alpha);
    group = group.createGroup("alpha");
    dataset = group.createDataSet( "real", type_double, space_alpha);
    dataset.write(alpha_real.data(), type_double);
    dataset = group.createDataSet( "imag", type_double, space_alpha);
    dataset.write(alpha_imag.data(), type_double);
}
