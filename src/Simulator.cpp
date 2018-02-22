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
               tdmc_time, tdmc_split, trmm, trmm_estimator );
    
    if(ksearch){
        mode = "k-eigenvalue";
        k_estimator = std::make_shared<EstimatorK>(Ncycle, Ncycle-Npassive, 
                                                   Nsample);
    }
    if(tdmc) { mode = "time-dependent"; }
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
        if (tally) { for ( auto& E : Estimators ) { E->end_cycle(); } }
        if (ksearch){ 
            k_estimator->report_cycle(tally);
            k = k_estimator->k;
        }
    
    } // All cycles are done, end of simulation loop
    for ( auto& E : Estimators ) { E->end_simulation(); }
}

//=============================================================================
// Set TRM
//=============================================================================
void Simulator_t::set_TRM()
{;
}


//=============================================================================
// Report results
//=============================================================================

void Simulator_t::report()
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
    unsigned long long trm_N = trmm_estimator[0]->tally_size()/2;
    Eigen::MatrixXd TRM,TRM_real;
    TRM = Eigen::MatrixXd(trm_N,trm_N);

    for( int f = 0; f < trm_N; f++ ){
        for( int i = 0; i < trm_N; i++ ){
            if( i == f ){
                TRM(i,i) = -trmm_estimator[0]->tally(i).mean /
                    trmm_estimator[0]->tally(i+trmm_estimator[0]->idx_factor[0]).mean;
                TRM(i,i) += trmm_estimator[1]->tally(i+i*trm_N).mean /
                    trmm_estimator[0]->tally(i+trmm_estimator[0]->idx_factor[0]).mean;
                TRM(i,i) += trmm_estimator[2]->tally(i+i*trm_N).mean /
                    trmm_estimator[0]->tally(i+trmm_estimator[0]->idx_factor[0]).mean;
            } else{
                TRM(f,i) = trmm_estimator[1]->tally(f+i*trm_N).mean /
                    trmm_estimator[0]->tally(i+trmm_estimator[0]->idx_factor[0]).mean;
                TRM(f,i) += trmm_estimator[2]->tally(f+i*trm_N).mean /
                    trmm_estimator[0]->tally(i+trmm_estimator[0]->idx_factor[0]).mean;
            }
        }
    }
    TRM_real = TRM.real();

    // Solve eigenvalue fo TRM
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
