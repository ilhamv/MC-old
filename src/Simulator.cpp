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
    std::cout<<"Simulation done!\n\n";
}


//=============================================================================
// Report results
//=============================================================================

void Simulator_t::report()
{
    // H5 output
    std::cout<<"Creating output.h5...\n";
    io_dir += "output.h5";
    H5std_string FILE_NAME(io_dir);
    H5::H5File output(FILE_NAME, H5F_ACC_TRUNC);
    H5::DataSet dataset;
    H5::Group group;
    hsize_t dims[1]; dims[0] = 1;
    H5::DataSpace data_space(1,dims);
    H5::DataType type_ull    = H5::PredType::NATIVE_ULLONG;
    H5::DataType type_double = H5::PredType::NATIVE_DOUBLE;
    H5::DataType type_string = H5::StrType(H5::PredType::C_S1, 20);

    // Summary
    group = output.createGroup("/summary");
    dataset = group.createDataSet( "Ncycle", type_ull, data_space );
    dataset.write(&Ncycle, type_ull);
    dataset = group.createDataSet( "Nsample", type_ull, data_space );
    dataset.write(&Nsample, type_ull);
    dataset = group.createDataSet( "Npassive",type_ull, data_space );
    dataset.write(&Npassive, type_ull);
    dataset = group.createDataSet( "Ntrack",type_ull, data_space );
    dataset.write(&Ntrack, type_ull);
    dataset = group.createDataSet( "mode", type_string, data_space );
    dataset.write(mode.c_str(), type_string);
    dataset = group.createDataSet( "cut_off-E", type_double, data_space);
    dataset.write(&Ecut_off, type_double);
    dataset = group.createDataSet( "cut_off-t", type_double, data_space);
    dataset.write(&tcut_off, type_double);
    group = output.createGroup("/summary/survival_roulette");
    dataset = group.createDataSet( "wr", type_double, data_space);
    dataset.write(&wr, type_double);
    dataset = group.createDataSet( "ws", type_double, data_space);
    dataset.write(&ws, type_double);
    if(tdmc){
        group = output.createGroup("/summary/tdmc");
        dataset = group.createDataSet( "split", type_ull, data_space);
        dataset.write(&tdmc_split, type_ull);
        hsize_t dimsv[1]; dimsv[0] = tdmc_time.size();
        H5::DataSpace data_spacev(1,dims);
        dataset = group.createDataSet( "time", type_double, data_spacev);
        dataset.write(tdmc_time.data(), type_ull);
    }


    // Report estimators
    for ( auto& E : Estimators ) { E->report( output ); }
	
    std::cout<<"Output creation done!\n\n";
}
