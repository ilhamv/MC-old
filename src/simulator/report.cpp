#include <vector>  

#include "simulator.h"

//=============================================================================
// Report results
//=============================================================================

void Simulator::report( const std::string io_dir )
{
    // H5 tools
    H5::DataSet dataset;
    H5::Group group;
    H5::DataSpace space_scalar(H5S_SCALAR);
    H5::DataType type_ull    = H5::PredType::NATIVE_ULLONG;
    H5::DataType type_double = H5::PredType::NATIVE_DOUBLE;
    H5::StrType type_str(0, H5T_VARIABLE);
    
    // HDF5 output
    H5std_string output_name( io_dir + "output.h5" );
    H5::H5File output(output_name, H5F_ACC_TRUNC);
    
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

        // Initial precursor
        std::vector<double> C_init(J,0.0);
        for( int j = 0; j < J; j++ ){
            for( int g = 0; g < G; g++ ){
                C_init[j] += trmm_estimator_simple->tally(2*G+J*G+j*G+g).mean;
            }
        }
        dimsv[0] = J;
        data_spacev = H5::DataSpace(1,dimsv);
        dataset = output.createDataSet("C_initial",type_double,data_spacev);
        dataset.write(C_init.data(), type_double);

        // Initial flux
        std::vector<double> psi_init(G);
        for( int g = 0; g < G; g++ ){
            psi_init[g] = trmm_estimator_simple->tally(g+G).mean;
        }
        dimsv[0] = G;
        data_spacev = H5::DataSpace(1,dimsv);
        dataset = output.createDataSet("psi_initial",type_double,data_spacev);
        dataset.write(psi_init.data(), type_double);
    }
}
