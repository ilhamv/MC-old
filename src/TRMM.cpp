#include "TRMM.h"


//=============================================================================
// Constructor: Set TRM
//=============================================================================

TRMM::TRMM( const Simulator& S )
{
    // Sizes
    // N_flux is divided by 8 due to the number of scores
    N_flux = S.trmm_estimator_simple->tally_size()
             / S.trmm_estimator_simple->score_size();
    TRM = Eigen::MatrixXd(N_flux,N_flux);

    // Block M
    for( int f = 0; f < N_flux; f++ ){
        for( int i = 0; i < N_flux; i++ ){
            if( i == f ){
                TRM(i,i)  = -S.trmm_estimator_simple->tally(i).mean;
                TRM(i,i) +=  S.trmm_estimator_scatter->tally(i+i*N_flux).mean;
                TRM(i,i) +=  S.trmm_estimator_fission->tally(i+i*N_flux).mean;
                TRM(i,i) /=  S.trmm_estimator_simple->tally(i+N_flux).mean;
            } else{
                TRM(f,i)  = S.trmm_estimator_scatter->tally(f+i*N_flux).mean;
                TRM(f,i) += S.trmm_estimator_fission->tally(f+i*N_flux).mean;
                TRM(f,i) /= S.trmm_estimator_simple->tally(i+N_flux).mean;
            }
        }
    }
}


//=============================================================================
// Solve TRMM
//=============================================================================

void TRMM::solve()
{
    // Solve eigenvalue of TRM
    Eigen::EigenSolver<Eigen::MatrixXd> eSolve(TRM);
    phi_mode = eSolve.eigenvectors();
    alpha    = eSolve.eigenvalues();
    alpha_real.resize(N_flux);
    alpha_imag.resize(N_flux);
    for( int i = 0; i < N_flux; i++ ){
        alpha_real[i] = alpha[i].real();
        alpha_imag[i] = alpha[i].imag();
    }

    // Solve coefficients via initial condition
    Eigen::VectorXcd phi0;
    Eigen::VectorXcd A;
    phi0 = Eigen::VectorXcd::Zero(N_flux);
    phi0(N_flux-1) = 13831.5926439 * std::sqrt( 14.1E6 ) * 100.0;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> dec(phi_mode);
    A = dec.solve(phi0);

    // Construct solution in time
    t = std::vector<double>{0.0, 3E-8, 15E-8, 4E-6, 1E-4};
    phi = Eigen::MatrixXcd::Zero(t.size(),N_flux);
    phi_real.resize(N_flux*t.size());

    unsigned long long idx = 0;
    for( int i = 0; i < t.size(); i++){
        for(int g = 0; g < N_flux; g++){
            for(int n = 0; n < N_flux; n++){
                phi(i,g) += A(n) * phi_mode(g,n) * std::exp( alpha[n] * t[i] );
            }
            phi_real[idx] = phi(i,g).real();
            idx++;
        }
    }
}


//=============================================================================
// Report TRMM
//=============================================================================

void TRMM::report( H5::H5File& output )
{
    // H5 tools
    H5::DataSet dataset;
    H5::Group group;
    H5::DataSpace space_scalar(H5S_SCALAR);
    H5::DataType type_ull    = H5::PredType::NATIVE_ULLONG;
    H5::DataType type_double = H5::PredType::NATIVE_DOUBLE;
    H5::StrType type_str(0, H5T_VARIABLE);

    // TRMM results
    group = output.createGroup("/TRMM");
    hsize_t dimsM[2]; dimsM[0] = N_flux; dimsM[1] = N_flux;
    H5::DataSpace data_spaceM(2,dimsM);
    dataset = group.createDataSet( "TRM", type_double, data_spaceM);
    dataset.write(TRM.data(), type_double);
    hsize_t dims[2]; dims[0] = t.size(); dims[1] = N_flux;
    H5::DataSpace data_spacev(2,dims);
    dataset = group.createDataSet( "flux", type_double, data_spacev);
    dataset.write(phi_real.data(), type_double);
    
    hsize_t dims_alpha[1]; dims_alpha[0] = N_flux;
    H5::DataSpace space_alpha(1,dims_alpha);
    group = group.createGroup("alpha");
    dataset = group.createDataSet( "real", type_double, space_alpha);
    dataset.write(alpha_real.data(), type_double);
    dataset = group.createDataSet( "imag", type_double, space_alpha);
    dataset.write(alpha_imag.data(), type_double);
}
