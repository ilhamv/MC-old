#include <iostream>
#include <Eigen/Dense>
#include <cstring>
#include <cmath>

#include "H5Cpp.h"
#include "eigen3-hdf5.hpp"


int main( int argc, char* argv[] )
{
    // I/O directory and TRM HDF5 file name
    std::string file_name = std::string(argv[1]);
    size_t last = file_name.find_last_of('/');
    std::string io_dir = file_name.substr(0,last+1);

    //=========================================================================
    // Load TRM
    //=========================================================================

    Eigen::MatrixXd TRM;
    
    const H5std_string FILE_NAME( file_name );
    H5::H5File file( FILE_NAME, H5F_ACC_RDONLY );
    EigenHDF5::load(file, "TRM", TRM);
    const int N_TRM = std::sqrt(TRM.size());
    const int J = 6;
    const int G = N_TRM - J;
    
    //=========================================================================
    // Solve eigen-pairs of TRM
    //=========================================================================
    
    Eigen::MatrixXcd phi_mode;
    Eigen::VectorXcd alpha;

    Eigen::EigenSolver<Eigen::MatrixXd> eSolve(TRM);
    phi_mode = eSolve.eigenvectors();
    alpha    = eSolve.eigenvalues();

    //=========================================================================
    // Set TRM -> Adjoint, then solve the eigen-pairs
    //=========================================================================

    // Set the adjoint TRM
    Eigen::VectorXd speed_inv;
    EigenHDF5::load(file, "inverse_speed", speed_inv);
    for( int i = 0; i < G; i++ ){
        for( int j = 0; j < N_TRM; j++ ){
            TRM(i,j) *= speed_inv[i];
        }
    }
    TRM.transposeInPlace();
    for( int i = 0; i < G; i++ ){
        for( int j = 0; j < N_TRM; j++ ){
            TRM(i,j) /= speed_inv[i];
        }
    }
    
    // Solve the eigen-pairs
    Eigen::MatrixXcd phi_mode_adj;
    Eigen::VectorXcd alpha_adj;

    eSolve = Eigen::EigenSolver<Eigen::MatrixXd>(TRM);
    phi_mode_adj = eSolve.eigenvectors();
    alpha_adj    = eSolve.eigenvalues();

    //=========================================================================
    // Output eigen-pairs
    //=========================================================================

    const H5std_string FILE_NAME_OUT( io_dir+"output_TRMM.h5" );
    H5::H5File output( FILE_NAME_OUT, H5F_ACC_TRUNC );

    EigenHDF5::save(output, "alpha", alpha);
    EigenHDF5::save(output, "alpha_adj", alpha_adj);
    EigenHDF5::save(output, "phi_mode", phi_mode);
    EigenHDF5::save(output, "phi_mode_adj", phi_mode_adj);

    return 0;
}
