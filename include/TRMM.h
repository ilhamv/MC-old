#ifndef TRMM_H
#define TRMM_H

#include <Eigen/Dense>
#include "H5Cpp.h"

#include "Simulator.h"

class TRMM
{
    private:
        unsigned long long  N_flux;
        unsigned            N_precursor;
        Eigen::MatrixXd     TRM;
        Eigen::MatrixXcd    phi_mode;
        Eigen::VectorXcd    alpha;
        Eigen::MatrixXcd    phi;
        std::vector<double> t;
        
        std::vector<double> alpha_real;
        std::vector<double> alpha_imag;
        std::vector<double> phi_real;

    public:
        TRMM( const Simulator& S );
        ~TRMM() {};

        void solve();
        void report( H5::H5File& output );
};


#endif // TRMM_H
