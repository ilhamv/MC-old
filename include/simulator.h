#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>  
#include <cstring> 
#include <memory>  

#include "H5Cpp.h"

#include "Constants.h"              
#include "Distribution.h"
#include "Particle.h"
#include "Estimator.h"
#include "Source.h"
#include "Geometry.h"
#include "Material.h"
#include "Nuclide.h"
#include "Entropy.h"


class Simulator
{
    private:
        // Distributions
        DistributionIsotropicDirection isotropic_direction;
        
        // Banks
        SourceBank            Fbank;  // Fission bank
        SourceBank            TDbank; // Time dependent bank
        SourceBank            Sbank;  // Sample bank
        std::vector<Particle> Pbank;  // Particle bank

        //=====================================================================
        // setup
        //=====================================================================

        void set_nuclide( const std::string name, const std::string label, 
                          std::shared_ptr<Nuclide>& Nuc );
        template<typename T>
        std::shared_ptr<T> find_by_name
                           ( const std::vector<std::shared_ptr<T>>& vec, 
                             const std::string name );

        //=====================================================================
        // general functions
        //=====================================================================
 
        void random_walk( Particle& P );
        bool test_point( const Point& p, const std::shared_ptr<Cell>& C );
        std::shared_ptr<Cell> search_cell( const Point& p );
        double collision_distance( const Particle& P );
	std::pair<std::shared_ptr<Surface>, double> 
                                     surface_intersect( const Particle& P );
        void move_particle( Particle& P, const double l );
        void surface_hit( Particle& P, const std::shared_ptr<Surface>& S );
        void push_particle_bank( const Particle& P );
        void collision( Particle& P );
        
        //=====================================================================
        // ksearch functions
        //=====================================================================
        
        void add_fission_source( const std::shared_ptr<Source>& S, 
                                 const double i );
        void implicit_fission_ksearch( const Particle& P, const double bank_nu, 
                                const std::shared_ptr<Nuclide>& N_fission );

        //=====================================================================
        // fixed source functions
        //=====================================================================
        
        void implicit_fission_fixed_source( const Particle& P, 
            const double bank_nu, const std::shared_ptr<Nuclide>& N_fission );

        //=====================================================================
        // time dependent functions
        //=====================================================================

        void time_hit( Particle& P );
        Particle forced_decay( const Particle& P, 
                               const std::shared_ptr<Nuclide>& N, 
                               const double initial, const double interval,
                               const int p_tdmc);

        //=====================================================================
        // Population control
        //=====================================================================

        // Weight roulette
        double wr = 0.001; // Weight rouletting
        double ws = 1.0;  // Survival weight

        // Particle comb
        bool comb = false;
        int comb_teeth, bank_max;
        
        void weight_roulette( Particle& P );
        void cell_importance( Particle& P, std::vector<Particle>& Pbank );
        void particle_comb( std::vector<Particle>& Pbank );


    public:
        std::string simulation_name;
        
        unsigned long long Nsample;     
        unsigned long long Ncycle   = 1; 
        unsigned long long Npassive = 0;
        unsigned long long icycle, isample;
        unsigned long long Ntrack = 0;

        // Mode flags
        std::string mode = "fixed source";
        bool ksearch     = false; 
        bool tdmc        = false;
        bool trmm        = false;
        
        // THE OBJECTS
        std::vector<std::shared_ptr<Surface>>  Surfaces;
        std::vector<std::shared_ptr<Cell>>       Cells;
        std::vector<std::shared_ptr<Nuclide>>  Nuclides;
        std::vector<std::shared_ptr<Material>> Materials;
        std::vector<std::shared_ptr<Estimator>>  Estimators;
        
        // ksearch specific
        std::shared_ptr<EstimatorK> k_estimator;
        double k = 1.0;
        bool tally = false; // Estimator activation toggle

        // TDMC specific
        std::vector<double> tdmc_time;
        std::vector<double> tdmc_interval;

        // TRMM specific
        std::shared_ptr<Estimator> trmm_estimator_simple;
        std::shared_ptr<Estimator> trmm_estimator_scatter;
        std::shared_ptr<Estimator> trmm_estimator_fission_prompt;
        std::vector<std::shared_ptr<Estimator>> trmm_estimator_fission_delayed;
        
        // User-defined distributions
        std::vector<std::shared_ptr<Distribution<double>>>Distribution_Double;
        std::vector<std::shared_ptr<Distribution<Point>>>Distribution_Point;

        // Constructor: Set up the simulator with XML parser
        Simulator( const std::string input_dir );
        ~Simulator() {};

        // Start simulation
        void start();
        void report( const std::string io_dr );
};

#endif // SIMULATOR_H
