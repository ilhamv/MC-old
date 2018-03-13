#include <vector>  
#include <cstring> 
#include <memory>  
#include <stack>   

#include <Eigen/Dense>

#include "Constants.h"              
#include "Distribution.h"
#include "Particle.h"
#include "Estimator.h"
#include "Source.h"
#include "Geometry.h"
#include "Material.h"
#include "Nuclide.h"


class Simulator
{
    private:
        void set_nuclide( const std::string name, const std::string label, 
                          std::shared_ptr<Nuclide>& Nuc );
        template<typename T>
        std::shared_ptr<T> find_by_name
                           ( const std::vector<std::shared_ptr<T>>& vec, 
                             const std::string name );

        bool test_point( const Point& p, const std::shared_ptr<Cell>& C );
        std::shared_ptr<Cell> search_cell( const Point& p );
        void move_particle( Particle& P, const double l );
        void collision( Particle& P );
        void surface_hit( Particle& P, const std::shared_ptr<Surface>& S );
        void cut_off( Particle& P );
        double collision_distance( const Particle& P );
	std::pair<std::shared_ptr<Surface>, double> 
            surface_intersect( const Particle& P );

        void add_fission_source( const std::shared_ptr<Source>& S, 
                                 const double i );
        void push_particle_bank( const Particle& P );

        DistributionIsotropicDirection isotropic;

    public:
        std::string simulation_name;
        std::string io_dir;
        
        unsigned long long Nsample;     
        unsigned long long Ncycle   = 1; 
        unsigned long long Npassive = 0;
        unsigned long long icycle, isample;
        unsigned long long Ntrack = 0;

        // Mode flags
        std::string mode = "fixed source";
        bool ksearch = false; 
        bool tdmc    = false;
        bool trmm    = false;
        
        // Cut off
        double Ecut_off = 0.0;
        double tcut_off = MAX_float;

        // Survival rouletting
        double wr = 0.25; // Weight rouletting
        double ws = 1.0;  // Survival weight
        
        // Banks
        SourceBank           Fbank; // Fission bank
        SourceBank           Sbank; // Sample bank
        std::stack<Particle> Pbank; // Particle bank

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
        unsigned long long  tdmc_split = 1.0;

        // TRMM specific
        std::shared_ptr<Estimator> trmm_estimator_collision;
        std::shared_ptr<Estimator> trmm_estimator_scatter;
        std::shared_ptr<Estimator> trmm_estimator_fission;
        
        // User-defined distributions
        std::vector<std::shared_ptr<Distribution<double>>>Distribution_Double;
        std::vector<std::shared_ptr<Distribution<Point>>>Distribution_Point;

        // Constructor: Set up the simulator with XML parser
        Simulator( const std::string input_dir );
        ~Simulator() {};

        // Start simulation
        void start();
        void report();
};
