#include <vector>  // vector
#include <cstring> // string, strcmp
#include <memory>  // shared_ptr, make_shared
#include <stack>   // stack

#include "VReduction.h"         // Split_Roulette
#include "Const.h"              // MAX
#include "pugixml.hpp"
#include "Geometry.h"
#include "Particle.h"
#include "Distribution.h"
#include "Source.h"
#include "Nuclide.h"
#include "Material.h"
#include "Reaction.h"
#include "Estimator.h"


class Simulator_t
{
    public:
        std::string simulation_name;
        std::string io_dir;
        
        unsigned long long Nsample;     
        unsigned long long Ncycle   = 1; 
        unsigned long long Npassive = 0;
        unsigned long long icycle, isample;

        // Mode flags
        bool ksearch = false; // ksearch mode flag
        bool tally   = false; // Estimator activation toggle for ksearch
        bool tdmc    = false; // TDMC mode flag
        
        // Cut off
        double Ecut_off = 0.0;
        double tcut_off = MAX;

        // Survival rouletting
        double wr = 0.25; // Weight rouletting
        double ws = 1.0;  // Survival weight

        unsigned long long  tracks = 0; // # of particle tracks
        
        Source_Bank            Fbank; // Fission bank
        Source_Bank            Sbank; // Sample bank
        std::stack<Particle_t> Pbank; // Particle bank

        // Universe
        std::vector<std::shared_ptr<Surface_t>>  Surfaces;
        std::vector<std::shared_ptr<Cell>>       Cells;
        std::vector<std::shared_ptr<Nuclide_t>>  Nuclides;
        std::vector<std::shared_ptr<Material_t>> Materials;
       
        // Estimators
        std::vector<std::shared_ptr<Estimator>> Estimators;
        std::shared_ptr<Estimator>              k_C;
        
        // ksearch  specific
        std::shared_ptr<EstimatorK> k_estimator;
        double k = 1.0;

        // TDMC specific
        std::vector<double> tdmc_time;
        int                 tdmc_split = 1.0;
        
        // User-defined distributions
        std::vector<std::shared_ptr<Distribution_t<double>>> Distribution_Double;
        std::vector<std::shared_ptr<Distribution_t<Point_t>>> Distribution_Point;

        // Constructor: Set up the simulator with XML parser
        Simulator_t( const std::string input_dir );
        ~Simulator_t() {};

        // Start simulation
        void start();
        // Report results
        void report();

        void move_particle( Particle_t& P, const double l );
        void cut_off( Particle_t&P );
        void collision( Particle_t& P );
};
