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
        std::string                                  simName;          // Simulation name
        
        unsigned long long                           nSample;          // Number of particle samples per cycle
        unsigned long long                           nCycle   = 1;     // Number of iteration/cycle
        unsigned long long                           nPassive = 0;     // Number of passive cycle

        bool                                         ksearch  = false; // ksearch mode flag
        bool                                         tally    = false; // Estimator activation toggle
        
        std::vector<double>                          k_cycle;          // Criticality estimate at each cycle
        double                                       k = 1.0;
        std::vector<double>                          entropy_cycle;    // Entropy at each cycle
        
        double                                       Ecut_off = 0.0;   // Energy cut-off
        double                                       tcut_off = MAX;   // Time cut-off

        unsigned long long                           tracks   = 0;     // # of particle tracks
        double                                       wr       = 0.25;  // Weight rouletting
        double                                       ws       = 1.0;  // Survival weight
        
        Source_Bank                                  Fbank;            // Fission bank
        Source_Bank                                  Sbank;            // Sample bank
        std::stack  < Particle_t >                   Pbank;            // Particle bank

        std::vector < std::shared_ptr<Surface_t>   > Surface;          // Surfaces
        std::vector < std::shared_ptr<Cell_t>      > Cell;             // Cells
        std::vector < std::shared_ptr<Nuclide_t>   > Nuclide;          // Nuclides
        std::vector < std::shared_ptr<Material_t>  > Material;         // Materials
        std::vector < std::shared_ptr<Estimator_t> > Estimator;        // Estimators  	
        
        // User-defined distributions
        std::vector < std::shared_ptr<Distribution_t<double>> > Distribution_Double;
        std::vector < std::shared_ptr<Distribution_t<Point_t>>> Distribution_Point;

        // Constructor: Set up the simulator with XML parser
        Simulator_t( const std::string input_dir );
        ~Simulator_t() {};

        // Start simulation
        void start();
        // Report results
        void report();
};
