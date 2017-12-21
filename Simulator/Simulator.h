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
        std::string                                  simName;              // Simulation name
        unsigned long long                           nhist;                // Number of particle samples
        double                                       Ecut_off  = 0.0;      // Energy cut-off
        double                                       tcut_off  = MAX;      // Time cut-off
        unsigned long long                           trackTime = 0;        // "Computation time" (particle track) for variance reduction
        Source_Bank                                  Sbank;                // Source bank
        std::stack  < Particle_t >                   Pbank;                // Particle bank
        std::vector < std::shared_ptr<Surface_t>   > Surface;              // Surfaces
        std::vector < std::shared_ptr<Cell_t>      > Cell;                 // Cells
        std::vector < std::shared_ptr<Nuclide_t>   > Nuclide;              // Nuclides
        std::vector < std::shared_ptr<Material_t>  > Material;             // Materials
        std::vector < std::shared_ptr<Estimator_t> > Estimator;            // Estimators  	
        // User-defined distributions
        std::vector < std::shared_ptr<Distribution_t<double>> > double_distributions;
        std::vector < std::shared_ptr<Distribution_t<int>>    > int_distributions;
        std::vector < std::shared_ptr<Distribution_t<Point_t>>> point_distributions;

        // Constructor: Set up the simulator with XML parser
        Simulator_t( const std::string input_dir );
        ~Simulator_t() {};

        // Start simulation
        void start();
        // Report results
        void report();
};
