#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include <cmath>    // sqrt
#include <iostream> // cout
#include <iomanip>
#include <vector>   // vector
#include <fstream>  // ofstream
#include <memory>   // shared_ptr
#include <cstring>  // string
#include <sstream>  // ostringstream

#include "Particle.h"
#include "Geometry.h"
#include "Distribution.h"
#include "H5Cpp.h"

class Particle_t;


//=============================================================================
// Scoring Kernel
//   - Supports: Neutron, Track Length, Collision
//=============================================================================

class ScoreKernel
{
    public:
         ScoreKernel() {};
        ~ScoreKernel() {};

        virtual double score( const Particle_t& P, const double l ) = 0;
};
// Neutron
class ScoreKernelNeutron : public ScoreKernel
{
    public:
         ScoreKernelNeutron() {};
        ~ScoreKernelNeutron() {};

        double score( const Particle_t& P, const double l );
};
// Track Length
class ScoreKernelTrackLength : public ScoreKernel
{
    public:
         ScoreKernelTrackLength() {};
        ~ScoreKernelTrackLength() {};

        double score( const Particle_t& P, const double l );
};
// Collision
class ScoreKernelCollision : public ScoreKernel
{
    public:
         ScoreKernelCollision() {};
        ~ScoreKernelCollision() {};

        double score( const Particle_t& P, const double l );
};
// Velocity
class ScoreKernelVelocity : public ScoreKernel
{
    public:
         ScoreKernelVelocity() {};
        ~ScoreKernelVelocity() {};

        double score( const Particle_t& P, const double l );
};
// Track Length Velocity
class ScoreKernelTrackLengthVelocity : public ScoreKernel
{
    public:
         ScoreKernelTrackLengthVelocity() {};
        ~ScoreKernelTrackLengthVelocity() {};

        double score( const Particle_t& P, const double l );
};


//=============================================================================
// Score
//   - Specified Scoring Kernel
//   - Supports: Flux, Absorption, Scatter, Capture, Fission,
//               NuFission, Total
//=============================================================================

class Score
{
    private:
        const std::string s_name;

    protected:
        std::shared_ptr<ScoreKernel> s_kernel;

    public:
    	Score( const std::string n, const std::shared_ptr<ScoreKernel>& k )
            : s_name(n), s_kernel(k) {};
        ~Score() {};

	// Get score name
        virtual std::string name() final {return s_name;}

        // Get score to be added at event
	virtual double score( const Particle_t& P, const double l ) = 0;
};

// Flux
class ScoreFlux : public Score
{
    public:
	ScoreFlux( const std::string n,
                    const std::shared_ptr<ScoreKernel>& k )
            : Score(n,k) {};
	~ScoreFlux() {};

	double score( const Particle_t& P, const double l );
};
// Absorption
class ScoreAbsorption : public Score
{
    public:
	ScoreAbsorption( const std::string n,
                         const std::shared_ptr<ScoreKernel>& k )
            : Score(n,k) {};
	~ScoreAbsorption() {};

	double score( const Particle_t& P, const double l );
};
// Scatter
class ScoreScatter : public Score
{
    public:
	ScoreScatter( const std::string n,
                      const std::shared_ptr<ScoreKernel>& k )
            : Score(n,k) {};
	~ScoreScatter() {};

	double score( const Particle_t& P, const double l );
};
// Capture
class ScoreCapture : public Score
{
    public:
	ScoreCapture( const std::string n,
                      const std::shared_ptr<ScoreKernel>& k )
            : Score(n,k) {};
	~ScoreCapture() {};

	double score( const Particle_t& P, const double l );
};
// Fission
class ScoreFission : public Score
{
    public:
	ScoreFission( const std::string n,
                      const std::shared_ptr<ScoreKernel>& k )
            : Score(n,k) {};
	~ScoreFission() {};

	double score( const Particle_t& P, const double l );
};
// NuFission
class ScoreNuFission : public Score
{
    public:
	ScoreNuFission( const std::string n,
                        const std::shared_ptr<ScoreKernel>& k )
            : Score(n,k) {};
	~ScoreNuFission() {};

	double score( const Particle_t& P, const double l );
};
// Total
class ScoreTotal : public Score
{
    public:
	ScoreTotal( const std::string n,
                         const std::shared_ptr<ScoreKernel>& k )
            : Score(n,k) {};
	~ScoreTotal() {};

	double score( const Particle_t& P, const double l );
};


//=============================================================================
// Tally: the estimated mean with its uncertainty
//=============================================================================

class Tally
{
    public:
        double hist    = 0.0; // History sum
	double sum     = 0.0; // Sum over all histories
	double squared = 0.0; // Sum of history sum squared over all histories
		
        double mean  = 0.0;   // Estimated mean
	double uncer = 0.0;   // Uncertainty of the estimated mean
	
	 Tally() {};
	~Tally() {};
};


//=============================================================================
// Filter
//   - Supports: recently crossed Surface, Cell, Energy, Time
//=============================================================================

class Filter
{
    private:
        const std::string f_name;
        const std::string f_unit;

    protected:
        const std::vector<double> f_grid;
        const int                 f_Nbin;

    public:
         Filter( const std::string n, const std::string u, 
                 const std::vector<double> g): f_name(n), f_unit(u), 
                                               f_grid(g), f_Nbin(g.size()-1){};
        ~Filter() {};

        // Get the index and the corresponding track length to be scored
        virtual std::vector<std::pair<int,double>> idx_l( Particle_t& P,
                                                          const double l ) = 0;
        // Getters
        virtual int                 size() = 0;
        virtual std::vector<double> grid() final { return f_grid; }
        virtual std::string         name() final { return f_name; }
        virtual std::string         unit() final { return f_unit; }
};
// Surface
class FilterSurface : public Filter
{
    public:
         FilterSurface( const std::vector<double> g ): 
             Filter("surface","id#",g) {};
        ~FilterSurface() {};

        // Get the index and the corresponding track length to be scored
        std::vector<std::pair<int,double>> idx_l( Particle_t& P,
                                                  const double l );
        virtual int size() final { return f_grid.size(); }
};
// Cell
class FilterCell : public Filter
{
    public:
         FilterCell( const std::vector<double> g ): 
             Filter("cell","id#",g) {};
        ~FilterCell() {};

        // Get the index and the corresponding track length to be scored
        std::vector<std::pair<int,double>> idx_l( Particle_t& P,
                                                  const double l );
        virtual int size() final { return f_grid.size(); }
};
// Energy
class FilterEnergy : public Filter
{
    public:
         FilterEnergy( const std::vector<double> g ): 
             Filter("energy","eV",g) {};
        ~FilterEnergy() {};

        // Get the index and the corresponding track length to be scored
        std::vector<std::pair<int,double>> idx_l( Particle_t& P,
                                                   const double l );
        virtual int size() final { return f_grid.size()-1; }
};
// Time bin
class FilterTime : public Filter
{
    public:
         FilterTime( const std::vector<double> g ): 
             Filter("time","s",g) {};
        ~FilterTime() {};

        // Get the index and the corresponding track length to be scored
        std::vector<std::pair<int,double>> idx_l( Particle_t& P,
                                                   const double l );
        virtual int size() final { return f_grid.size()-1; }
};
// Time - TDMC
class FilterTDMC : public Filter
{
    public:
         FilterTDMC( const std::vector<double> g ): 
             Filter("time","s",g) {};
        ~FilterTDMC() {};

        // Get the index and the corresponding track length to be scored
        std::vector<std::pair<int,double>> idx_l( Particle_t& P,
                                                   const double l );
        virtual int size() final { return f_grid.size(); }
};


//=============================================================================
/// Estimator
//=============================================================================

// Basic Estimator
class Estimator
{
    private:
	const std::string                    e_name;  
        std::vector<std::shared_ptr<Score>>  e_scores;
        std::vector<std::shared_ptr<Filter>> e_filters;
        std::vector<Tally>                   e_tally;
        // Tally structure: [score][filter1][filter2]... casted into 1D
        //   filter1 are the geometries that the estimator is attached to
        
        // Indexes and track length to be scored on each individual filter
        std::vector<std::vector<std::pair<int,double>>> e_idx_l;

        // Multiplication of filter size with index > i
        std::vector<double> idx_factor;

        const unsigned long long e_Nsample, e_Nactive;

    public:
	 Estimator( const std::string n, const unsigned long long Ns,
                    const unsigned long long Na ) : e_name(n), e_Nsample(Ns),
                                                    e_Nactive(Na) {};
	~Estimator() {};

        // Initialization
	void add_score( const std::shared_ptr<Score>& S );
	void add_filter( const std::shared_ptr<Filter>& F );
        void initialize_tallies();
        
        void score( Particle_t& P, const double l );	

	// Loop closeouts
	void end_history();              
	void end_cycle( const double tracks );
        void end_simulation();
	void report( std::ostringstream& output, H5::H5File& output_H5 );

        Tally tally( const int i );
};

// k-eigenvalue
class EstimatorK
{
    private:
        double              k_sum_C  = 0.0;
        double              k_sum_TL = 0.0;
        double              k_sq_C   = 0.0;
        double              k_sq_TL  = 0.0;
        std::vector<double> k_cycle, k_avg, k_uncer;
        unsigned long long  Navg   = 0;
        unsigned long long  icycle = 0;
        double              k_C  = 0.0;
        double              k_TL = 0.0;
        unsigned long long Nsample, Nactive;

        // Accumulator
        double uncer_sq_accumulator = 0.0;
        double mean_accumulator  = 0.0;
    
    public:
        EstimatorK( const unsigned long long Ncycle, 
                    const unsigned long long Na,
                    const unsigned long long Ns );
        ~EstimatorK() {};

        void end_history();
        void report_cycle( const bool tally );
        //void report();
        void estimate_C( const Particle_t& P );
        void estimate_TL( const Particle_t& P, const double l );

        double k;
};

#endif // ESTIMATOR_H
