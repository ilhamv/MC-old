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


//==============================================================================
// Tally: the estimated mean with its uncertainty
//==============================================================================

class Tally
{
    public:
        double hist    = 0.0; // History sum
	double sum     = 0.0; // Sum over all histories
	double squared = 0.0; // Sum of history sum squared over all histories
		
        double mean  = 0.0;   // Estimated mean
	double uncer = 0.0;   // Uncertainty of the estimated mean
	double uncer_rel;     // Relative uncertainty of estimated mean
	double FOM   = 0.0;   // Figure of merit		
	
	 Tally() {};
	~Tally() {};
};


//==============================================================================
// Filter
//   - Supports: ID(Surface, Cell), Energy, Time
//==============================================================================

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
                                               f_grid(g), f_Nbin(g.size()-1) {};
        ~Filter() {};

        // Get the index and the corresponding track length to be scored
        virtual std::vector<std::pair<int,double>> idx_tl( const Particle_t& P,
                                                           const double l,
                                                           const double id) = 0;
        // Getters
        virtual std::vector<double> grid() final { return f_grid; }
        virtual std::string name() final { return f_name; }
        virtual std::string unit() final { return f_unit; }
};
// ID
class FilterID : public Filter
{
    public:
         FilterID( const std::string n, const std::string u, 
                        const std::vector<double> b ): Filter(n,u,b) {};
        ~FilterID() {};

        // Get the index and the corresponding track length to be scored
        std::vector<std::pair<int,double>> idx_tl( const Particle_t& P,
                                                   const double l,
                                                   const double id);
};
// Energy
class FilterEnergy : public Filter
{
    public:
         FilterEnergy( const std::string n, const std::string u, 
                        const std::vector<double> b ): Filter(n,u,b) {};
        ~FilterEnergy() {};

        // Get the index and the corresponding track length to be scored
        std::vector<std::pair<int,double>> idx_tl( const Particle_t& P,
                                                   const double l,
                                                   const double id);
};
// Time
class FilterTime : public Filter
{
    public:
         FilterTime( const std::string n, const std::string u, 
                        const std::vector<double> b ): Filter(n,u,b) {};
        ~FilterTime() {};

        // Get the index and the corresponding track length to be scored
        std::vector<std::pair<int,double>> idx_tl( const Particle_t& P,
                                                   const double l,
                                                   const double told);
};

// Bin base class 
class Bin_t
{
	public:
		std::vector<std::vector<Tally>>     tally;  // Bin tallies ( indexing --> [bin#][score#] )
		const int                             Nbin;   // # of bins
		std::vector<std::shared_ptr<Score>> scores; // Things to be scored
		const int                             Nscore; // # of scores
		const std::string                     unit;   // result unit, for report
	
	public:
		// Constructor: pass grid points and construct bin and their tallies
		Bin_t( const std::vector<double>& grid, const std::vector<Tally> total_tally, std::vector<std::shared_ptr<Score>>& s, 
				 const std::string str ) : Nbin(grid.size() - 1), unit(str), Nscore(s.size())
		{ 
			// Store scores
			scores = s;
			
			// Each bin is set with the same # of score tallies as the total tally of the estimator
			tally.resize( Nbin, total_tally );
		}
		~Bin_t() {};

		// Score bin
		virtual void score( const Particle_t& P, const std::vector<double>& grid, const double told, const double track = 0.0 ) = 0;
};

// Energy bin (multi score)
class Energy_Bin : public Bin_t
{
	public:
		// Constructor: pass grid points
		 Energy_Bin( const std::vector<double>& grid, const std::vector<Tally> total_tally, std::vector<std::shared_ptr<Score>>& s ) : 
			 Bin_t(grid,total_tally,s,"eV") {};
		~Energy_Bin() {};

		void score( const Particle_t& P, const std::vector<double>& grid, const double told, const double track = 0.0 );
}; 

// Time bin (multi score)
class Time_Bin : public Bin_t
{
	public:
		// Constructor: pass grid points
		 Time_Bin( const std::vector<double>& grid, const std::vector<Tally> total_tally, std::vector<std::shared_ptr<Score>>& s ) : 
			 Bin_t(grid,total_tally,s,"sec") {};
		~Time_Bin() {};

		void score( const Particle_t& P, const std::vector<double>& grid, const double told, const double track = 0.0 );
};



/////////////////
/// Estimator ///
/////////////////

// Estimator base class
class Estimator_t
{
	protected:
		const std::string    e_name;     // Estimator name
		unsigned long long   nhist = 0;  // # of histories estimated
		unsigned long long   ncycle = 0;  // # of histories estimated

	public:
		// Constructor: pass the estimator name and bin grid
		 Estimator_t( const std::string n ) : e_name(n) {};
		~Estimator_t() {};

		// Add thing to be scored
		virtual void addScore( const std::shared_ptr<Score>& S ) = 0;

		// Set bin grid and corresponding tallies
		virtual void setBin( const std::string type, std::vector<double> bin ) = 0;
		
		// Score at events
		virtual void score( const Particle_t& P, const double told, const double track = 0.0  ) = 0;
		
		// Closeout
		virtual void endHistory() = 0;              
		virtual void endCycle( const double tracks ) = 0;              
		
		// Report results
		virtual void report( std::ostringstream& output ) = 0;
		
                std::vector<Tally>                  total_tally; // Total tallies [Nscore]
};


// Generic estimator
////////////////////

class Generic_Estimator : public Estimator_t
{
	public:
		std::vector<std::shared_ptr<Score>> scores;      // Things to be scored
		int                                   Nscore = 0;  // # of scores (things to be scored)
		
		std::vector<double>                   grid;        // Bin grid
		std::shared_ptr<Bin_t>                bin;         // Estimator bin
		int                                   Nbin = 0;    // # of bins

		// Constructor: pass the estimator name
		 Generic_Estimator( const std::string n ) : Estimator_t(n) {};
		~Generic_Estimator() {};

		// Add thing to be scored and push new total tally
		void addScore( const std::shared_ptr<Score>& S );

		// Set bin
		virtual void setBin( const std::string type, const std::vector<double> gr );
		
		// Tally operations
		void tally_endHistory( Tally& T );
                void tally_endCycle( Tally& T, const double tracks );
		void tally_average( Tally& T );
		
		// Score at events
		virtual void score( const Particle_t& P, const double told, const double track = 0.0 );

		// Closeout
                virtual void endHistory();
                virtual void endCycle( const double tracks );

		// Report results
		virtual void report( std::ostringstream& output );
};


// Homogenized MG Constant Generator
////////////////////////////////////
/*
class MGXS_Estimator : public Generic_Estimator
{
	protected:
		std::vector<std::vector<std::shared_ptr<Bin_t>>> tensor_bin; // Legendre components tensor ( N x G x G )
		                                                             //   consists of GxG scattering matrix for each legendre order n=[0,N]
		const unsigned int                               N;          // Legendre order considered
		std::vector<double>                              Pl;         // To store legendre polynomials value		
		// Simple group constants are handled by generic estimator bin
	
	public:
		// Constructor: pass the estimator name and energy grids
		MGXS_Estimator( const std::string n, const unsigned int pn ) : Generic_Estimator(n), N(pn) {};
		~MGXS_Estimator() {};
		
		// Set bin (or group structure)
		void setBin( const std::string type, const std::vector<double> gr );

		// Calculate Legendre polynomials at mu
		void calculatePl( const double mu );

		// Score at events
		void score( const Particle_t& P, const double told, const double track = 0.0 );
		
		// Closeout history
		// Update the sum and sum of squared, and restart history sum of all tallies
		void endHistory();
                void endCycle( const double tracks );

		// Report results
		void report( std::ostringstream& output );
};
*/

#endif // ESTIMATOR_H

