#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

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

// Forward declaration
class Particle_t;

///////////////
/// Scoring ///
///////////////

// Score base class
class Score_t
{
	private:
		const std::string s_name;
	public:
		 Score_t( const std::string n ) : s_name(n) {};
		~Score_t() {};

		// Get score type name
		virtual std::string name() final { return s_name; }

		// Get score to be added at event
		virtual double score( const Particle_t& P, const double l ) = 0;
};

// Current (event)
class Current_Score : public Score_t
{
	public:
		 Current_Score() : Score_t( "Current" ) {};
		~Current_Score() {};

		double score( const Particle_t& P, const double l );
};

// Flux (path length)
class Flux_Score : public Score_t
{
	public:
		 Flux_Score() : Score_t( "Flux" ) {};
		~Flux_Score() {};

		double score( const Particle_t& P, const double l );
};

// Absorption (path length)
class Absorption_Score : public Score_t
{
	public:
		 Absorption_Score() : Score_t( "Abs. Rate" ) {};
		~Absorption_Score() {};

		double score( const Particle_t& P, const double l );
};

// Scatter (path length)
class Scatter_Score : public Score_t
{
	public:
		 Scatter_Score() : Score_t( "Scat. Rate" ) {};
		~Scatter_Score() {};

		double score( const Particle_t& P, const double l );
};

// Capture (path length)
class Capture_Score : public Score_t
{
	public:
		 Capture_Score() : Score_t( "Capt. Rate" ) {};
		~Capture_Score() {};

		double score( const Particle_t& P, const double l );
};

// Fission (path length)
class Fission_Score : public Score_t
{
	public:
		 Fission_Score() : Score_t( "Fis. Rate" ) {};
		~Fission_Score() {};

		double score( const Particle_t& P, const double l );
};


// Nu Fission (path length)
class nuFission_Score : public Score_t
{
	public:
		 nuFission_Score() : Score_t( "Prod. Rate" ) {};
		~nuFission_Score() {};

		double score( const Particle_t& P, const double l );
};


// Total (path length)
class Total_Score : public Score_t
{
	public:
		 Total_Score() : Score_t( "Tot. Rate" ) {};
		~Total_Score() {};

		double score( const Particle_t& P, const double l );
};



/////////////
/// Tally ///
/////////////

class Tally_t
{
	public:
		double hist    = 0.0; // History sum
		double sum     = 0.0; // Sum over all histories
	       	double squared = 0.0; // Sum of history sum squared over all histories
		
                double mean  = 0.0;   // Estimated mean
		double uncer = 0.0;   // Uncertainty of the estimated mean
		double uncer_rel;     // Relative uncertainty of estimated mean
		double FOM   = 0.0;   // Figure of merit		
	
	public:
		 Tally_t() {};
		~Tally_t() {};
};



///////////
/// Bin ///
///////////

// Bin base class 
class Bin_t
{
	public:
		std::vector<std::vector<Tally_t>>     tally;  // Bin tallies ( indexing --> [bin#][score#] )
		const int                             Nbin;   // # of bins
		std::vector<std::shared_ptr<Score_t>> scores; // Things to be scored
		const int                             Nscore; // # of scores
		const std::string                     unit;   // result unit, for report
	
	public:
		// Constructor: pass grid points and construct bin and their tallies
		Bin_t( const std::vector<double>& grid, const std::vector<Tally_t> total_tally, std::vector<std::shared_ptr<Score_t>>& s, 
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
		 Energy_Bin( const std::vector<double>& grid, const std::vector<Tally_t> total_tally, std::vector<std::shared_ptr<Score_t>>& s ) : 
			 Bin_t(grid,total_tally,s,"eV") {};
		~Energy_Bin() {};

		void score( const Particle_t& P, const std::vector<double>& grid, const double told, const double track = 0.0 );
}; 

// Time bin (multi score)
class Time_Bin : public Bin_t
{
	public:
		// Constructor: pass grid points
		 Time_Bin( const std::vector<double>& grid, const std::vector<Tally_t> total_tally, std::vector<std::shared_ptr<Score_t>>& s ) : 
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
		virtual void addScore( const std::shared_ptr<Score_t>& S ) = 0;

		// Set bin grid and corresponding tallies
		virtual void setBin( const std::string type, std::vector<double> bin ) = 0;
		
		// Score at events
		virtual void score( const Particle_t& P, const double told, const double track = 0.0  ) = 0;
		
		// Closeout
		virtual void endHistory() = 0;              
		virtual void endCycle( const double tracks ) = 0;              
		
		// Report results
		virtual void report( std::ostringstream& output ) = 0;
		
                std::vector<Tally_t>                  total_tally; // Total tallies [Nscore]
};


// Generic estimator
////////////////////

class Generic_Estimator : public Estimator_t
{
	public:
		std::vector<std::shared_ptr<Score_t>> scores;      // Things to be scored
		int                                   Nscore = 0;  // # of scores (things to be scored)
		
		std::vector<double>                   grid;        // Bin grid
		std::shared_ptr<Bin_t>                bin;         // Estimator bin
		int                                   Nbin = 0;    // # of bins

		// Constructor: pass the estimator name
		 Generic_Estimator( const std::string n ) : Estimator_t(n) {};
		~Generic_Estimator() {};

		// Add thing to be scored and push new total tally
		void addScore( const std::shared_ptr<Score_t>& S );

		// Set bin
		virtual void setBin( const std::string type, const std::vector<double> gr );
		
		// Tally operations
		void tally_endHistory( Tally_t& T );
                void tally_endCycle( Tally_t& T, const double tracks );
		void tally_average( Tally_t& T );
		
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


#endif

