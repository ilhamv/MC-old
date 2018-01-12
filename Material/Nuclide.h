#ifndef _NUCLIDE_HEADER_
#define _NUCLIDE_HEADER_

#include <memory>  // shared_ptr
#include <vector>  // vector
#include <cstring> // strcmp

#include "Reaction.h"

// Forward declaration
class Reaction_t;

// Nuclide base class
class Nuclide_t
{
	private:
    		std::string                                  n_name;        // Name
		const double                                 n_A;           // Nuclide mass
		std::vector< std::shared_ptr< Reaction_t > > reactions;     // Reactions
		std::shared_ptr< Reaction_t >                scatter;       // Scattering reaction (for simulating in MGXS)
		unsigned long long                           idx_help  = 0; // Helping index to fasten table look-up
		double                                       E_current = 0; 
		std::shared_ptr< std::vector<double> >       E_table = std::make_shared< std::vector<double> >(); // Energy grids on XS table

	public:
		 Nuclide_t( const std::string str, const double a ) : n_name(str), n_A(a) {};
           	~Nuclide_t() {};

		// Getters
		std::string name(); // Name
		// microXs
		double sigmaT( const double E );
		double sigmaS( const double E );
		double sigmaA( const double E );
		double sigmaC( const double E );
		double sigmaF( const double E );
		double nusigmaF( const double E );

		// Set energy grids for table look-up XS
		void setTable ( const std::shared_ptr< std::vector<double> >& Evec );

		// Check energy at cross secton call
		//   if it's another different energy, search the location on the table --> idx_help
		void checkE( const double E );
		
		// Sample Chi spectrum
    		double Chi( const double E );

		// Add reaction
		void addReaction( const std::shared_ptr< Reaction_t >& C );
		
		// Randomly sample a reaction type from the nuclide
		std::shared_ptr<Reaction_t> reaction_sample( const double E, const bool ksearch );
		
		// Simulate scattering for scattering matrix MGXS
		void simulate_scatter( Particle_t& P );
};


#endif
