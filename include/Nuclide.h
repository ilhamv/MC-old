#ifndef _NUCLIDE_H
#define _NUCLIDE_H

#include <memory>  // shared_ptr
#include <vector>  // vector
#include <cstring> // strcmp

#include "Reaction.h"

class Reaction_t;

//==============================================================================
// Nuclide
//==============================================================================

class Nuclide
{
    private:
    	const std::string n_name;
	const double n_A;
	std::vector<std::shared_ptr<Reaction_t>> n_reactions;

        // To avoid redundant table look-up
	unsigned long long idx_help  = 0;
	double             E_current = 0; 
	
        // xs table energy grid
        std::shared_ptr<std::vector<double>> E_table = 
            std::make_shared<std::vector<double>>();
	
        // Check energy at cross secton call
	//   if it's another different energy, search the location on the table --> idx_help
	void checkE( const double E );
	

    public:
	Nuclide( const std::string str, const double a ) : n_name(str), 
                                                             n_A(a) {};
        ~Nuclide() {};

	// Getters
	std::string name();
	double sigmaT( const double E );
	double sigmaS( const double E );
	double sigmaA( const double E );
	double sigmaC( const double E );
	double sigmaF( const double E );
	double nusigmaF( const double E );

	// Set energy grids for table look-up XS
	void setTable ( const std::shared_ptr< std::vector<double> >& Evec );

	// Sample Chi spectrum
	double Chi( const double E );

	// Add reaction
	void addReaction( const std::shared_ptr< Reaction_t >& C );
	
	// Randomly sample a reaction type from the nuclide
	std::shared_ptr<Reaction_t> reaction_sample( const double E, const bool ksearch );
        
        // Reaction
	std::shared_ptr<Reaction_t> scatter = nullptr;
};


#endif // NUCLIDE_H
