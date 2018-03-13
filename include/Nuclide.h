#ifndef _NUCLIDE_H
#define _NUCLIDE_H

#include <memory>  
#include <vector>  
#include <cstring> 

#include "Reaction.h"


//==============================================================================
// Nuclide
//==============================================================================

class Nuclide
{
    private:
    	const std::string n_name;
	const double n_A;

        // To avoid redundant table look-up
	unsigned long long idx  = 0;
	double             E_current = 0; 
	
        // Check energy at cross secton call
	//   if it's another different energy, search the location on the table 
        //     --> idx
	void checkE( const double E );
        
        // Reactions and XS
	const std::shared_ptr<Reaction> n_capture;
	const std::shared_ptr<Reaction> n_absorb;
	const std::shared_ptr<Reaction> n_total;
	const std::shared_ptr<ReactionScatter> n_scatter;
	const std::shared_ptr<ReactionFission> n_fission;
        const std::vector<double> n_E;
	

    public:
	Nuclide( const std::string str, const double a, 
                 const std::shared_ptr<Reaction> rc, 
                 const std::shared_ptr<ReactionScatter> rs, 
                 const std::shared_ptr<ReactionFission> rf,
                 const std::shared_ptr<Reaction> ra, 
                 const std::shared_ptr<Reaction> rt, 
                 const std::vector<double> Et):
            n_name(str), n_A(a), n_capture(rc), n_scatter(rs), n_fission(rf),
            n_absorb(ra), n_total(rt), n_E(Et) {};
        ~Nuclide() {};

	std::string name();
        std::shared_ptr<ReactionScatter> scatter();
        std::shared_ptr<ReactionFission> fission();
	double sigmaT( const double E );
	double sigmaS( const double E );
	double sigmaA( const double E );
	double sigmaC( const double E );
	double sigmaF( const double E );
	double nusigmaF( const double E );
	double beta( const double E );
};


#endif // NUCLIDE_H
