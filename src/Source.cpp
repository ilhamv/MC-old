#include <cmath>
#include <iostream>

#include "Distribution.h"
#include "Random.h"
#include "Particle.h"
#include "Source.h"
#include "Point.h"


//=============================================================================
// Particle Source
//=============================================================================

Particle SourceDelta::get_source()
{
    return s_P;
}
Particle SourcePoint::get_source()
{
    return Particle( s_pos, s_dir->sample(), s_energy->sample(),
                     0.0, 1.0, 0, s_cell );
}


//=============================================================================
// Source Bank
//=============================================================================

void SourceBank::add_source( const std::shared_ptr<Source> S, 
                            const double prob )
{
    sources.push_back( std::make_pair( S, prob ) );
    total += prob;
}
void SourceBank::reset()
{
    sources.clear();
    total = 0.0;
}
Particle SourceBank::get_source()
{
    const double xi = total * Urand();
    double s  = 0.0;
    for( const auto& So : sources ) {
        s += So.second;
        if ( s > xi ) { 
            return So.first->get_source();
        }
    }
    if( total == 0 ){
        std::cout<< "[ERROR] Source bank is empty...\n";
        std::exit(EXIT_FAILURE);
    }
    else{
        std::cout<< "[ERROR] Source particle lost...\n";
        std::exit(EXIT_FAILURE);
    }
}
