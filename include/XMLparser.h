#ifndef _XMLPARSER_HEADER_
#define _XMLPARSER_HEADER_

#include <vector>       // vector
#include <iostream>     // cout
#include <cstring>      // strcmp
#include <memory>       // shared_ptr, make_shared
#include <stack>        // stack
#include <cmath>        // exp
#include <sstream>      // istringstream
#include <fstream>      // file stream
#include <dirent.h>     // open a folder

#include "VReduction.h" // Split_Roulette
#include "Const.h"      // MAX
#include "pugixml.hpp"
#include "Geometry.h"
#include "Particle.h"
#include "Distribution.h"
#include "Source.h"
#include "Nuclide.h"
#include "Material.h"
#include "Reaction.h"
#include "Estimator.h"
#include "XSec.h"


// Function that returns an item from a vector of objects of type T by name provided
// the object class has a string and a method called name() allowing for it to be returned
template< typename T >
std::shared_ptr< T > findByName( const std::vector< std::shared_ptr< T > >& vec, const std::string name );


// Function that set nuclide based on library availability
// Current approximation: isotropic in C.O.M.
void setNuclide( const std::string name, const std::string label, std::shared_ptr<Nuclide_t>& Nuc );


// XML input pasrese
void XML_input
( 
    std::string                                              file_dir,
    std::string&                                             simName,
    unsigned long long&                                      nSample,          
    bool&                                                    ksearch,
    bool&                                                    tdmc,
    unsigned long long&                                      nCycle,          
    unsigned long long&                                      nPassive,          
    double&                                                  Ecut_off,
    double&                                                  tcut_off,
    Source_Bank&                                             Sbank,
    std::vector < std::shared_ptr<Surface_t>   >&            Surface,     
    std::vector < std::shared_ptr<Cell_t>    >&              Cell,    
    std::vector < std::shared_ptr<Nuclide_t>   >&            Nuclide,   
    std::vector < std::shared_ptr<Material_t>  >&            Material, 
    std::vector < std::shared_ptr<Estimator> >&            estimator,
    std::vector < std::shared_ptr<Distribution_t<double>> >& Distribution_Double,
    std::vector < std::shared_ptr<Distribution_t<Point_t>>>& Distribution_Point,
    std::vector<double>& tdmc_time,
    int& tdmc_split
);


#endif
