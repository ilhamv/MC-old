#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>  
#include <limits> 


const double PI             = std::acos( -1.0 );
const double PI_2           = 2.0 * PI;
const double PI_half        = 0.5 * PI;
const double PI_sqrt        = std::sqrt(PI);
const double EPSILON_float  = std::numeric_limits<float>::epsilon();
const double MAX_float      = std::numeric_limits<float>::max();
const double MAX_float_less = 0.9 * MAX_float; 


#endif // CONSTANTS_H
