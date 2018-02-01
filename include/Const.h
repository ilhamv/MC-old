#ifndef _CONST_HEADER_
#define _CONST_HEADER_

#include <cmath>  // acos, sqrt
#include <limits> // epsilon(), max()


const double PI       = std::acos( -1.0 );
const double PI2      = 2.0 * PI;
const double PI_half  = 0.5 * PI;
const double PI_sqrt  = std::sqrt(PI);
const double EPSILON  = std::numeric_limits<float>::epsilon();
const double MAX      = std::numeric_limits<float>::max();
const double MAX_less = 0.9 * MAX; 


#endif

