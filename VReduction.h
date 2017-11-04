#ifndef _SOLVER_HEADER_
#define _SOLVER_HEADER_

#include <cmath>      // floor
#include <memory>     // shared_ptr
#include <stack>      // stack
#include <iostream>

#include "Random.h"
#include "Particle.h"
#include "Geometry.h"


// Perform Splitting & Rouletting variance reduction technique
// arguments: Old working cell, working particle, particle bank
void Split_Roulette( std::shared_ptr<Cell_t>& C, Particle_t& P, std::stack<Particle_t>& Pbank );


#endif
