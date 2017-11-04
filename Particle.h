#ifndef _PARTICLE_HEADER_
#define _PARTICLE_HEADER_

#include <memory> // shared_ptr
#include <vector> // vector

#include "Point.h"


// Forward declarations
class Cell_t;


// Particle base class
class Particle_t
{
    private:
        Point_t                 p_pos;          // Position
        Point_t                 p_dir;          // Direction
        bool                    p_alive = true; // Alive/dead flag
        double                  p_weight;       // Weight
        double                  p_time;         // Elapsed time, s
        std::shared_ptr<Cell_t> p_cell;         // Cell particle currently in
        double                  p_energy;       // Energy, eV
        double                  p_speed;        // Speed, m/s

    public:
        // Constructor: Pass position and direction and initialize particle
        Particle_t( const Point_t& p1, const Point_t& p2, const double E = 2e6, const double t = 0.0, const double w = 1.0 ) :
            p_pos(p1), p_dir(p2), p_time(t), p_weight(w) { setEnergy(E); }
        ~Particle_t() {};

        // Getters
        Point_t                 pos()    const; // Position
        Point_t                 dir()    const; // Direction 
        bool                    alive()  const; // Alive/dead flag
        double                  weight() const; // Weight
        double                  time()   const; // Time
        double                  energy() const; // Energy
        double                  speed()  const; // Speed
        std::shared_ptr<Cell_t> cell()   const; // Particle cell

        // Setters
        void setDirection( const Point_t& p );                 // Direction
        void setWeight   ( const double w );                   // Weight
        void setCell     ( const std::shared_ptr<Cell_t>& C ); // Particle cell
        void setTime     ( const double t );                   // Elapsed time
        void setEnergy   ( const double E );                   // Energy, and speed
        void setSpeed    ( const double v );                   // Speed, and energy

        // Kill particle: live/die flag --> False, weight = 0
        void kill();		                                                   
		
        // Move particle a distance dmove along its current trajectory
        void move( const double dmove );                                           

        // Scatter with scattering angle mu0, with nucleus having relative atomic mass of A
        void scatter( const double mu0, const double A );                                          
		
        // Search and set particle cell (Where am I now?)
        void searchCell( const std::vector<std::shared_ptr<Cell_t>>& Cell ); 
};


#endif
