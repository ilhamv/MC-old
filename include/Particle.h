#ifndef PARTICLE_H
#define PARTICLE_H

#include <memory> // shared_ptr
#include <vector> // vector

#include "Point.h"

class Cell;
class Surface;

class Particle
{
    private:
        Point p_pos; // cm
        Point p_dir;
        bool p_alive;
        double p_energy; // eV
        double p_energy_old;
        double p_speed; // cm/s
        double p_weight;
        double p_time; //s
        double p_time_old;
        int    p_tdmc;
        std::shared_ptr<Cell> p_cell;
        std::shared_ptr<Cell> p_cell_old;
        std::shared_ptr<Surface> p_surface_old;

    public:
        Particle( const Point& p1, const Point& p2, const double E, 
                  const double t, const double w, const int td, 
                  const std::shared_ptr<Cell> c ) :
            p_alive(true), p_pos(p1), p_dir(p2), p_time(t), p_weight(w), 
            p_tdmc(td), p_cell(c) { set_energy(E); }
        ~Particle() {};

        // Getters
        Point pos() const;
        Point dir() const;
        bool alive() const;
        double weight() const;
        double time() const;
        double time_old() const;
        double energy() const;
        double energy_old() const;
        double speed() const;
        std::shared_ptr<Cell> cell() const;
        std::shared_ptr<Cell> cell_old() const;
        std::shared_ptr<Surface> surface_old() const;
        int tdmc() const;

        // Setters
        void set_direction( const Point& p );
        void set_weight( const double w );
        void set_cell( const std::shared_ptr<Cell> C );
        void set_energy( const double E );
        void set_speed( const double v );
        void set_surface_old( const std::shared_ptr<Surface> S );

        // Modifiers
        void kill(); 
        void move( const double dmove );                                      
        void increment_tdmc();
};


#endif // PARTICLE_H
