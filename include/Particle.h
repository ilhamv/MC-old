#ifndef _PARTICLE_HEADER_
#define _PARTICLE_HEADER_

#include <memory> // shared_ptr
#include <vector> // vector

#include "Point.h"
#include "Geometry.h"

class Cell_t;
class Surface_t;

class Particle_t
{
    private:
        Point_t p_pos;         
        Point_t p_dir;         
        bool p_alive = true;
        double p_weight;      
        double p_time;        
        double p_time_old;
        std::shared_ptr<Cell_t> p_cell;        
        std::shared_ptr<Cell_t> p_cell_old;
        std::shared_ptr<Surface_t> p_surface_old;
        double p_energy;      
        double p_speed;       

    public:
        Particle_t( const Point_t& p1, const Point_t& p2, const double E, 
                    const double t, const double w ) :
            p_pos(p1), p_dir(p2), p_time(t), p_weight(w) { setEnergy(E); }
        ~Particle_t() {};

        // Getters
        Point_t pos() const;
        Point_t dir() const;
        bool alive() const;
        double weight() const;
        double time() const;
        double time_old() const;
        double energy() const;
        double speed() const;
        std::shared_ptr<Cell_t> cell() const;
        std::shared_ptr<Cell_t> cell_old() const;
        std::shared_ptr<Surface_t> surface_old() const;

        // Setters
        void setDirection( const Point_t& p );
        void setWeight( const double w );
        void setCell( const std::shared_ptr<Cell_t>& C );
        void setTime( const double t );
        void setEnergy( const double E );
        void setSpeed( const double v );
        void set_surface_old( const std::shared_ptr<Surface_t>& S );

        // Modifiers
        void kill();		                                                   
        void move( const double dmove );                                           
        void searchCell( const std::vector<std::shared_ptr<Cell_t>>& Cell ); 
};


#endif
