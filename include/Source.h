#ifndef SOURCE_H
#define SOURCE_H

#include <vector>     
#include <memory>     

#include "Particle.h"
#include "Point.h"
#include "Distribution.h"

class Cell;


//=============================================================================
// Particle Source
//=============================================================================

class Source
{
    public:
 	 Source() {};
	~Source() {};

	virtual Particle get_source() = 0;
};
class SourceDelta : public Source
{
    private:
        const Particle s_P;

    public:
        SourceDelta( const Particle& P ): s_P(P) {};
        ~SourceDelta() {};

        Particle get_source();
};
class SourcePoint : public Source
{
    private:
	const Point s_pos;
        const std::shared_ptr<Cell> s_cell;
    	const std::shared_ptr<Distribution<Point>>  s_dir;
    	const std::shared_ptr<Distribution<double>> s_energy;

    public:
	SourcePoint( const Point p, const std::shared_ptr<Cell> cell,
                     const std::shared_ptr< Distribution<Point>> dir,
		     const std::shared_ptr< Distribution<double>> enrg ) :
            s_pos(p), s_cell(cell), s_dir(dir), s_energy(enrg) {};
	~SourcePoint() {};

	Particle get_source();
};


//=============================================================================
// Source Bank
//=============================================================================

class SourceBank
{
    public:
        std::vector<std::pair<std::shared_ptr<Source>,double>> sources;
        double total = 0.0;
        std::vector<double> p; // Source probability
	
        SourceBank() {};
        ~SourceBank() {};
		
        void add_source( const std::shared_ptr<Source> S, const double prob );
        void reset();
        Particle get_source();
        void set_up();
};


#endif // SOURCE_H
