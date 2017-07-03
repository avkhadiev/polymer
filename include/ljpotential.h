// 2017 Artur Avkhadiev
/*! \file ljpotential.h
*/
#ifndef POLYMER_LJPOTENTIAL_H
#define POLYMER_LJPOTENTIAL_H
#include "vector.h"
#include "atom.h"
class LJPotential {
    private:
        double _epsilon;
        double _sigma;
        double _sigma6;
        double _sigma12;
    public:
        double get_epsilon();
        double get_sigma();
        void set_epsilon(double epsilon);
        void set_sigma(double sigma);
        Vector calculate_force(Atom &on_atom, Atom &from_atom);
        // constructors and a destructor
        LJPotential();
        LJPotential(double epsilon, double sigma);
        ~LJPotential();
};
#endif
