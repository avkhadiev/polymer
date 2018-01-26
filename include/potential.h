// 2017 Artur Avkhadiev
/*! \file potential.h
*/
#ifndef POLYMER_POTENTIAL_H
#define POLYMER_POTENTIAL_H
#include <string>
#include "vector.h"
class Potential {
    public:
        virtual double get_epsilon() const = 0;
        virtual double get_sigma() const = 0;
        // get a string representation of the potential
        virtual std::string get_str() const = 0;
        // writes out the simulation's parameters into file
        // outdir/sim_name_potential.cfg
        virtual void writeout_parameters_to_file(std::string outdir,
            std::string sim_name) = 0;
        virtual void zero_observables() = 0;
        // calculates the separation vector between
        // atom at ri and atom at rj
        virtual Vector rij(Vector ri, Vector rj) const = 0;
        // calculates pair potential between atom i and atom j
        virtual double vij(Vector ri, Vector rij) const = 0;
        // calculates force on atom i from atom j
        virtual Vector fij(Vector ri, Vector rj,
            bool calculate_observables = false) = 0;
        Potential();
        ~Potential();
};
::std::ostream& operator<<(::std::ostream& os, const Potential& potential);
#endif
