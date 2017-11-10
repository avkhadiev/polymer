// 2017 Artur Avkhadiev
/*! \file potential.h
*/
#ifndef POLYMER_POTENTIAL_H
#define POLYMER_POTENTIAL_H
#include <string>
#include "vector.h"
class Potential {
    public:
        // get a string representation of the potential
        virtual std::string get_str() const {return "Base Potential Class (no interaction)\n";}
        // writes out the simulation's parameters into file
        // outdir/sim_name_potential.cfg
        virtual void writeout_parameters_to_file(std::string outdir,
            std::string sim_name) {};
        // calculates the pair potential energy of interaction between
        // atom at r1 and atom at r2.
        virtual double calculate_pair_potential(double inv_rijsq){return 0.0;};
        // calculates the negative of the pair virial function between r1 and r2
        virtual double calculate_neg_pair_virial(double inv_rijsq){return 0.0;};
        // calculates fij,u
        //  where f = fij * rij is the force on atom i from atom j
        virtual double calculate_fstrength_over_r(double inv_rijsq) {return 0.0;};
        Potential();
        ~Potential();
};
::std::ostream& operator<<(::std::ostream& os, const Potential& potential);
#endif
