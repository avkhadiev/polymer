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
        // calculates the separation vector between
        // atom at ri and atom at rj
        virtual Vector rij(Vector ri, Vector rj){
            return subtract(ri, rj);
        };
        // calculates the pair potential energy of interaction between
        // atom at ri and atom at rj
        virtual double vij(Vector ri, Vector rj){return 0.0;};
        // calculates the pair virial function between ri and rj
        virtual double wij(Vector ri, Vector rj){return 0.0;};
        // calculates force on atom i from atom j
        virtual Vector fij(Vector ri, Vector rj){
            return vector(0.0, 0.0, 0.0);
        };
        Potential();
        ~Potential();
};
::std::ostream& operator<<(::std::ostream& os, const Potential& potential);
#endif
