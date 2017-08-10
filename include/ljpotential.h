// 2017 Artur Avkhadiev
/*! \file ljpotential.h
*/
#ifndef POLYMER_LJPOTENTIAL_H
#define POLYMER_LJPOTENTIAL_H
#include <string>
#include "vector.h"
class LJPotential {
    private:
        double _epsilon;
        double _sigma;
    public:
        double get_epsilon() const;
        double get_sigma() const;
        // get a string representation of the potential
        std::string get_str() const;
        // writes out the simulation's parameters into file
        // outdir/sim_name_potential.cfg
        void writeout_parameters_to_file(std::string outdir,
            std::string sim_name);
        // calculates the pair potential energy of interaction between
        // atom at r1 and atom at r2.
        double calculate_pair_potential(double inv_rijsq);
        // calculates the negative of the pair virial function between r1 and r2
        double calculate_neg_pair_virial(double inv_rijsq);
        // calculates fij,u
        //  where f = fij * rij is the force on atom i from atom j
        double calculate_fstrength_over_r(double inv_rijsq);
        LJPotential();
        LJPotential(double epsilon,
            double sigma,
            bool calculate_prefactors = true);
        ~LJPotential();
};
::std::ostream& operator<<(::std::ostream& os, const LJPotential& potential);
#endif
