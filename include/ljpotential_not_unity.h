// 2017 Artur Avkhadiev
/*! \file ljpotential.h
* does NOT measure energy, distance in units of epsilon, sigma
* (they are included in the calculations)
*/
#ifndef POLYMER_LJPOTENTIAL_H
#define POLYMER_LJPOTENTIAL_H
#include <string>
#include "vector.h"
class LJPotential {
    private:
        double _epsilon;
        double _epsilon4;
        double _epsilon24;
        double _sigma;
        double _sigma6;
        double _sigma12;
        bool _calculate_prefactors;
    public:
        double get_epsilon() const;
        double get_sigma() const;
        bool get_calculate_prefactors();
        void set_epsilon(double epsilon);
        void set_sigma(double sigma);
        void set_calculate_prefactors(bool calculate_prefactors);
        // get a string representation of the potential
        std::string get_str() const;
        // writes out the simulation's parameters into file
        // outdir/sim_name_potential.cfg
        void writeout_parameters_to_file(std::string outdir,
            std::string sim_name);
        // reads in simulation's parameterse from file
        //  indir/sim_name_potential.cfg and sets epsilon and sigma accordingly
        void read_parameters_from_file(std::string indir,
            std::string sim_name);
        // calculates the pair potential energy of interaction between
        // atom at r1 and atom at r2.
        // if calculate_prefactors = true, calculates with the proper prefactor
        // 4 * epsilon; calculates without the prefactor otherwise.
        double calculate_pair_potential(double inv_rijsq);
        // calculates the negative of the pair virial function between r1 and r2
        // if calculate_prefactors = true, calculates with the proper prefactor
        // 24 * epsilon; calculates without the prefactor otherwise.
        double calculate_neg_pair_virial(double inv_rijsq);
        // calculates fij,
        //  where f = fij * rij is the force on atom i from atom j
        // if calculate_prefactors = true, calculates with the proper prefactor
        // 24 * epsilon; calculates without the prefactor otherwise.
        double calculate_fstrength_over_r(double inv_rijsq);
        LJPotential();
        LJPotential(double epsilon,
            double sigma,
            bool calculate_prefactors = true);
        ~LJPotential();
};
::std::ostream& operator<<(::std::ostream& os, const LJPotential& potential);
#endif
