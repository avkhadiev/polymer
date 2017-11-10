// 2017 Artur Avkhadiev
/*! \file ljpotential.h
*/
#ifndef POLYMER_LJPOTENTIAL_H
#define POLYMER_LJPOTENTIAL_H
#include <string>
#include "vector.h"
#include "potential.h"
class LJPotential :
    public Potential {
    protected:
        double _epsilon;
        double _sigma;
        double _sigmasq;
    public:
        double get_epsilon() const;
        double get_sigma() const;
        // get a string representation of the potential
        virtual std::string get_str() const;
        // writes out the simulation's parameters into file
        // outdir/sim_name_potential.cfg
        virtual void writeout_parameters_to_file(std::string outdir,
            std::string sim_name);
        // calculates the pair potential energy of interaction between
        // atom at r1 and atom at r2.
        virtual double calculate_pair_potential(double inv_rijsq);
        // calculates the negative of the pair virial function between r1 and r2
        virtual double calculate_neg_pair_virial(double inv_rijsq);
        // calculates fij,u
        //  where f = fij * rij is the force on atom i from atom j
        virtual double calculate_fstrength_over_r(double inv_rijsq);
        LJPotential();
        LJPotential(double epsilon,
            double sigma);
        ~LJPotential();
};
class AdjustedLJPotential :
    public LJPotential {
    protected:
        double _rc;                         /**> cutoff in sigma      */
        double _rcsq;                       /**> square of the cutoff */
        double _corr1;
        double _corr2;
        /** calculate _corr1 and _corr2 according to Stoddard and Ford, 1973 */
        double _calculate_corr1();
        double _calculate_corr2();
    public:
        double get_rc() const {return _rc;};
        // get a string representation of the potential
        virtual std::string get_str() const;
        // writes out the simulation's parameters into file
        // outdir/sim_name_potential.cfg
        virtual void writeout_parameters_to_file(std::string outdir,
            std::string sim_name);
        // calculates the pair potential energy of interaction between
        // atom at r1 and atom at r2.
        virtual double calculate_pair_potential(double inv_rijsq);
        // calculates the negative of the pair virial function between r1 and r2
        virtual double calculate_neg_pair_virial(double inv_rijsq);
        // calculates fij,u
        //  where f = fij * rij is the force on atom i from atom j
        virtual double calculate_fstrength_over_r(double inv_rijsq);
        AdjustedLJPotential(double cutoff);
        AdjustedLJPotential(double epsilon, double sigma, double cutoff);
        ~AdjustedLJPotential();
};
::std::ostream& operator<<(::std::ostream& os, const LJPotential& potential);
#endif
