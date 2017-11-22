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
        virtual double rijsq(Vector rij);
    public:
        double get_epsilon() const;
        double get_sigma() const;
        virtual std::string get_str() const;
        virtual void writeout_parameters_to_file(std::string outdir,
            std::string sim_name);
        virtual Vector rij(Vector ri, Vector rj){return subtract(ri, rj);};
        virtual double vij(Vector ri, Vector rj);
        virtual double wij(Vector ri, Vector rj);
        // calculates force on atom i from atom j
        virtual Vector fij(Vector ri, Vector rj);
        LJPotential();
        LJPotential(double epsilon, double sigma);
        ~LJPotential();
};
class AdjustedLJPotential :
    public LJPotential {
    protected:
        double _rc;                         /**> cutoff in sigma      */
        double _rcsq;                       /**> square of the cutoff */
        double _box;                        /**> box size             */
        double _corr1;
        double _corr2;
        /** calculate _corr1 and _corr2 according to Stoddard and Ford, 1973 */
        double _calculate_corr1();
        double _calculate_corr2();
        virtual double rijsq(Vector rij);
    public:
        double get_rc() const {return _rc;};
        virtual std::string get_str() const;
        virtual void writeout_parameters_to_file(std::string outdir,
            std::string sim_name);
        virtual Vector rij(Vector ri, Vector rj);
        virtual double vij(Vector ri, Vector rj);
        virtual double wij(Vector ri, Vector rj);
        // calculates force on atom i from atom j
        virtual Vector fij(Vector ri, Vector rj);
        AdjustedLJPotential(double cutoff, double box);
        AdjustedLJPotential(double epsilon, double sigma, double cutoff, double box);
        ~AdjustedLJPotential();
};
::std::ostream& operator<<(::std::ostream& os, const LJPotential& potential);
#endif
