// 2017 Artur Avkhadiev
/*! \file ljpotential.h
*/
#ifndef POLYMER_LJPOTENTIAL_H
#define POLYMER_LJPOTENTIAL_H
#include <string>
#include "vector.h"
#include "potential.h"
#include "observable.h"
class LJPotential :
    public Potential {
    protected:
        double _epsilon;
        double _sigma;
        double _sigmasq;
        virtual double rijsq(Vector rij) const;
        typedef struct obs_t {
            bool is_set;
            Observable *ptr;
        } ObservableStruct;
        void _setup_observable(ObservableStruct* obs, Observable* obs_ptr);
    private:
        ObservableStruct _pe, _w;
        // calculates the pair virial
        double _wij(Vector ri, Vector rj) const;
        // add pair potential to the  potential energy observable
        void _update_vij(Vector ri, Vector rj);
    public:
        virtual double get_epsilon() const;
        virtual double get_sigma() const;
        virtual void zero_observables();
        virtual std::string get_str() const;
        virtual void writeout_parameters_to_file(std::string outdir,
            std::string sim_name);
        virtual Vector rij(Vector ri, Vector rj) const{
            return subtract(ri, rj);
        };
        virtual double vij(Vector ri, Vector rij) const;
        // calculates force on atom i from atom j
        virtual Vector fij(Vector ri, Vector rj,
            bool calculate_observables = false);
        LJPotential(Observable *pe = NULL, Observable *virial = NULL);
        LJPotential(double epsilon, double sigma,
                    Observable *pe = NULL, Observable *virial = NULL);
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
        double _calculate_corr1() const;
        double _calculate_corr2() const;
        virtual double rijsq(Vector rij) const;
    private:
        ObservableStruct _pe_shftd, _pe_unshftd, _w_shftd, _w_unshftd;
        // calculates the pair virial
        double _wij(Vector ri, Vector rj, bool calculate_observables);
        // add pair potential to the  potential energy observable
        void _update_vij(Vector ri, Vector rj);
    public:
        double get_rc() const {return _rc;};
        virtual void zero_observables();
        virtual std::string get_str() const;
        virtual void writeout_parameters_to_file(std::string outdir,
            std::string sim_name);
        virtual Vector rij(Vector ri, Vector rj) const;
        virtual double vij(Vector ri, Vector rij) const;
        // calculates force on atom i from atom j
        virtual Vector fij(Vector ri, Vector rj,
            bool calculate_observables = false);
        AdjustedLJPotential();
        AdjustedLJPotential(double cutoff, double box,
            Observable *pe_shftd = NULL, Observable* _pe_unshftd = NULL,
            Observable *w_shftd = NULL, Observable* w_unshftd = NULL);
        AdjustedLJPotential(double epsilon, double sigma,
            double cutoff, double box,
            Observable *pe_shftd = NULL, Observable* _pe_unshftd = NULL,
            Observable *w_shftd = NULL, Observable* w_unshftd = NULL);
        ~AdjustedLJPotential();
};
::std::ostream& operator<<(::std::ostream& os, const LJPotential& potential);
#endif
