// 2017 Artur Avkhadiev
/*! \file dynamic_observables.h
*/
#ifndef POLYMER_SOLVENT_OBSERVABLES_H
#define POLYMER_SOLVENT_OBSERVABLES_H
#include "observable.h"
#include <stdexcept>
namespace solvent {
    /**
    * kinetic energy, updated in integrator
    */
    class KE :
        public Observable {
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        virtual void update(const simple::Solvent& molecule);
        ~KE(){};
        KE( bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
    };
    /**
    * potential energy, updated in force loop
    */
    class PE :
        public Observable {
    private:
        bool _is_shifted;
    public:
        virtual bool update_method_specified() const {return false;};
        ~PE(){};
        PE( bool shifted,          // is this for shifted or unshifted potential
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
    };
    /**
    * virial, updated in force loop
    */
    class Virial :
        public Observable {
    public:
        virtual bool update_method_specified() const {return false;};
        ~Virial(){};
        Virial( bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
    };
    /**
    * a component of linear momentum of the solvent, updated in main loop
    */
    class MomComponent :
        public Observable {
    public:
        typedef enum component_t {X, Y, Z} Component;
    private:
        Component _component;
        void _update(const simple::Solvent& molecule, component_t component);
    public:
        int component() const {return _component;};
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        virtual void update(const simple::Solvent& molecule);
        MomComponent(Component component,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false // use precision notation when writing)
        );
        ~MomComponent(){};
    };
    /**
    * magnitude of linear momentum of the solvent squared, updated in main loop
    */
    class MomMagSq :
        public Observable {
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        virtual void update(const simple::Solvent& molecule);
        MomMagSq(bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
        ~MomMagSq(){};
    };
    /**
    * kinetic temperature of solvents, updated in main loop
    */
    class KinTemp :
        public Observable {
    private:
        const size_t _nc;       // number of constraints (= 3 on com motion)
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        KinTemp(
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
        ~KinTemp(){};
    };
} // namespace solvent
#endif
