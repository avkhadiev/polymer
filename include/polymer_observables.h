// 2017 Artur Avkhadiev
/*! \file polymer_observables.h
*/
#ifndef POLYMER_POLYMER_OBSERVABLES_H
#define POLYMER_POLYMER_OBSERVABLES_H
#include "observable.h"
#include "general_observables.h"
#include "geodesic_observables.h"
namespace polymer{
    /**
    * kinetic energy of a polymer, updated in integrator
    */
    class KE :
        public observable::KE {
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        void update(const simple::AtomPolymer& molecule);
        ~KE(){};
        KE( bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
    };
    /**
    * potential energy of a polymer, updated in force loop
    */
    class PE :
        public observable::PE {
    public:
        virtual bool update_method_specified() const {return false;};
        ~PE(){};
        PE( bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
    };
    /**
    * polymer virial due to potential forces, updated in force loop
    */
    class Virial :
        public observable::Virial {
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
    * polymer virial due to constraint forces, updated in RATTLE loop
    */
    class ConstraintVirial :
        public observable::Virial {
    public:
        virtual bool update_method_specified() const {return false;};
        ~ConstraintVirial(){};
        ConstraintVirial( bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
    };
    /**
    * kinetic temperature of the polymer
    */
    class KineticTemperature :
        public observable::KineticTemperature {
    protected:
        size_t _ndof;
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        size_t ndof() const {return _ndof;};
        void update(const simple::AtomPolymer& polymer);
        ~KineticTemperature(){};
        KineticTemperature(
            bool remove_linear_momentum,
            bool remove_angular_momentum,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
    };
    /**
    * a component of linear momentum of a polymer, updated in main loop
    */
    class LinMomComponent :
        public observable::LinMomComponent {
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        void update(const simple::AtomPolymer& polymer);
        LinMomComponent(Vector component,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,  // print out inst value in log & state?
            bool e_format = false // use precision notation when writing)
        );
        ~LinMomComponent(){};
    };
    /**
    * a component of angular momentum of a polymer, updated in main loop
    */
    class AngMomComponent :
        public observable::AngMomComponent {
    private:
        void _update(Vector momentum, Vector component);
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        void update(const simple::AtomPolymer& polymer);
        AngMomComponent(Vector component,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing)
        );
        ~AngMomComponent(){};
    };
    /**
    * a component of CM of a polymer, updated in main loop
    */
    class RCMComponent :
        public observable::RCMComponent {
    private:
        virtual Vector rcm(const simple::AtomState& state);
        void _update(Vector rcm, Vector component);
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        RCMComponent(Vector component,
            size_t molecule_index,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing)
        );
        ~RCMComponent(){};
    };
    /**
    * the length of a given bond of the polymer
    */
    class BondLength :
        public geodesic::LinkObservable {
    private:
        size_t _molecule_index;
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        BondLength(size_t molecule_index,
            size_t bond_index,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = true  // use precision notation when writing)
        );
        ~BondLength(){};
    };
    /**
    * a component of an atomic position of the polymer
    */
    class PolAtomPosComponent :
        public Observable {
    private:
        Vector _component;
        size_t _molecule_index;
        size_t _atom_index;
        void _update(Vector r, Vector component);
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        void amend_names();            // indicate which component in the name
        PolAtomPosComponent(Vector component,
            size_t molecule_index,
            size_t atom_index,
            bool calculate_mean = false,
            bool calculate_error = false,
            bool print_inst_val = false,// print out inst value in log & state?
            bool e_format = false       // use precision notation when writing?
        );
        ~PolAtomPosComponent(){};
    };
} // namespace polymer
#endif
