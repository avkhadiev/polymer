// 2017 Artur Avkhadiev
/*! \file general_observables.h
*/
#ifndef GENERAL_POLYMER_OBSERVABLES_H
#define GENERAL_POLYMER_OBSERVABLES_H
#include "observable.h"
namespace observable {
    extern Vector x_vec;
    extern Vector y_vec;
    extern Vector z_vec;
    /**
    * Time, updated in main loop
    */
    class Time:
        public Observable {
    public:
        virtual bool update_method_specified() const {return true;};
        virtual void update(const simple::AtomState& state);
        ~Time(){};
        Time();
    };
    /**
    * kinetic energy, updated in integrator
    */
    class KE :
        public Observable {
    public:
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        void update(const simple::Atom& atom, double mass);
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
    * kinetic temperature, updated in main loop
    */
    class KineticTemperature :
    public Observable {
    protected:
        bool _remove_linear_momentum;
        bool _remove_angular_momentum;
        size_t _ndof;
    public:
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        void update(const simple::Atom& atom, double mass, size_t ndof);
        size_t ndof() const {return _ndof;};
        bool remove_linear_momentum() const {return  _remove_linear_momentum;};
        bool remove_angular_momentum() const{return _remove_angular_momentum;};
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
    * virial
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
    * a component of linear momentum, updated in main loop
    */
    class LinMomComponent :
        public Observable {
    private:
        Vector _component;
        void _update(Vector momentum, Vector component);
    public:
        Vector component() const {return _component;};
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        void update(const simple::Atom& atom, double mass);
        LinMomComponent(Vector component,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false // use precision notation when writing)
        );
        ~LinMomComponent(){};
    };
    /**
    * a component of angular momentum, updated in main loop
    */
    class AngMomComponent :
        public Observable {
    private:
        Vector _component;
        void _update(Vector momentum, Vector component);
    public:
        Vector component() const {return _component;};
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        void update(const simple::Atom& atom, double mass);
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
        public Observable {
    protected:
        Vector _component;
        size_t _molecule_index;     // which polymer does this RCM refer to?
        virtual Vector rcm(const simple::AtomState& state){
            return vector(0.0, 0.0, 0.0);
        };
        void _update(Vector rcm, Vector component);
    public:
        Vector component() const {return _component;};
        size_t molecule_index() const {return _molecule_index;};
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        RCMComponent(Vector component,
            size_t molecule_index,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing)
        );
        ~RCMComponent(){};
    };
} // namespace observable
#endif
