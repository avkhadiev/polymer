// 2017 Artur Avkhadiev
/*! \file dynamic_observables.h
*/
#ifndef POLYMER_DYNAMIC_OBSERVABLES_H
#define POLYMER_DYNAMIC_OBSERVABLES_H
#include "observable.h"
namespace dynamic{
    /**
    * calculates kinetic energy of the polymer
    * value = internal kinetic energy (motion wrt R_CM)
    */
    class IntKE :
        public TimeLogObservable {
    private:
        double _m;              /**> atomic mass */
    public:
        virtual void update(const simple::Atom &atom);
        virtual void update(const simple::AtomPolymer& polymer);
        virtual void update(const simple::AtomState& state);
        double K(const simple::AtomState& state);
        IntKE(double m);
        ~IntKE(){};
    };
    /**
    * calculates angular momentum of the polymer
    * value = norm of angular momentum
    */
    class LNorm :
        public TimeLogObservable {
    private:
        double _m;              /**> atomic mass */
        Vector _acc;            /**> intermediate value stored here */
    public:
        virtual void zero() {Observable::zero(); _acc = vector(0.0, 0.0, 0.0);};
        virtual void update(const simple::Atom &atom);
        virtual void update(const simple::AtomPolymer& polymer);
        virtual void update(const simple::AtomState& state);
        // double L(const simple::AtomState& state);
        LNorm(double m);
        ~LNorm(){};
    };
    /**
    * calculates angular momentum's projection on the specified axis
    */
    class LProj :
        public TimeLogObservable {
    private:
        double _m;              /**> atomic mass                    */
        Vector _acc;            /**> intermediate value stored here */
        Vector _axis;           /**> projected onto this vector     */
    public:
        virtual std::string to_string() const;
        Vector L(const simple::AtomState& state);
        virtual void update(std::string obs_string);
        virtual void update(const simple::Atom &atom);
        virtual void update(const simple::AtomPolymer& polymer);
        virtual void update(const simple::AtomState& state);
        // double L(const simple::AtomState& state);
        LProj(double m, Vector axis = vector(0.0, 0.0, 1.0));
        ~LProj(){};
    };
} // namespace dynamic
#endif
