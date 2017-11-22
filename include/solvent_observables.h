// 2017 Artur Avkhadiev
/*! \file dynamic_observables.h
*/
#ifndef POLYMER_SOLVENT_OBSERVABLES_H
#define POLYMER_SOLVENT_OBSERVABLES_H
#include "observable.h"
namespace solvent {
    /**
    * calculates potential energy of the solvent
    * ForceLoop knows how to work with this observable, no update functions
    */
    class V :
        public TimeLogObservable {
    public:
        V();
        ~V(){};
    };
    /**
    * calculates average potential energy of the solvent
    * ForceLoop knows how to work with this observable, no update functions
    */
    class AvgV :
        public TimeLogObservable,
        public AvgObservable {
    public:
        AvgV(V& ve, double acc = 0.0);
        ~AvgV(){};
    };
    /**
    * calculates kinetic energy of the solvent
    */
    class KE :
        public TimeLogObservable {
    public:
        virtual void update(const simple::Solvent& molecule);
        virtual void update(const simple::AtomState& state);
        double K(const simple::AtomState& state);
        KE();
        ~KE(){};
    };
    /**
    * calculates average kinetic energy of the solvent
    */
    class AvgKE :
        public TimeLogObservable,
        public AvgObservable {
    public:
        AvgKE(KE& ke, double acc = 0.0);
        ~AvgKE(){};
    };
    /**
    * calculates magnitude of linear momentum of the solvent
    */
    class Momentum :
        public TimeLogObservable {
    private:
        Vector _acc;            /**> intermediate value stored here */
    public:
        virtual void zero() {Observable::zero(); _acc = vector(0.0, 0.0, 0.0);};
        virtual void update(const simple::Solvent& molecule);
        virtual void update(const simple::AtomState& state);
        double p(const simple::AtomState& state);
        Momentum();
        ~Momentum(){};
    };
} // namespace solvent
#endif
