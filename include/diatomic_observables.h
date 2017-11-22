// 2017 Artur Avkhadiev
/*! \file diatomic_observables.h
*/
#ifndef POLYMER_DIATOMIC_OBSERVABLES_H
#define POLYMER_DIATOMIC_OBSERVABLES_H
#include "observable.h"
#include "polymer_observables.h"
namespace diatomic{
    /**
    * assumes the polymer is a diatomic and calculates the length of its
    * bond in the atomic representation
    */
    class BondPosLength :
        public TimeLogObservable {
    public:
        virtual void update(const simple::AtomPolymer& polymer);
        virtual void update(const simple::AtomState& state);
        BondPosLength();
        ~BondPosLength(){};
    };
    /**
    * assumes the polymer is a diatomic and calculates the length of its
    * bond in the atomic representation
    */
    class BondVelLength :
        public TimeLogObservable {
    public:
        virtual void update(const simple::AtomPolymer& polymer);
        virtual void update(const simple::AtomState& state);
        BondVelLength();
        ~BondVelLength(){};
    };
    /**
    * assumes the polymer is a diatomic and calculates the length of its bond's
    * projection on the specified axis
    */
    class BondPosProj:
        public TimeLogObservable {
    protected:
        Vector _axis;
    public:
        virtual std::string to_string() const;
        virtual void update(std::string obs_string);
        virtual void update(const simple::AtomPolymer& polymer);
        virtual void update(const simple::AtomState& state);
        BondPosProj(Vector axis = vector(1.0, 0.0, 0.0));
        ~BondPosProj(){};
    };
    /**
    * assumes the polymer is a diatomic and calculates the length of its bond's
    * velocity projection on the specified axis
    */
    class BondVelProj:
        public TimeLogObservable {
    protected:
        Vector _axis;
    public:
        virtual std::string to_string() const;
        virtual void update(std::string obs_string);
        virtual void update(const simple::AtomPolymer& polymer);
        virtual void update(const simple::AtomState& state);
        BondVelProj(Vector axis = vector(0.0, -1.0, 0.0));
        ~BondVelProj(){};
    };
} // namespace diatomic
#endif
