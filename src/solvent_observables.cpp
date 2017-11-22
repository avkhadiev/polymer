// 2017 Artur Avkhadiev
/*! \file solvent_observables.cpp
*/
#include "../include/solvent_observables.h"
#include "../include/simple_solvent.h"
namespace solvent{
    /**************************************************************************
    * LJ Potential Energy
    **************************************************************************/
    V::V() :
        Observable("Solvent Potential", "svlj", "\\varepsilon", ""),
        TimeLogObservable(){}
    /**************************************************************************
    * Average LJ Potential Energy
    **************************************************************************/
    AvgV::AvgV(V& ve, double acc ) :
        Observable("Average LJ Potential", "savgvlj", "\\varepsilon", ""),
        TimeLogObservable(),
        AvgObservable(ve, acc){}
    /**************************************************************************
    * Kinetic Energy
    **************************************************************************/
    KE::KE() :
        Observable("Translational Kinetic Energy", "sk", "\\varepsilon", ""),
        TimeLogObservable(){}
    void KE::update(const simple::Solvent &molecule){
        _value += simple::Solvent::m() * normsq(molecule.v()) / 2;
    }
    void KE::update(const simple::AtomState &state){
        for(const simple::Solvent& solvent : state.solvents){
            update(solvent);
        }
    }
    double KE::K(const simple::AtomState &state){
        double K = 0.0;
        for(const simple::Solvent& solvent : state.solvents){
            K += simple::Solvent::m() * normsq(solvent.v()) / 2;
        }
        return K;
    }
    /**************************************************************************
    * Average Kinetic Energy
    **************************************************************************/
    AvgKE::AvgKE(KE& ke, double acc) :
        Observable("Average Kinetic Energy", "savgik",
            "\\varepsilon",
            ""),
        TimeLogObservable(),
        AvgObservable(ke, acc){}
    /**************************************************************************
    * Kinetic Energy
    **************************************************************************/
    Momentum::Momentum() :
        Observable("Linear Momentum", "sp", "\\varepsilon", ""),
        TimeLogObservable(),
        _acc(vector(0.0, 0.0, 0.0)) {}
    void Momentum::update(const simple::Solvent &molecule){
        _acc += multiply(molecule.v(), simple::Solvent::m());
    }
    void Momentum::update(const simple::AtomState &state){
        _acc = vector(0.0, 0.0, 0.0);
        for(const simple::Solvent& solvent : state.solvents){
            update(solvent);
        }
        Observable::update(norm(_acc));
    }
    double Momentum::p(const simple::AtomState &state){
        Vector p = vector(0.0, 0.0, 0.0);
        for(const simple::Solvent& solvent : state.solvents){
            p += multiply(solvent.v(), simple::Solvent::m());
        }
        return norm(p);
    }
} // namespace solvent
