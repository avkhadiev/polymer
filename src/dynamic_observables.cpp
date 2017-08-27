// 2017 Artur Avkhadiev
/*! \file dynamic_observables.cpp
*/
#include "../include/dynamic_observables.h"
namespace dynamic{
    /**************************************************************************
    * LJ Potential Energy
    **************************************************************************/
    V::V() :
        Observable("LJ Potential", "vlj", "\\varepsilon", ""),
        TimeLogObservable(){}
    /**************************************************************************
    * Average LJ Potential Energy
    **************************************************************************/
    AvgV::AvgV(V& ve, double acc ) :
        Observable("Average LJ Potential", "avgvlj", "\\varepsilon", ""),
        TimeLogObservable(),
        AvgObservable(ve, acc){}
    /**************************************************************************
    * Negative Virial
    **************************************************************************/
    NegW::NegW() :
        Observable("Negative Virial", "w", "\\varepsilon", ""),
        TimeLogObservable(){}
    /**************************************************************************
    * Average Negative Virial
    **************************************************************************/
    AvgNegW::AvgNegW(NegW& w, double acc ) :
        Observable("Average Negative Virial", "avgw", "\\varepsilon", ""),
        TimeLogObservable(),
        AvgObservable(w, acc){}
    /**************************************************************************
    * Negative Constraint Virial
    **************************************************************************/
    NegWC::NegWC() :
        Observable("Negative Constraint Virial", "wc", "\\varepsilon", ""),
        TimeLogObservable(){}
    /**************************************************************************
    * Average Negative Virial
    **************************************************************************/
    AvgNegWC::AvgNegWC(NegWC& wc, double acc) :
        Observable("Average Negative Constraint Virial",
            "avgwc", "\\varepsilon", ""),
        TimeLogObservable(),
        AvgObservable(wc, acc){}
    /**************************************************************************
    * Internal Kinetic Energy
    **************************************************************************/
    IntKE::IntKE(double m) :
        Observable("Internal Kinetic Energy", "ik", "\\varepsilon", ""),
        TimeLogObservable(),
        _m(m){}
    void IntKE::update(const simple::Atom &atom){
        _value += _m * normsq(atom.velocity) / 2;
    }
    void IntKE::update(const simple::AtomPolymer &polymer){
        for(const simple::Atom& atom : polymer.atoms){
            update(atom);
        }
    }
    void IntKE::update(const simple::AtomState &state){
        // assumes there is only one molecule in the state
        update(state.polymers.at(0));
    }
    double IntKE::K(const simple::AtomState &state){
        double K = 0.0;
        for(const simple::Atom& atom : state.polymers.at(0).atoms){
            K += _m * normsq(atom.velocity) / 2;
        }
        return K;
    }
    /**************************************************************************
    * Average Internal Kinetic Energy
    **************************************************************************/
    AvgIntKE::AvgIntKE(IntKE& ke, double acc) :
        Observable("Average Internal Kinetic Energy", "avgik",
            "\\varepsilon",
            ""),
        TimeLogObservable(),
        AvgObservable(ke, acc){}
    /**************************************************************************
    * Angular Momentum Norm
    **************************************************************************/
    LNorm::LNorm(double m) :
        Observable("Angular Momentum Norm", "ln",
            "m * \\sigma^2 / \\tau",
            ""),
        TimeLogObservable(),
        _m(m),
        _acc(vector(0.0, 0.0, 0.0)) {}
    void LNorm::update(const simple::Atom &atom){
        _acc += cross(atom.position, multiply(atom.velocity, _m));
    }
    void LNorm::update(const simple::AtomPolymer &polymer){
        _acc = vector(0.0, 0.0, 0.0);
        for (const simple::Atom& atom : polymer.atoms){
            update(atom);
        }
        Observable::update(norm(_acc));
    }
    void LNorm::update(const simple::AtomState &state){
        _acc = vector(0.0, 0.0, 0.0);
        // assumes there is only one molecule in the state
        update(state.polymers.at(0));
    }
    /**************************************************************************
    * Angular Momentum Projection
    **************************************************************************/
    LProj::LProj(double m, Vector axis) :
        Observable("Angular Momentum Projection on " + vector_to_string(axis),
            "lp", "", ""),
        TimeLogObservable(),
        _m(m),
        _acc(vector(0.0, 0.0, 0.0)),
        _axis(divide(axis, norm(axis))){}
    std::string LProj::to_string() const {
        return "lproj " + vector_to_string(_axis) + " " + print_value();
    }
    void LProj::update(std::string obs_string){
        std::istringstream ss(obs_string.c_str());
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> words(begin, end);
        _axis.x = atof(words.at(1).c_str());
        _axis.y = atof(words.at(2).c_str());
        _axis.z = atof(words.at(3).c_str());
        Observable::update(atof(words.back().c_str()));
    }
    void LProj::update(const simple::Atom &atom){
        _acc += cross(atom.position, multiply(atom.velocity, _m));
    }
    void LProj::update(const simple::AtomPolymer &polymer){
        _acc = vector(0.0, 0.0, 0.0);
        for (const simple::Atom& atom : polymer.atoms){
            update(atom);
        }
        Observable::update(dot(_acc, _axis)/norm(_acc));
    }
    void LProj::update(const simple::AtomState &state){
        _acc = vector(0.0, 0.0, 0.0);
        // assumes there is only one molecule in the state
        update(state.polymers.at(0));
    }
    Vector LProj::L(const simple::AtomState &state){
        Vector L = vector(0.0, 0.0, 0.0);
        for (const simple::Atom& atom : state.polymers.at(0).atoms){
            L += cross(atom.position, multiply(atom.velocity, _m));
        }
        return L;
    }
} // namespace dynamic
