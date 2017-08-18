// 2017 Artur Avkhadiev
/*! \file diatomic_observables.cpp
*/
#include <cmath>                                    /**> fabs() */
#include "../include/diatomic_observables.h"
#include "../include/simple_polymer.h"
namespace diatomic{
    BondPosLength::BondPosLength() :
        Observable("Bond Position Length", "bpl", "", ""),
        TimeLogObservable(){}
    void BondPosLength::update(const simple::AtomPolymer &polymer){
        // assumes the molecule is a diatomic
        Observable::update(norm(subtract(
                polymer.atoms.at(0).position,
                polymer.atoms.at(1).position)));
    }
    void BondPosLength::update(const simple::AtomState &state){
        // assumes there is only one molecule in the state
        update(state.polymers.at(0));
    }
    BondVelLength::BondVelLength() :
        Observable("Bond Velocity Length", "bvl", "", ""),
        TimeLogObservable(){}
    void BondVelLength::update(const simple::AtomPolymer &polymer){
        // assumes the molecule is a diatomic
        Observable::update(norm(subtract(
                polymer.atoms.at(0).velocity,
                polymer.atoms.at(1).velocity)));
    }
    void BondVelLength::update(const simple::AtomState &state){
        // assumes there is only one molecule in the state
        update(state.polymers.at(0));
    }
    BondPosProj::BondPosProj(Vector axis) :
        Observable("Bond Position Projection on " + vector_to_string(axis),
            "bpp", "1", ""),
        TimeLogObservable(),
        _axis(divide(axis, norm(axis))){}
    std::string BondPosProj::to_string() const {
        return "bondposproj " + vector_to_string(_axis) + " " + print_value();
    }
    void BondPosProj::update(std::string obs_string){
        std::istringstream ss(obs_string.c_str());
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> words(begin, end);
        _axis.x = atof(words.at(1).c_str());
        _axis.y = atof(words.at(2).c_str());
        _axis.z = atof(words.at(3).c_str());
        Observable::update(atof(words.back().c_str()));
    }
    void BondPosProj::update(const simple::AtomPolymer &polymer){
        // assumes the molecule is a diatomic
        Observable::update(
            dot(subtract(
                    polymer.atoms.at(1).position,
                    polymer.atoms.at(0).position),
                _axis)
            / norm(subtract(
                    polymer.atoms.at(1).position,
                    polymer.atoms.at(0).position)));
    }
    void BondPosProj::update(const simple::AtomState &state){
        // assumes there is only one molecule in the state
        update(state.polymers.at(0));
    }
    BondVelProj::BondVelProj(Vector axis) :
        Observable("Bond Velocity Projection on " + vector_to_string(axis),
            "bvp", "1", ""),
        TimeLogObservable(),
        _axis(divide(axis, norm(axis))){}
    std::string BondVelProj::to_string() const {
        return "bondvelproj " + vector_to_string(_axis) + " " + print_value();
    }
    void BondVelProj::update(std::string obs_string){
        std::istringstream ss(obs_string.c_str());
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> words(begin, end);
        _axis.x = atof(words.at(1).c_str());
        _axis.y = atof(words.at(2).c_str());
        _axis.z = atof(words.at(3).c_str());
        Observable::update(atof(words.back().c_str()));
    }
    void BondVelProj::update(const simple::AtomPolymer &polymer){
        // assumes the molecule is a diatomic
        Observable::update(
            dot(subtract(
                    polymer.atoms.at(1).velocity,
                    polymer.atoms.at(0).velocity),
                _axis)
            / norm(subtract(
                    polymer.atoms.at(1).velocity,
                    polymer.atoms.at(0).velocity)));
    }
    void BondVelProj::update(const simple::AtomState &state){
        // assumes there is only one molecule in the state
        update(state.polymers.at(0));
    }
} // namespace diatomic
