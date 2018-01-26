// 2017 Artur Avkhadiev
/*! \file polymer_observables.cpp
*/
#include "../include/polymer_observables.h"
namespace polymer{
    /**************************************************************************
    * Kinetic Energy
    **************************************************************************/
    KE::KE( bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    observable::KE(calculate_mean, calculate_error, print_inst_val, e_format){
        _name.full += " (Polymer)";
        _name.abridged += "p";
        _name.latex += "^{\\mathrm{pol}}";
    }
    void KE::update(const simple::AtomPolymer& molecule){
        double mass = molecule.m();
        for(const simple::Atom& atom : molecule.atoms){
            observable::KE::update(atom, mass);
        }
    }
    void KE::update(const simple::AtomState &state){
        for(const simple::AtomPolymer& polymer : state.polymers){
            update(polymer);
        }
    }
    /**************************************************************************
    * LJ Potential Energy
    **************************************************************************/
    PE::PE(bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    observable::PE(calculate_mean, calculate_error, print_inst_val, e_format){
        _name.full += " (Polymer)";
        _name.abridged += "p";
        _name.latex += "^{\\mathrm{pol}}";
    }
    /**************************************************************************
    * Kinetic Temperature
    **************************************************************************/
    KineticTemperature::KineticTemperature(
        bool remove_linear_momentum,
        bool remove_angular_momentum,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    observable::KineticTemperature(
        remove_linear_momentum, remove_angular_momentum,
        calculate_mean, calculate_error, print_inst_val, e_format),
    _ndof(std::max(0,
            3 * ((simple::BaseState::nm())
                - (int)_remove_linear_momentum
                - (int)_remove_angular_momentum)
            + 2 * (simple::BaseState::nm() * simple::BasePolymer::nb())
        ))
    {
        _name.full += " (Polymer)";
        _name.abridged += "p";
        _name.latex += "^{\\mathrm{pol}}";
        if (DEBUG) {
            fprintf(stderr, "%s %zu\n", "KineticTemperature: ndof is", _ndof);
        }
    }
    void KineticTemperature::update(const simple::AtomPolymer& molecule){
        double mass = molecule.m();
        for(const simple::Atom& atom : molecule.atoms){
            observable::KineticTemperature::update(atom, mass, ndof());
        }
    }
    void KineticTemperature::update(const simple::AtomState &state){
        for(const simple::AtomPolymer& polymer : state.polymers){
            update(polymer);
        }
    }
    /**************************************************************************
    * Virial from Potential Forces
    **************************************************************************/
    Virial::Virial(bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    observable::Virial(calculate_mean,calculate_error,print_inst_val, e_format){
        _name.full += " (Polymer)";
        _name.abridged += "p";
        _name.latex += "^{\\mathrm{pol}}_{\\mathrm{LJ}}";
    }
    /**************************************************************************
    * Virial from Constraint Forces
    **************************************************************************/
    ConstraintVirial::ConstraintVirial(bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    observable::Virial(calculate_mean,calculate_error,print_inst_val, e_format){
        _name.full += " from Constraint Forces (Polymer)";
        _name.abridged += "cp";
        _name.latex += "^{\\mathrm{pol}}_{\\mathrm{c}}";
        // constraint virial is calculated by RATTLE
        _update_time = observable::update_time_t::INTEGRATION_STEP;
    }
    /**************************************************************************
    * Linear Momentum Component
    **************************************************************************/
    LinMomComponent::LinMomComponent(Vector component,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    observable::LinMomComponent(component,
        calculate_mean, calculate_error,
        print_inst_val, e_format){
            _name.full += " (Polymer)";
            _name.abridged += "p";
            _name.latex += "^{\\mathrm{pol}}";
    }
    void LinMomComponent::update(const simple::AtomPolymer& molecule){
        double mass = molecule.m();
        for(const simple::Atom& atom : molecule.atoms){
            observable::LinMomComponent::update(atom, mass);
        }
    }
    void LinMomComponent::update(const simple::AtomState &state){
        for(const simple::AtomPolymer& polymer : state.polymers){
            update(polymer);
        }
    }
    /**************************************************************************
    * Angular Momentum Component
    **************************************************************************/
    AngMomComponent::AngMomComponent(Vector component,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    observable::AngMomComponent(component,
        calculate_mean, calculate_error,
        print_inst_val, e_format){
            _name.full += " (Polymer)";
            _name.abridged += "p";
            _name.latex += "^{\\mathrm{pol}}";
    }
    void AngMomComponent::update(const simple::AtomPolymer& molecule){
        double mass = molecule.m();
        for(const simple::Atom& atom : molecule.atoms){
            observable::AngMomComponent::update(atom, mass);
        }
    }
    void AngMomComponent::update(const simple::AtomState &state){
        for(const simple::AtomPolymer& polymer : state.polymers){
            update(polymer);
        }
    }
    /**************************************************************************
    * RCM Component
    **************************************************************************/
//    RCMComponent::RCMComponent(Vector component,
//        size_t molecule_index,
//        bool calculate_mean,
//        bool calculate_error,
//        bool print_inst_val,
//        bool e_format) :
//    observable::RCMComponent(component, molecule_index,
//        calculate_mean, calculate_error,
//        print_inst_val, e_format){
//            _name.full += " (Polymer)";
//            _name.abridged += "p";
//            _name.latex += "(\\mathrm{pol})";
//    }
//    Vector RCMComponent::rcm(const simple::AtomState& state){
//        Vector rcm = vector(0.0, 0.0, 0.0);
//        if(molecule_index() > simple::BaseState::nm()){
//            fprintf(stderr, "%s\n", "Polymer RCMComponent: molecule index exceeds number of polymer molecules in the state. Returning a zero vector.");
//        }
//        else {
//            Vector ri;
//            size_t natoms = simple::BasePolymer::nb() + 1;
//            for(size_t ia = 0; ia < natoms; ++ia){
//                ri = state.polymers.at(molecule_index()).atoms.at(ia).position;
//                add(rcm, ri);
//            }
//            divide(rcm, natoms);
//        }
//        return rcm;
//    }
//    void RCMComponent::update(const simple::AtomState& state){
//        _update(rcm(state), component());
//    }
} // namespace polymer
