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
        value = 0.0;
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
        value = 0.0;
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
        value = 0.0;
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
        value = 0.0;
        for(const simple::AtomPolymer& polymer : state.polymers){
            update(polymer);
        }
    }
    /**************************************************************************
    * RCM Component
    **************************************************************************/
    RCMComponent::RCMComponent(Vector component,
        size_t molecule_index,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    observable::RCMComponent(component, molecule_index,
        calculate_mean, calculate_error,
        print_inst_val, e_format){
            _name.full += " (Polymer)";
            _name.abridged += "p";
            _name.latex += "(\\mathrm{pol})";
    }
    Vector RCMComponent::rcm(const simple::AtomState& state){
        Vector rcm = vector(0.0, 0.0, 0.0);
        if(molecule_index() > simple::BaseState::nm()){
            fprintf(stderr, "%s\n", "Polymer RCMComponent: molecule index exceeds number of polymer molecules in the state. Returning a zero vector.");
        }
        else {
            rcm = state.get_polymers().at(molecule_index()).rcm();
            //fprintf(stderr, "%s\n", "RCM");
            //fprintf(stderr, "%s\n", vector_to_string(rcm).c_str());
        }
        return rcm;
    }
    void RCMComponent::update(const simple::AtomState& state){
        observable::RCMComponent::_update(rcm(state), component());
    }
    BondLength::BondLength(
        size_t molecule_index,
        size_t bond_index,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    geodesic::LinkObservable({.full = "Link Length",
                .abridged = "omega_norm",
                .latex = "\\left|{\\hat{\\Omega}_{k}}\\right|",
                .units = "1"},
                bond_index + 1,
                observable::MAIN_LOOP,
                calculate_mean, calculate_error,
                print_inst_val, e_format),
    _molecule_index(molecule_index){
        if(_molecule_index > simple::BaseState::nm()){
            fprintf(stderr, "%s\n", "Polymer BondLength: molecule index exceeds number of polymer molecules in the state.");
        }
    }
    void BondLength::update(const simple::AtomState& state){
        simple::Atom atom1 = state.get_polymers().at(_molecule_index).get_atoms().at(_link_number - 1);
        simple::Atom atom2 = state.get_polymers().at(_molecule_index).get_atoms().at(_link_number);
        value = norm(subtract(atom1.position, atom2.position));
    }
    /**************************************************************************
    * Polymer Atom Position Component
    **************************************************************************/
    PolAtomPosComponent::PolAtomPosComponent(
    Vector component,
        size_t molecule_index,
        size_t atom_index,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format ):
    Observable({.full = "Polymer Atom Position",
                .abridged = "rpol_atom",
                .latex = "\\vec{r}",
                .units = "\\sigma"},
                observable::MAIN_LOOP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format),
    _component(component),
    _molecule_index(molecule_index),
    _atom_index(atom_index)
    {
        // check indices for existence
        if (_molecule_index > simple::BaseState::nm()){
            fprintf(stderr, "%s\n", "Polymer Atom Position: molecule index exceeds number of polymer molecules in the state.");
        }
        if (_atom_index > simple::BasePolymer::nb() + 1){
            fprintf(stderr, "%s\n", "Polymer Atom Position: atom index exceeds number of atoms in the polymer.");
        }
        amend_names();
    }
    void PolAtomPosComponent::amend_names(){
        std::string mol_index_str = std::to_string(_molecule_index + 1);
        std::string atom_index_str = std::to_string(_atom_index + 1);
        std::string component_str;
        if (_component == vector(1.0, 0.0, 0.0)){
            component_str = "x";
        }
        else if (_component == vector(0.0, 1.0, 0.0)){
            component_str = "y";
        }
        else if (_component == vector(0.0, 0.0, 1.0)){
            component_str = "z";
        }
        else {
            component_str
                = std::to_string(_component.x) + ","
                + std::to_string(_component.y) + ","
                + std::to_string(_component.z);
        }
        _name.full += " " + mol_index_str + " " + atom_index_str
            + " (" + component_str + ")";
        _name.abridged += "_" + mol_index_str
                        + "_" +  atom_index_str
                        + "_" + component_str;
        _name.latex +=
            "^{(" + mol_index_str + ", " + atom_index_str + ")}_{"
            + component_str + "}";
    }
    void PolAtomPosComponent::_update(Vector r, Vector component){
        value = dot(r, component);
    }
    void PolAtomPosComponent::update(const simple::AtomState& state){
        simple::Atom atom = state.get_polymers().at(_molecule_index).get_atoms().at(_atom_index);
        _update(atom.position, _component);
    }
} // namespace polymer
