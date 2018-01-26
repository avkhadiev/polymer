// 2017 Artur Avkhadiev
/*! \file general_observables.cpp
*/
#include <algorithm>    /**> std::max */
#include "../include/general_observables.h"
namespace observable{
    Vector x_vec = vector(1.0, 0.0, 0.0);
    Vector y_vec = vector(0.0, 1.0, 0.0);
    Vector z_vec = vector(0.0, 0.0, 1.0);
    /**************************************************************************
    * TIME
    **************************************************************************/
    Time::Time() :
    Observable({.full = "Time",
                .abridged = "t",
                .latex = "\\tau",
                .units = "\\tau_{\\mathrm{LJ}}"},
                observable::MAIN_LOOP,
                {.mean = false,
                 .meansq = false},
                 false, false){}
    void Time::update(const simple::AtomState& state){
        value = state.time();
    }
    /**************************************************************************
    * Kinetic Energy
    **************************************************************************/
    KE::KE( bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Kinetic Energy",
                .abridged = "ke",
                .latex = "\\mathcal{K}",
                .units = "\\varepsilon"},
                observable::INTEGRATION_STEP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format){}
    void KE::update(const simple::Atom& atom, double mass){
        value += 0.5 * mass * normsq(atom.velocity);
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
    Observable({.full = "Kinetic Temperature",
                .abridged = "temp_kin",
                .latex = "\\mathcal{T}_{\\mathrm{kin}}",
                .units = "\\varepsilon"},
                observable::MAIN_LOOP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format),
    _remove_linear_momentum(remove_linear_momentum),
    _remove_angular_momentum(remove_angular_momentum),
    _ndof(std::max(0,
            3 * ((simple::BaseState::nm() + simple::BaseState::nsolvents())
                - (int)_remove_linear_momentum
                - (int)_remove_angular_momentum)
            + 2 * (simple::BaseState::nm() * simple::BasePolymer::nb())
        ))
    {}
    void KineticTemperature::update(const simple::Atom& atom,
        double mass, size_t ndof){
        if (ndof > 0) {
            value += (mass * normsq(atom.velocity) / ((double)_ndof));
        }
        else {
            value += 0.0;
        }
    }
    /**************************************************************************
    * LJ Potential Energy
    **************************************************************************/
    PE::PE(bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Potential Energy",
                .abridged = "vlj",
                .latex = "\\mathcal{U}_{\\mathrm{LJ}}",
                .units = "\\varepsilon"},
                observable::FORCE_LOOP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format){};
    /**
    * Virial
    */
    Virial::Virial(bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Virial",
                .abridged = "w",
                .latex = "\\mathcal{W}",
                .units = "\\varepsilon"},
                observable::INTEGRATION_STEP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format){};
    /**
    * Linear Momentum Component
    */
    LinMomComponent::LinMomComponent(Vector component,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Linear Momentum",
                .abridged = "p",
                .latex = "p",
                .units = "M\\sigma\\tau_{\\mathrm{LJ}}^{-1}"},
                observable::MAIN_LOOP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format),
    _component(divide(component, norm(component)))
    {
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
                = std::to_string(component.x) + ","
                + std::to_string(component.y) + ","
                + std::to_string(component.z);
        }
        _name.full += " (" + component_str + ")";
        _name.abridged += "_" + component_str;
        _name.latex += "_{" + component_str + "}";
    };
    void LinMomComponent::_update(Vector momentum, Vector component){
        value += dot(momentum, component);
    }
    void LinMomComponent::update(const simple::Atom &atom, double mass){
        Vector momentum = multiply(atom.velocity, mass);
        _update(momentum, _component);
    }
    /**
    * Angular Momentum Component
    */
    AngMomComponent::AngMomComponent(Vector component,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Angular Momentum",
                .abridged = "L",
                .latex = "L",
                .units = "M\\sigma^{2}\\tau_{\\mathrm{LJ}}^{-1}"},
                observable::MAIN_LOOP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format),
    _component(divide(component, norm(component)))
    {
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
                = std::to_string(component.x) + ","
                + std::to_string(component.y) + ","
                + std::to_string(component.z);
        }
        _name.full += " (" + component_str + ")";
        _name.abridged += "_" + component_str;
        _name.latex += "_{" + component_str + "}";
    };
    void AngMomComponent::_update(Vector momentum, Vector component){
        value += dot(momentum, component);
    }
    void AngMomComponent::update(const simple::Atom &atom, double mass){
        Vector momentum = cross(atom.position, multiply(atom.velocity, mass));
        _update(momentum, _component);
    }
    /**
    * Angular Momentum Component
    */
    RCMComponent::RCMComponent(Vector component,
        size_t molecule_index,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "COM Position",
                .abridged = "rcm",
                .latex = "R^{" + std::to_string(molecule_index) + "}",
                .units = "\\sigma"},
                observable::MAIN_LOOP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format),
    _component(divide(component, norm(component))),
    _molecule_index(molecule_index)
    {
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
                = std::to_string(component.x) + ","
                + std::to_string(component.y) + ","
                + std::to_string(component.z);
        }
        _name.full += " (" + component_str + ")";
        _name.abridged += "_" + component_str;
        _name.latex += "_{\\mathrm{CM}" + component_str + "}";
    }
    void RCMComponent::_update(Vector rcm, Vector component){
        value = dot(rcm, component);
    }
} // namespace observable
