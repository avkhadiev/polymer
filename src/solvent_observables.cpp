// 2017 Artur Avkhadiev
/*! \file solvent_observables.cpp
*/
#include "../include/solvent_observables.h"
#include "../include/simple_solvent.h"
namespace solvent{
    /**************************************************************************
    * Kinetic Energy
    **************************************************************************/
    KE::KE( bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Solvent Kinetic Energy",
                .abridged = "ske",
                .latex = "\\T_{\\mathrm{s}}",
                .units = "\\varepsilon"},
                observable::INTEGRATION_STEP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format){}
    void KE::update(const simple::Solvent &molecule){
        value += simple::Solvent::m() * normsq(molecule.v()) / 2;
    }
    void KE::update(const simple::AtomState &state){
        for(const simple::Solvent& solvent : state.solvents){
            update(solvent);
        }
    }
    /**************************************************************************
    * LJ Potential Energy
    **************************************************************************/
    PE::PE( bool shifted,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Solvent Potential",
                .abridged = "svlj",
                .latex = "",
                .units = "\\varepsilon"},
                observable::FORCE_LOOP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format),
    _is_shifted(shifted){
        if (_is_shifted){
            _name.full += ", Shifted";
            _name.abridged += "_shft";
            _name.latex = "U^{\\mathrm{s shft}}_{\\mathrm{LJ}}";
        }
        else {
            _name.full += ", Unshifted";
            _name.abridged += "_unshft";
            _name.latex = "U^{\\mathrm{s unshft}}_{\\mathrm{LJ}}";
        }
    };
    /**
    * Virial
    */
    Virial::Virial(bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Solvent Internal Virial",
                .abridged = "sw",
                .latex = "\\mathcal{W}_{\\mathrm{s}}",
                .units = "\\varepsilon"},
                observable::INTEGRATION_STEP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format){};
    /**
    * Momentum Component
    */
    MomComponent::MomComponent(Component component,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Solvent Momentum",
                .abridged = "sp",
                .latex = "",
                .units = ""},
                observable::MAIN_LOOP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format),
    _component(component)
    {
        switch (_component) {
            case X:
                _name.full += " (X)";
                _name.abridged += "_x";
                _name.latex = "p^{\\mathrm{s}}_{x}";
                break;
            case Y:
                _name.full += " (Y)";
                _name.abridged += "_y";
                _name.latex = "p^{\\mathrm{s}}_{y}";
                break;
            case Z:
                _name.full += " (Z)";
                _name.abridged += "_z";
                _name.latex = "p^{\\mathrm{s}}_{z}";
                break;
            default:
                std::string err_msg
                    = "solvent momentum: component argument is  = 1, 2, or 3!";
                throw std::invalid_argument(err_msg);
        }
    };
    void MomComponent::_update(const simple::Solvent &molecule,
        Component component){
        switch (_component) {
        case X:
            value += molecule.m() * molecule.v().x;
            break;
        case Y:
            value += molecule.m() * molecule.v().y;
            break;
        case Z:
            value += molecule.m() * molecule.v().z;
        }
    }
    void MomComponent::update(const simple::Solvent &molecule){
        _update(molecule, _component);
    }
    void MomComponent::update(const simple::AtomState &state){
        for(const simple::Solvent& solvent : state.solvents){
            update(solvent);
        }
    }
    /**
    * Momentum Magnitude Squared
    */
    MomMagSq::MomMagSq(bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
    Observable({.full = "Solvent Momentum Magnitude Squared",
                .abridged = "spsq",
                .latex = "\\left\\lVert \\vec{p} \\right\\rVert^2",
                .units = ""},
                observable::MAIN_LOOP,
                {.mean = calculate_mean || calculate_error,
                 .meansq = calculate_error},
                 print_inst_val, e_format){};
    void MomMagSq::update(const simple::Solvent &molecule){
        value += molecule.m() * normsq(molecule.v());
    }
    void MomMagSq::update(const simple::AtomState &state){
        for(const simple::Solvent& solvent : state.solvents){
            update(solvent);
        }
    }
     /**
     * Kinetic Temperature
     */
     KinTemp::KinTemp(bool calculate_mean,
         bool calculate_error,
         bool print_inst_val,
         bool e_format) :
     Observable({.full = "Kinetic Temperature",
                 .abridged = "sT",
                 .latex = "\\mathcal{T}^{\\mathrm{s}}_{\\mathrm{kin}}",
                 .units = "\\varepsilon"},
                 observable::MAIN_LOOP,
                 {.mean = calculate_mean || calculate_error,
                  .meansq = calculate_error},
                  print_inst_val, e_format),
    _nc(3){};
    void KinTemp::update(const simple::AtomState &state){
        size_t ndof = 3 * (state.nsolvents() - _nc);
        value = 0.0;
        for(const simple::Solvent& solvent : state.solvents){
            value += normsq(solvent.v());
        }
        value = value * state.solvent_mass() / ndof;
    }
} // namespace solvent
