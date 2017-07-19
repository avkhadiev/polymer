// 2017 Artur Avkhadiev
/*! \file rattle_integrator.cpp
*/
#include <vector>
#include <string>
#include <algorithm>                                /* std::fill */
#include <cmath>                                    /* pow */
#include "../include/ljpotential.h"
#include "../include/force_updater.h"
#include "../include/simulation.h"
#include "../include/molecule.h"
#include "../include/observable_container.h"
#include "../include/integrator.h"
#include "../include/verlet_integrator.h"
#include "../include/rattle_integrator.h"
RattleIntegrator::RattleIntegrator(ForceUpdater force_updater,
    double tol,
    double rvtol,
    double tiny,
    int maxiter,
    double *kinetic_energy_acc,
    double *neg_constraint_virial_acc) :
    VerletIntegrator(force_updater, NULL),
    _maxiter (maxiter),
    _tol (tol),
    _rvtol (rvtol),
    _tiny (tiny),
    _tol2 (2 * _tol),
    _rvtol2 (2 * _rvtol),
    _is_kinetic_energy_acc_set (false),
    _is_neg_constraint_virial_acc_set (false)
{
    if(kinetic_energy_acc != NULL){
        _is_kinetic_energy_acc_set = true;
        _kinetic_energy_acc = kinetic_energy_acc;
    }
    if(neg_constraint_virial_acc != NULL){
        _is_neg_constraint_virial_acc_set = true;
        _neg_constraint_virial_acc = neg_constraint_virial_acc;
    }
    _moving = {};
    _moved = {};
}
RattleIntegrator::~RattleIntegrator(){};
double RattleIntegrator::get_tol() const {
    return _tol;
}
double RattleIntegrator::get_rvtol() const {
    return _rvtol;
}
double RattleIntegrator::get_tiny() const {
    return _tiny;
}
double *RattleIntegrator::get_kinetic_energy_acc(){
    return _kinetic_energy_acc;
}
double *RattleIntegrator::get_neg_constraint_virial_acc(){
    return _neg_constraint_virial_acc;
}
void RattleIntegrator::set_tol(double tol){
    _tol = tol;
    _tol2 = _tol * 2.0;
}
void RattleIntegrator::set_rvtol(double rvtol){
    _rvtol = rvtol;
    _rvtol2 = _rvtol * 2.0;
}
void RattleIntegrator::set_tiny(double tiny){
    _tiny = tiny;
}
void RattleIntegrator::set_kinetic_energy_acc(double *kinetic_energy_acc){
    if (kinetic_energy_acc != NULL){
        _is_kinetic_energy_acc_set = true;
    }
    else {
        _is_kinetic_energy_acc_set = false;
    }
    _kinetic_energy_acc = kinetic_energy_acc;
}
void RattleIntegrator::set_neg_constraint_virial_acc(double *neg_constraint_virial_acc){
    if (neg_constraint_virial_acc != NULL){
        _is_neg_constraint_virial_acc_set = true;
    }
    else {
        _is_neg_constraint_virial_acc_set = false;
    }
    _neg_constraint_virial_acc = neg_constraint_virial_acc;
}
bool RattleIntegrator::_is_constraint_within_tol(Bond *bond, Molecule& molecule){
    bool is_contraint_within_tol = false;
    Atom *atom1 = &(molecule.atoms.at(bond->atom1));
    Atom *atom2 = &(molecule.atoms.at(bond->atom2));
    if (atom1->position.second != atom2->position.second){
        std::string err_msg = "is_contraint_within_tol: atomic position times ("
        + std::to_string(atom1->position.second)
        + " and "
        + std::to_string(atom2->position.second)
        + ") do not coincide";
        throw std::invalid_argument(err_msg);
    }
    else{
        double atomic_time = atom1->position.second;
        if (bond->current_length_sq.second != atomic_time) {
            update_bond_length(bond, molecule.atoms);
        }
        double difference_sq
            = std::abs(bond->current_length_sq.first - bond->fixed_length_sq);
        if (difference_sq < _tol2 * bond->fixed_length_sq){
            is_contraint_within_tol = true;
        }
        else {
            is_contraint_within_tol = false;
        }
    }
    return is_contraint_within_tol;
}
bool RattleIntegrator::_is_constraint_derivative_within_rvtol(Bond *bond,
    Molecule &molecule){
    bool is_contraint_derivative_within_rvtol = false;
    return is_contraint_derivative_within_rvtol;
}
Molecule RattleIntegrator::_move_correct_half_step(Molecule molecule_last_step,
    Molecule molecule_half_step_unconstrained){
    _set_up_correction_bookkeeping(molecule_half_step_unconstrained);
    bool still_correcting = true;
    int iter = 0;
    while (still_correcting && (iter < _maxiter)){
        /***********************************************************/
        iter += 1;
        still_correcting = false;
        /***********************************************************/
        // check the constrained bonds in the molecule
        try
        {
            check_bonds(molecule.bonds, molecule.atoms);
        }
        catch (std::invalid_argument &e)
        {
            throw;
        }
        // if all constrained bonds are legal, iterate over constrained bonds
        Bond *bond;
        for(int ib = 0; ib < molecule.nb; ++ib){
            if (_moved.at(ib) == false){
                continue;
            }
            else{
                bond = &(molecule.bonds.at(ib));
                if(_is_constraint_within_tol(bond, molecule) == false){
                    /**********************************************************/
                    still_correcting = true;
                    _moving.at(ib) = true;
                    /**********************************************************/
                }
            }
        }
    }
}
void RattleIntegrator::_move_correct_full_step(Molecule& molecule, bool calculate_observables){
    _set_up_correction_bookkeeping(molecule);
}
void RattleIntegrator::_zero_accumulators(){
    if(_is_kinetic_energy_acc_set){
        *_kinetic_energy_acc = 0.0;
    }
    if(_is_neg_constraint_virial_acc_set){
        *_neg_constraint_virial_acc = 0.0;
    }
}
void RattleIntegrator::_correct_accumulators(){
    if(_is_kinetic_energy_acc_set){
        *_kinetic_energy_acc = *_kinetic_energy_acc / 2.0;
    }
    if(_is_neg_constraint_virial_acc_set){
        *_neg_constraint_virial_acc
            = *_neg_constraint_virial_acc * 2.0 / 3.0 / _timestep;
    }
}
void RattleIntegrator::_set_up_correction_bookkeeping(Molecule& molecule){
    // all constrained bonds in the molecule have been moved
    _moved.resize(molecule.nb);
    std::fill(_moved.begin(), _moved.end(), true);
    // none of the bounds are moving
    _moving.resize(molecule.nb);
    std::fill(_moved.begin(), _moved.end(), false);
}
void RattleIntegrator::move(double timestep, State &state, bool calculate_observables){
    if(!is_time_consistent(state)){
        std::string err_msg = "move: state times not consistent before the move.";
        throw std::invalid_argument(err_msg);
    }
    else {
        // if necessary, zero the counters of accumulators
        if(calculate_observables){
            _zero_accumulators();
        }
        // set the time step and perform integration
        _set_timestep(timestep);
        // loop over molecules, perform verlet half-step and correction
        for (Molecule& molecule : state.molecules){
            _move_verlet_half_step(molecule);
            _move_correct_half_step(molecule);
        }
        // r(t + dt) has been calculated and corrected
        // now update forces in all atoms from f(t) to f(t + dt)
        // if calculate_observables = true, accumulators of observables that
        // should be updated in the force loop will be updated in this iteration
        _force_updater.update_forces(state, calculate_observables);
        // loop over molecules, perfrom verlet full-step and correction
        for (Molecule& molecule : state.molecules){
            _move_verlet_full_step(molecule, calculate_observables);
            _move_correct_full_step(molecule, calculate_observables);
        }
        // advance state time
        // now velocity, positions, and forces of all atoms have to be at t + dt
        state.time += _timestep;
        // multiply accumulators by whatever factors necessary
        if(calculate_observables){
            _correct_accumulators();
        }
        // check state for consistent time
        if(!is_time_consistent(state)){
            std::string err_msg = "move: state times not consistent after the move.";
            throw std::invalid_argument(err_msg);
        }
    }
}
