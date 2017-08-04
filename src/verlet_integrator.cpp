// 2017 Artur Avkhadiev
/*! \file verlet_integrator.cpp
*/
#include <vector>
#include <string>
#include <cmath>              /* pow */
#include "../include/ljpotential.h"
#include "../include/force_updater.h"
#include "../include/simulation.h"
#include "../include/molecule.h"
#include "../include/observable_container.h"
#include "../include/integrator.h"
#include "../include/verlet_integrator.h"
VerletIntegrator::VerletIntegrator(ForceUpdater force_updater,
    double *kinetic_energy_acc) :
    _force_updater (force_updater),
    _timestep (0.001),
    _halfstep (_timestep * 0.5),
    _is_kinetic_energy_acc_set (false) {
    if(kinetic_energy_acc != NULL){
        _is_kinetic_energy_acc_set = true;
        _kinetic_energy_acc = kinetic_energy_acc;
    }
}
VerletIntegrator::~VerletIntegrator(){};
ForceUpdater& VerletIntegrator::get_force_updater() {
    return _force_updater;
}
double *VerletIntegrator::get_kinetic_energy_acc(){
    return _kinetic_energy_acc;
}
void VerletIntegrator::_set_timestep(double timestep){
    _timestep = timestep;
    _halfstep = timestep * 0.5;
}
void VerletIntegrator::set_force_updater(ForceUpdater force_updater) {
    _force_updater = force_updater;
}
void VerletIntegrator::set_kinetic_energy_acc(double *kinetic_energy_acc){
    if (kinetic_energy_acc != NULL){
        _is_kinetic_energy_acc_set = true;
    }
    else {
        _is_kinetic_energy_acc_set = false;
    }
    _kinetic_energy_acc = kinetic_energy_acc;
}
void VerletIntegrator::_zero_accumulators(){
    if(_is_kinetic_energy_acc_set){
        *_kinetic_energy_acc = 0.0;
    }
}
void VerletIntegrator::_correct_accumulators(){
    if(_is_kinetic_energy_acc_set){
        *_kinetic_energy_acc = *_kinetic_energy_acc / 2.0;
    }
}
/**
* ::STATE STRUCTURES (LEGACY)
*/
Molecule VerletIntegrator::_move_verlet_half_step(Molecule molecule) {
    // loop over all atoms
    Atom *atom;
    for(int j = 0; j < molecule.na; ++j){
        atom = &(molecule.atoms.at(j));
        // check that force in atom is calculated at the same time as its
        // position
        if(atom->force.second != atom->position.second){
            std::string err_msg = "move_verlet_half_step: force in an atom recorded for time " + std::to_string(atom->force.second) + ", but expected time " + std::to_string(atom->position.second) + ".";
            throw std::invalid_argument(err_msg);
        }
        else{
            // v(t + 0.5 dt) = v(t) + 1/2 dt a(t)
            Vector acceleration = divide(atom->force.first, atom->mass);
            atom->velocity.first += multiply(acceleration, _halfstep);
            // will not add halfstep to velocity time; will add full step
            // later --- this is to avoid floating pointing arithmetic
            // issues.
            // r(t + dt) = r(t) + dt * v(t + 0.5 dt)
            atom->position.first += multiply(atom->velocity.first, _timestep);
            atom->position.second += _timestep;
        }
    }
    // now all atoms have their velocities at half step and their
    // positions at full step; we can update bond lengths.
    for (int ib = 0; ib < molecule.nb; ++ib){
        try {
            update_bond_length(&(molecule.bonds.at(ib)), molecule.atoms);
        }
        catch (std::invalid_argument &e) {
            throw;
        }
    }
    return molecule;
}
Molecule VerletIntegrator::_move_verlet_full_step(Molecule molecule, bool calculate_observables) {
    // loop over all atoms
    Atom *atom;
    for(int j = 0; j < molecule.na; ++j){
        atom = &(molecule.atoms.at(j));
        // v(t + dt) = v(t + 0.5 dt) + 0.5 dt * a(t + dt)
        Vector acceleration = divide(atom->force.first, atom->mass);
        atom->velocity.first += multiply(acceleration, _halfstep);
        // adding full time step (half step was not applied) --- this is
        // to avoid floating pointing arithmetic issues.
        atom->velocity.second += _timestep;
        // update kinetic energy accumulator if necessary and possible
        if(calculate_observables){
            if(_is_kinetic_energy_acc_set){
                double k_increase = atom->mass*normsq(atom->velocity.first);
                *_kinetic_energy_acc += k_increase;
            }
        }
    }
    return molecule;
}
void VerletIntegrator::move(double timestep, State &state, bool calculate_observables){
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
        Molecule *molecule;
        // loop over molecules, perform verlet half-step
        for (int im = 0; im < state.nm; ++im){
            molecule = &(state.molecules.at(im));
            *molecule = _move_verlet_half_step(*molecule);
        }
        // now update forces in all atoms from f(t) to f(t + dt)
        // if calculate_observables = true, accumulators of observables that
        // should be updated in the force loop will be updated in this iteration
        _force_updater.update_forces(state, calculate_observables);
        // loop over molecules, perfrom verlet full-step
        for (int im = 0; im < state.nm; ++im){
            molecule = &(state.molecules.at(im));
            *molecule = _move_verlet_full_step(*molecule, calculate_observables);
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
/**
* SIMPLE::STATE STRUCTURES
*/
simple::AtomPolymer VerletIntegrator::_move_verlet_half_step(simple::AtomPolymer molecule) {
    int na = molecule.nb() + 1;
    double m = molecule.m();
    // loop over all atoms
    for(int j = 0; j < na; ++j){
        simple::Atom& atom = molecule.atoms.at(j);
        // v(t + 0.5 dt) = v(t) + 1/2 dt a(t)
        Vector acceleration = divide(atom.force, m);
        atom.velocity += multiply(acceleration, _halfstep);
        atom.position += multiply(atom.velocity, _timestep);
    }
    // now all atoms have their velocities at half step and their
    // positions at full step
    return molecule;
}
simple::AtomPolymer VerletIntegrator::_move_verlet_full_step(simple::AtomPolymer molecule,
    bool calculate_observables) {
    int na = molecule.nb() + 1;
    double m = molecule.m();
    // loop over all atoms
    for(int j = 0; j < na; ++j){
        simple::Atom& atom = molecule.atoms.at(j);
        // v(t + dt) = v(t + 0.5 dt) + 0.5 dt * a(t + dt)
        Vector acceleration = divide(atom.force, m);
        atom.velocity += multiply(acceleration, _halfstep);
        // update kinetic energy accumulator if necessary and possible
        if(calculate_observables){
            if(_is_kinetic_energy_acc_set){
                double k_increase = m * normsq(atom.velocity);
                *_kinetic_energy_acc += k_increase;
            }
        }
    }
    return molecule;
}
void VerletIntegrator::move(double timestep, simple::AtomState& state,
    bool calculate_observables){
        // if necessary, zero the counters of accumulators
        if(calculate_observables){
            _zero_accumulators();
        }
        simple::AtomPolymer polymer;
        // set the time step and perform integration
        _set_timestep(timestep);
        // loop over molecules, perform verlet half-step
        for (int im = 0; im < state.nm(); ++im){
            // passing molecule to the helper function by value is important
            // for inheritance reasons --- RATTLE, which derives from Verlet,
            // will use both old information from time (t) as well as new
            // information at (t + dt) from unconstrained Verlet integration
            polymer = state.polymers.at(im);
            state.polymers.at(im) = _move_verlet_half_step(polymer);
        }
        // now update forces in all atoms from f(t) to f(t + dt)
        // if calculate_observables = true, accumulators of observables that
        // should be updated in the force loop will be updated in this iteration
        _force_updater.update_forces(state, calculate_observables);
        // loop over molecules, perfrom verlet full-step
        for (int im = 0; im < state.nm(); ++im){
            // passing molecule to the helper function by value is important
            // for inheritance reasons --- RATTLE, which derives from Verlet,
            // will use both old information from time (t) as well as new
            // information at (t + dt) from unconstrained Verlet integration
            polymer = state.polymers.at(im);
            state.polymers.at(im) = _move_verlet_full_step(polymer,
                calculate_observables);
        }
        // advance state time
        // now velocity, positions, and forces of all atoms have to be at t + dt
        state.advance_time(timestep);
        // multiply the accumulators by whatever factors necessary
        if(calculate_observables){
            _correct_accumulators();
        }
}
