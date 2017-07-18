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
    _is_kinetic_energy_set (false) {
    if(kinetic_energy_acc != NULL){
        _is_kinetic_energy_set = true;
        _kinetic_energy_acc = kinetic_energy_acc;
    }
}
VerletIntegrator::~VerletIntegrator(){};
ForceUpdater& VerletIntegrator::get_force_updater() {
    return _force_updater;
}
void VerletIntegrator::set_force_updater(ForceUpdater force_updater) {
    _force_updater = force_updater;
}
void VerletIntegrator::_move_verlet_half_step(double timestep, State& state) {
    // assumes forces have been calculated at time t
    double halfstep = 0.5 * timestep;
    double old_time = state.time;
    Molecule *molecule;
    Atom *atom;
    // loop over all molecules
    for(int i = 0; i < state.nm; ++i){
        molecule = &(state.molecules.at(i));
        // loop over all atoms
        for(int j = 0; j < molecule->na; ++j){
            atom = &(molecule->atoms.at(j));
            // check that force in atom is calculated at time t
            if(atom->force.second != old_time){
                std::string err_msg = "move_verlet_half_step: force in an atom recorded for time " + std::to_string(atom->force.second) + ", but expected time " + std::to_string(old_time) + ".";
                throw std::invalid_argument(err_msg);
            }
            else{
                // v(t + 0.5 dt) = v(t) + 1/2 dt a(t)
                Vector acceleration = divide(atom->force.first, atom->mass);
                atom->velocity.first += multiply(acceleration, halfstep);
                // will not add halfstep to velocity time; will add full step
                // later --- this is to avoid floating pointing arithmetic
                // issues.
                // r(t + dt) = r(t) + dt * v(t + 0.5 dt)
                atom->position.first += multiply(atom->velocity.first, timestep);
                atom->position.second += timestep;
            }
        }
    }
    // now all atoms have their velocities at half step and their
    // positions at full step
}
void VerletIntegrator::_move_verlet_full_step(double timestep, State& state, bool calculate_observables) {
    double halfstep = 0.5 * timestep;
    // now update forces in all atoms from f(t) to f(t + dt)
    // if calculate_observables = true, accumulators of observables that should
    // be updated in the force loop will be updated in this iteration
    _force_updater.update_forces(state, calculate_observables);
    // run the force loop
    Molecule *molecule;
    Atom *atom;
    // loop over all molecules
    for(int i = 0; i < state.nm; ++i){
        molecule = &(state.molecules.at(i));
        // loop over all atoms
        for(int j = 0; j < molecule->na; ++j){
            atom = &(molecule->atoms.at(j));
            // v(t + dt) = v(t + 0.5 dt) + 0.5 dt * a(t + dt)
            Vector acceleration = divide(atom->force.first, atom->mass);
            atom->velocity.first += multiply(acceleration, halfstep);
            // adding full time step (half step was not applied) --- this is
            // to avoid floating pointing arithmetic issues.
            atom->velocity.second += timestep;
            // update kinetic energy accumulator if necessary and possible
            if(calculate_observables){
                if(_is_kinetic_energy_set){
                    double k_increase = atom->mass*normsq(atom->velocity.first);
                    *_kinetic_energy_acc += k_increase;
                }
            }
        }
    }
    // now velocity, positions, and forces of all atoms have to be at t + dt
    // update state time
    state.time += timestep;
}
void VerletIntegrator::move(double timestep, State &state, bool calculate_observables){
    if(!is_time_consistent(state)){
        std::string err_msg = "move: state times not consistent before the move.";
        throw std::invalid_argument(err_msg);
    }
    else {
        // if necessary, zero the counters of accumulators
        if(calculate_observables){
            if(_is_kinetic_energy_set){
                *_kinetic_energy_acc = 0.0;
            }
        }
        _move_verlet_half_step(timestep, state);
        _move_verlet_full_step(timestep, state, calculate_observables);
        if(calculate_observables){
            if(_is_kinetic_energy_set){
                *_kinetic_energy_acc = *_kinetic_energy_acc / 2.0;
            }
        }
        // check state for consistent time
        if(!is_time_consistent(state)){
            std::string err_msg = "move: state times not consistent after the move.";
            throw std::invalid_argument(err_msg);
        }
    }
}
