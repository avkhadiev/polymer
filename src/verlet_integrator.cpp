// 2017 Artur Avkhadiev
/*! \file verlet_integrator.cpp
*/
#include <vector>
#include <string>
#include <cmath>              /* pow   */
#include <math.h>             /* round */
#include "../include/verlet_integrator.h"
VerletIntegrator::VerletIntegrator(ForceUpdater& force_updater,
    Observable *ke,
    double box) :
    _force_updater (force_updater),
    _timestep (0.001),
    _halfstep (_timestep * 0.5) {
    _ke.ptr = ke;
    _mirror_image.box = box;
    if(_ke.ptr == NULL){
        _ke.is_set = false;
    }
    else{
        _ke.is_set = true;
    }
    if(box == 0.0){
        _mirror_image.is_set = false;
    }
    else{
        _mirror_image.is_set = true;
    }
}
VerletIntegrator::~VerletIntegrator(){};
ForceUpdater& VerletIntegrator::get_force_updater() {
    return _force_updater;
}
void VerletIntegrator::_set_timestep(double timestep){
    _timestep = timestep;
    _halfstep = timestep * 0.5;
}
void VerletIntegrator::set_ke(Observable *ptr){
    if (ptr != NULL){
        _ke.is_set = true;
    }
    else {
        _ke.is_set = false;
    }
    _ke.ptr = ptr;
}
void VerletIntegrator::_zero_accumulators(){
    if(_ke.is_set) _ke.ptr->zero();
}
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
            if(_ke.is_set){
                _ke.ptr->update(atom);
            }
        }
    }
    return molecule;
}
// if mirror image is set, uses box to apply the mirror image convention
// outputs the adjusted position.
// uses function round() from math.h instead of ANINT (same behavior)
Vector VerletIntegrator::_get_mirror_image(Vector r){
    if (_mirror_image.is_set){
        double box = _mirror_image.box;
        r.x = r.x - box * round( r.x / box );
        r.y = r.y - box * round( r.y / box );
        r.z = r.z - box * round( r.z / box );
    }
    return r;
}
simple::Solvent
    VerletIntegrator::_move_verlet_half_step(simple::Solvent molecule) {
    double m = simple::Solvent::m();
    // v(t + 0.5 dt) = v(t) + 1/2 dt a(t)
    Vector acceleration = divide(molecule.f(), m);
    molecule.set_v( add(molecule.v(), multiply(acceleration, _halfstep)) );
    molecule.set_r( add(molecule.r(), multiply(molecule.v(), _timestep)) );
    // mirror image adjustment
    molecule.set_r( _get_mirror_image(molecule.r()) );
    return molecule;
}
simple::Solvent
    VerletIntegrator::_move_verlet_full_step(simple::Solvent molecule,
    bool calculate_observables) {
    double m = simple::Solvent::m();
    // v(t + dt) = v(t + 0.5 dt) + 0.5 dt * a(t + dt)
    Vector acceleration = divide(molecule.f(), m);
    molecule.set_v( add(molecule.v(), multiply(acceleration, _halfstep)) );
    // update kinetic energy accumulator if necessary and possible
    if(calculate_observables){
        if(_ke.is_set){
            _ke.ptr->update(molecule.atom);
        }
    }
    return molecule;
}
void VerletIntegrator::move(double timestep, simple::AtomState& state,
    bool calculate_observables){
        // if necessary, zero the counters of accumulators
        if(calculate_observables) _zero_accumulators();
        simple::AtomPolymer polymer;
        simple::Solvent solvent;
        // set the time step and perform integration
        _set_timestep(timestep);
        // loop over molecules, perform verlet half-step
        // polymer loop
        for (int im = 0; im < state.nm(); ++im){
            // passing molecule to the helper function by value is important
            // for inheritance reasons --- RATTLE, which derives from Verlet,
            // will use both old information from time (t) as well as new
            // information at (t + dt) from unconstrained Verlet integration
            polymer = state.polymers.at(im);
            state.polymers.at(im) = _move_verlet_half_step(polymer);
        }
        // solvent loop
        for (int is = 0; is < state.nsolvents(); ++is){
            solvent = state.solvents.at(is);
            state.solvents.at(is) = _move_verlet_half_step(solvent);
        }
        // now update forces in all atoms from f(t) to f(t + dt)
        // if calculate_observables = true, accumulators of observables that
        // should be updated in the force loop will be updated in this iteration
        _force_updater.update_forces(state, calculate_observables);
        // loop over molecules, perfrom verlet full-step
        // polymer loop
        for (int im = 0; im < state.nm(); ++im){
            // passing molecule to the helper function by value is important
            // for inheritance reasons --- RATTLE, which derives from Verlet,
            // will use both old information from time (t) as well as new
            // information at (t + dt) from unconstrained Verlet integration
            polymer = state.polymers.at(im);
            state.polymers.at(im) = _move_verlet_full_step(polymer,
                calculate_observables);
        }
        // solvent loop
        for (int is = 0; is < state.nsolvents(); ++is){
            solvent = state.solvents.at(is);
            state.solvents.at(is) = _move_verlet_full_step(solvent,
                calculate_observables);
        }
        // advance state time
        // now velocity, positions, and forces of all atoms have to be at t + dt
        state.advance_time(timestep);
        // multiply the accumulators by whatever factors necessary
}
