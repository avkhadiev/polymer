// 2017 Artur Avkhadiev
/*! \file rattle_integrator.cpp
*/
#include <vector>
#include <string>
#include <algorithm>                                /* std::fill */
#include <cmath>                                    /* pow */
#include <cassert>
#include "../include/rattle_integrator.h"
RattleIntegrator::RattleIntegrator(ForceUpdater force_updater,
    double tol,
    double tiny,
    int maxiter,
    Observable *ke,
    Observable *wc) :
    VerletIntegrator(force_updater, NULL),
    _maxiter (maxiter),
    _tol (tol),
    _tiny (tiny),
    _tol2 (2 * _tol),
    _nb (simple::BasePolymer::nb()),
    _dabsq (pow(simple::BasePolymer::d(), 2.0)),
    _rm (1.0/(simple::BasePolymer::m())),
    _inv_timestep (0.0)         /**> will be set in move() by _set_timestep()*/
{
    _moving.resize(_nb);
    _moved.resize(_nb);
    _ke.ptr = ke;
    _wc.ptr = wc;
    if(_ke.ptr == NULL){
        _ke.is_set = false;
    }
    else{
        _ke.is_set = true;
    }
    if(_wc.ptr == NULL){
        _wc.is_set = false;
    }
    else{
        _wc.is_set = true;
    }
}
RattleIntegrator::~RattleIntegrator(){};
double RattleIntegrator::get_tol() const {
    return _tol;
}
double RattleIntegrator::get_tiny() const {
    return _tiny;
}
void RattleIntegrator::_set_timestep(double timestep){
    _timestep = timestep;
    _halfstep = timestep * 0.5;
    _inv_timestep = 1 / _timestep;
}
void RattleIntegrator::set_tol(double tol){
    _tol = tol;
    _tol2 = _tol * 2.0;
}
void RattleIntegrator::set_tiny(double tiny){
    _tiny = tiny;
}
void RattleIntegrator::set_ke(Observable *ptr){
    if (ptr != NULL){
        _ke.is_set = true;
    }
    else {
        _ke.is_set = false;
    }
    _ke.ptr = ptr;
}
void RattleIntegrator::set_wc(Observable *ptr){
    if (ptr != NULL){
        _wc.is_set = true;
    }
    else {
        _wc.is_set = false;
    }
    _wc.ptr = ptr;
}
bool RattleIntegrator::_is_constraint_within_tol(double dabsq,
    double difference_of_squares){
    // a factor of two in the tolerance comes from a Taylor expansion:
    // instead of considering the difference of values, we are considering
    // the difference of squares of values.
    bool is_constraint_within_tol
        = (std::abs(difference_of_squares) < _tol2 * dabsq);
    return is_constraint_within_tol;
}
bool RattleIntegrator::_is_angle_okay(double dabsq, double rr_dot){
    // will be dividing by this dot product to find constraint forces
    // physically, given a reasonable timestep, the angle between a bond vector
    // at time t and the same bond vector at time t+dt should be small =>
    // cosine of the angle is close to 1
    bool is_angle_okay_tol = (rr_dot > _tiny * dabsq);
    return is_angle_okay_tol;
}
bool RattleIntegrator::_is_constraint_derivative_within_rvtol(double dabsq,
    double rv_dot){
    // the projection of bond velocity onto the bond vector should not
    // increase the bond length by more than the bondlength tolerance within
    // the given timestep
    bool is_constraint_derivative_within_rvtol =
        (std::abs(rv_dot) < _tol * dabsq);
    return is_constraint_derivative_within_rvtol;
}
void RattleIntegrator::_zero_accumulators(){
    if(_ke.is_set){
        _ke.ptr->zero();
    }
    if(_wc.is_set){
        _wc.ptr->zero();
    }
}
void RattleIntegrator::_correct_accumulators(){
    if(_wc.is_set){
        // ftp://ftp.dl.ac.uk/ccp5/ALLEN_TILDESLEY/F.09
        // here, a *negative* of the constraint virial is calculated
        _wc.ptr->update(_wc.ptr->value() * (-2.0) * _inv_timestep);
    }
}
void RattleIntegrator::_set_up_correction_bookkeeping(){
    // all constrained bonds in the molecule have been moved
    std::fill(_moved.begin(), _moved.end(), true);
    // none of the bounds are moving
    std::fill(_moving.begin(), _moving.end(), false);
}
simple::AtomPolymer RattleIntegrator::_move_correct_half_step(
    simple::AtomPolymer molecule_old,
    simple::AtomPolymer molecule_new){
    // set up "moved" and "moving" arrays
    _set_up_correction_bookkeeping();
    bool still_correcting = true;
    int iter = 0;
    double rr_dot;              /**> r_{AB}(t) * r^{i}_{AB}(t+dt) */
    double dabsq = _dabsq;      /**> fixed bond length squared */
    double diffsq;              /**> d^2_{AB} - r^{i}_{AB}(t+dt) */
    double gab;                 /**> proportional to constraint force */
    double rma = _rm;           /**> 1/mass of atom 1 */
    double rmb = _rm;           /**> 1/mass of atom 2 */
    double dt = _timestep;
    int nb = simple::BasePolymer::nb();
    while (still_correcting && (iter < _maxiter)){
        /***********************************************************/
        iter += 1;
        still_correcting = false;
        /***********************************************************/
        for(int ib = 0; ib < nb; ++ib){
            if (_moved.at(ib) == false){
                continue;
            }
            else{
                Vector& r_old_atom1 = molecule_old.atoms.at(ib).position;
                Vector& r_old_atom2 = molecule_old.atoms.at(ib+1).position;
                Vector& r_new_atom1 = molecule_new.atoms.at(ib).position;
                Vector& r_new_atom2 = molecule_new.atoms.at(ib+1).position;
                Vector& v_new_atom1 = molecule_new.atoms.at(ib).velocity;
                Vector& v_new_atom2 = molecule_new.atoms.at(ib+1).velocity;
                Vector rab_old = subtract(r_old_atom1, r_old_atom2);
                Vector rab_new = subtract(r_new_atom1, r_new_atom2);
                // Order of terms is important! Determines the sign of the
                // constraint force
                diffsq = dabsq - normsq(rab_new);
                if(_is_constraint_within_tol(dabsq, diffsq) == false){
                    // NEED TO CORRECT POSITIONS AND VELOCITIES
                    /**********************************************************/
                    still_correcting = true;
                    _moving.at(ib) = true;
                    /**********************************************************/
                    // check angle between r_AB(t) and r_AB(t+dt)
                    rr_dot = dot(rab_old, rab_new);
                    if(_is_angle_okay(dabsq, rr_dot) == false){
                        std::string err_msg = "move_correct_half_step: r_{AB}(t) and r_{AB}(t + dt) are nearly normal to each other!\n";
                        err_msg += "r_{AB}(t) = " + vector_to_string(rab_old);
                        err_msg += "\nr_{AB(t + dt) = " + vector_to_string(rab_new);
                        throw std::invalid_argument(err_msg);
                    }
                    // COMMENCE CORRECTION
                    // calculate constraint force
                    gab = diffsq/(2 * dt * (rma + rmb) * rr_dot);
                    // correct the atomic velocities
                    Vector dva = multiply(rab_old, rma * gab);
                    Vector dvb = multiply(-rab_old, rmb * gab);
                    v_new_atom1 += dva;
                    v_new_atom2 += dvb;
                    // correct atomic position
                    Vector dra = multiply(dva, dt);
                    Vector drb = multiply(dvb, dt);
                    r_new_atom1 += dra;
                    r_new_atom2 += drb;
                }
            }
        }
    // what was moving is moved, nothing is moving
    _moved = _moving;
    std::fill(_moving.begin(), _moving.end(), false);
    }
    if (iter >= _maxiter){
        fprintf(stderr, "%s (%d)\n",
            "RATTLE MOVE A: maximum number of iterations reached. Stopping correction",
            _maxiter);
    }
    return molecule_new;
}
simple::AtomPolymer RattleIntegrator::_move_correct_full_step(
    simple::AtomPolymer molecule,
    bool calculate_observables){
    // set up "moved" and "moving" arrays
    _set_up_correction_bookkeeping();
    bool still_correcting = true;
    int iter = 0;
    double rv_dot;              /**> r_{AB}(t) * r^{i}_{AB}(t+dt) */
    double dabsq = _dabsq;      /**> fixed bond length squared */
    double kab;                 /**> proportional to constraint force */
    double rma = _rm;           /**> 1/mass of atom 1 */
    double rmb = _rm;           /**> 1/mass of atom 2 */
    int nb = simple::BasePolymer::nb();
    while (still_correcting && (iter < _maxiter)){
        /***********************************************************/
        iter += 1;
        still_correcting = false;
        /***********************************************************/
        for(int ib = 0; ib < nb; ++ib){
            if (_moved.at(ib) == false){
                continue;
            }
            else{
                Vector& r_atom1 = molecule.atoms.at(ib).position;
                Vector& r_atom2 = molecule.atoms.at(ib + 1).position;
                Vector& v_atom1 = molecule.atoms.at(ib).velocity;
                Vector& v_atom2 = molecule.atoms.at(ib + 1).velocity;
                Vector rab = subtract(r_atom1, r_atom2);
                Vector vab = subtract(v_atom1, v_atom2);
                rv_dot = dot(rab, vab);
                if(_is_constraint_derivative_within_rvtol(dabsq, rv_dot) == false){
                    // NEED TO CORRECT VELOCITIES
                    /**********************************************************/
                    still_correcting = true;
                    _moving.at(ib) = true;
                    /**********************************************************/
                    // COMMENCE CORRECTION
                    // calculate constraint force
                    kab = - rv_dot /((rma + rmb) * dabsq);
                    // correct the atomic velocities
                    Vector dva = multiply(rab, rma * kab);
                    Vector dvb = multiply(-rab, rmb * kab);
                    v_atom1 += dva;
                    v_atom2 += dvb;
                    // increment negative constraint virial accumulator
                    // if neccessary.
                    if(calculate_observables && _wc.is_set){
                        _wc.ptr->update(_wc.ptr->value() + kab * dabsq);
                    }
                }
            }
        }
        // what was moving is moved, nothing is moving
        _moved = _moving;
        std::fill(_moving.begin(), _moving.end(), false);
    }
    // all the corrections have been performed; calculate kinetic energy
    // update kinetic energy accumulator if necessary and possible
    if(calculate_observables && _ke.is_set) _ke.ptr->update(molecule);
    if (iter >= _maxiter){
        fprintf(stderr, "%s (%d)\n",
            "RATTLE MOVE B: maximum number of iterations reached. Stopping correction.",
            _maxiter);
    }
    return molecule;
}
void RattleIntegrator::move(double timestep,
    simple::AtomState &state,
    bool calculate_observables){
    if(calculate_observables){
        _zero_accumulators();
        _force_updater.zero_observables();
    }
    _set_timestep(timestep);
    simple::AtomPolymer molecule_last_step;
    simple::AtomPolymer molecule_half_step_to_correct;
    // loop over molecules, perform verlet half-step and correction
    for (int im = 0; im < state.nm(); ++im){
        // molecule contains positions and velocities at last step
        molecule_last_step = state.polymers.at(im);
        // pass molecule at last step to verlet integrator by value
        // the integrator returns the molecule
        // with uncorrected full-step positions
        // and uncorrected half-step velocities
        molecule_half_step_to_correct = _move_verlet_half_step(molecule_last_step);
        // pass both last-step molecule and molecule from verlet integrator
        // to RATTLE for iterative correction
        // molecule with corrected full-step positions and
        // corrected half-step velocities is stored in the state
        state.polymers.at(im) = _move_correct_half_step(molecule_last_step, molecule_half_step_to_correct);
    }
    // r(t + dt) has been calculated and corrected
    // now update forces in all atoms from f(t) to f(t + dt)
    // if calculate_observables = true, accumulators of observables that
    // should be updated in the force loop will be updated in this iteration
    _force_updater.update_forces(state, calculate_observables);
    // loop over molecules, perfrom verlet full-step and correction
    for (int im = 0; im < state.nm(); ++im){
        // molecule contains corrected positions at full step,
        // corrected velocities at half step, and forces calculated at
        // full step
        simple::AtomPolymer molecule = state.polymers.at(im);
        // pass molecule to verlet integrator, which returns the molecule
        // with uncorrected velocities at full step.
        // store returned molecule in memory
        molecule = _move_verlet_full_step(molecule, calculate_observables);
        // pass the molecule to RATTLE to iteratively correct velocities
        // at full step. Store the returned molecule in memory.
        state.polymers.at(im) = _move_correct_full_step(molecule, calculate_observables);
    }
    // advance state time
    // now velocity, positions, and forces of all atoms have to be at t + dt
    state.advance_time(timestep);
    // multiply accumulators by whatever factors necessary
    if(calculate_observables) _correct_accumulators();
}
