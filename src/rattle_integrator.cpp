// 2017 Artur Avkhadiev
/*! \file rattle_integrator.cpp
*/
#include <vector>
#include <string>
#include <algorithm>                                /* std::fill */
#include <cmath>                                    /* pow */
#include <cassert>
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
    _inv_timestep (0.0),         /**> will be set in move() by _set_timestep()*/
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
void RattleIntegrator::_set_timestep(double timestep){
    _timestep = timestep;
    _halfstep = timestep * 0.5;
    _inv_timestep = 1 / _timestep;
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
bool RattleIntegrator::_is_constraint_within_tol(double dabsq,
    double difference_of_squares){
    bool is_constraint_within_tol
        = (std::abs(difference_of_squares) < _tol2 * dabsq);
    return is_constraint_within_tol;
}
bool RattleIntegrator::_is_angle_okay(double dabsq, double rr_dot){
    // will be dividing by this dot produduct to find constraint forces
    bool is_angle_okay_tol = (rr_dot > _tiny * dabsq);
    return is_angle_okay_tol;
}
bool RattleIntegrator::_is_constraint_derivative_within_rvtol(double dabsq,
    double rv_dot){
    bool is_constraint_derivative_within_rvtol = (rv_dot < _rvtol2 * dabsq);
    return is_constraint_derivative_within_rvtol;
}
Molecule RattleIntegrator::_move_correct_half_step(Molecule molecule_last_step,
    Molecule molecule_half_step_to_correct){
    // notational shorthand
    Molecule& molecule_old = molecule_last_step;
    Molecule& molecule_new = molecule_half_step_to_correct;
    // set up "moved" and "moving" arrays
    _set_up_correction_bookkeeping(molecule_new);
    bool still_correcting = true;
    int iter = 0;
    double rr_dot;              /**> r_{AB}(t) * r^{i}_{AB}(t+dt) */
    double dabsq;               /**> fixed bond length squared */
    double diffsq;              /**> d^2_{AB} - r^{i}_{AB}(t+dt) */
    double gab;                 /**> proportional to constraint force */
    double rma;                 /**> 1/mass of atom 1 */
    double rmb;                 /**> 1/mass of atom 2 */
    double dt = _timestep;
    double rdt = _inv_timestep;
    Atom *atom_old1;
    Atom *atom_old2;
    Atom *atom_new1;
    Atom *atom_new2;
    Bond *bond_old;
    Bond *bond_new;
    while (still_correcting && (iter < _maxiter)){
        /***********************************************************/
        iter += 1;
        still_correcting = false;
        /***********************************************************/
        assert(molecule_old.nb == molecule_new.nb);
        int nb = molecule_new.nb;
        for(int ib = 0; ib < nb; ++ib){
            if (_moved.at(ib) == false){
                continue;
            }
            else{
                bond_old = &(molecule_old.bonds.at(ib));
                bond_new = &(molecule_new.bonds.at(ib));
                assert(bond_old->fixed_length_sq == bond_new->fixed_length_sq);
                assert(bond_old->atom1 == bond_new->atom1);
                assert(bond_old->atom2 == bond_new->atom2);
                dabsq = bond_new->fixed_length_sq;
                atom_old1 = &(molecule_old.atoms.at(bond_old->atom1));
                atom_old2 = &(molecule_old.atoms.at(bond_old->atom2));
                atom_new1 = &(molecule_new.atoms.at(bond_new->atom1));
                atom_new2 = &(molecule_new.atoms.at(bond_new->atom2));
                Vector *r_old_atom1 = &(atom_old1->position.first);
                Vector *r_old_atom2 = &(atom_old2->position.first);
                Vector *r_new_atom1 = &(atom_new1->position.first);
                Vector *r_new_atom2 = &(atom_new2->position.first);
                Vector *v_new_atom1 = &(atom_new1->velocity.first);
                Vector *v_new_atom2 = &(atom_new2->velocity.first);
                rma = 1 / (atom_new1->mass);
                rmb = 1 / (atom_new2->mass);
                Vector rab_old = subtract(*r_old_atom1, *r_old_atom2);
                Vector rab_new = subtract(*r_new_atom1, *r_new_atom2);
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
                        std::string err_msg = "move_correct_half_step: r_{AB}(t) and r_{AB}(t + dt) are nearly normal to each other!";
                        throw std::invalid_argument(err_msg);
                    }
                    // COMMENCE CORRECTION
                    // calculate constraint force
                    gab = diffsq/(2 * rdt * (rma + rmb) * rr_dot);
                    // correct the atomic velocities
                    Vector dva = multiply(rab_new, rma * gab);
                    Vector dvb = multiply(-rab_new, rmb * gab);
                    *v_new_atom1 += dva;
                    *v_new_atom2 += dvb;
                    // correct atomic position
                    Vector dra = multiply(dva, dt);
                    Vector drb = multiply(dvb, dt);
                    *r_new_atom1 += dra;
                    *r_new_atom2 += drb;
                }
            }
        }
    // what was moving is moved, nothing is moving
    _moved = _moving;
    std::fill(_moving.begin(), _moving.end(), false);
    }
    return molecule_half_step_to_correct;
}
Molecule RattleIntegrator::_move_correct_full_step(Molecule molecule_full_step_to_correct,
    bool calculate_observables){
    // notational shorthand
    Molecule& molecule = molecule_full_step_to_correct;
    // set up "moved" and "moving" arrays
    _set_up_correction_bookkeeping(molecule);
    bool still_correcting = true;
    int iter = 0;
    double rv_dot;              /**> r_{AB}(t) * r^{i}_{AB}(t+dt) */
    double dabsq;               /**> fixed bond length squared */
    double kab;                 /**> proportional to constraint force */
    double rma;                 /**> 1/mass of atom 1 */
    double rmb;                 /**> 1/mass of atom 2 */
    Atom *atom1;
    Atom *atom2;
    Bond *bond;
    while (still_correcting && (iter < _maxiter)){
        /***********************************************************/
        iter += 1;
        still_correcting = false;
        /***********************************************************/
        int nb = molecule.nb;
        for(int ib = 0; ib < nb; ++ib){
            if (_moved.at(ib) == false){
                continue;
            }
            else{
                bond = &(molecule.bonds.at(ib));
                dabsq = bond->fixed_length_sq;
                atom1 = &(molecule.atoms.at(bond->atom1));
                atom2 = &(molecule.atoms.at(bond->atom2));
                Vector *r_atom1 = &(atom1->position.first);
                Vector *r_atom2 = &(atom2->position.first);
                Vector *v_atom1 = &(atom1->velocity.first);
                Vector *v_atom2 = &(atom2->velocity.first);
                rma = 1 / (atom1->mass);
                rmb = 1 / (atom2->mass);
                Vector rab = subtract(*r_atom1, *r_atom2);
                Vector vab = subtract(*v_atom1, *v_atom2);
                rv_dot = dot(rab, vab);
                if(_is_constraint_derivative_within_rvtol(dabsq, rv_dot) == false){
                    // NEED TO CORRECT VELOCITIES
                    /**********************************************************/
                    still_correcting = true;
                    _moving.at(ib) = true;
                    /**********************************************************/
                    // COMMENCE CORRECTION
                    // calculate constraint force
                    kab = - rv_dot /(2 * (rma + rmb) * dabsq);
                    // correct the atomic velocities
                    Vector dva = multiply(rab, rma * kab);
                    Vector dvb = multiply(-rab, rmb * kab);
                    *v_atom1 += dva;
                    *v_atom2 += dvb;
                    // increment negative constraint virial accumulator
                    // if neccessary.
                    if(calculate_observables){
                        if(_is_neg_constraint_virial_acc_set){
                            *_neg_constraint_virial_acc += kab * dabsq;
                        }
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
    if(calculate_observables){
        if(_is_kinetic_energy_acc_set){
            // loop over atoms in the molecule
            Atom *atom;
            for(int ia = 0; ia < molecule.na; ++ia){
                atom = &(molecule.atoms.at(ia));
                double k_increase = atom->mass*normsq(atom->velocity.first);
                *_kinetic_energy_acc += k_increase;
            }
        }
    }
    return molecule_full_step_to_correct;
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
            = *_neg_constraint_virial_acc * 2.0 * _inv_timestep / 3.0;
    }
}
void RattleIntegrator::_set_up_correction_bookkeeping(Molecule& molecule){
    // all constrained bonds in the molecule have been moved
    _moved.resize(molecule.nb);
    std::fill(_moved.begin(), _moved.end(), true);
    // none of the bounds are moving
    _moving.resize(molecule.nb);
    std::fill(_moving.begin(), _moving.end(), false);
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
        Molecule *molecule;
        Molecule molecule_last_step;
        Molecule molecule_half_step_to_correct;
        // loop over molecules, perform verlet half-step and correction
        for (int im = 0; im < state.nm; ++im){
            // get pointer to position of molecule in memory
            // molecule contains positions and velocities at last step
            molecule = &(state.molecules.at(im));
            molecule_last_step = *molecule;
            // pass molecule at last step to verlet integrator by value
            // the integrator returns the molecule
            // with uncorrected full-step positions
            // and uncorrected half-step velocities
            molecule_half_step_to_correct
                = _move_verlet_half_step(molecule_last_step);
            // pass both last-step molecule and molecule from verlet integrator
            // to RATTLE for iterative correction
            // molecule with corrected full-step positions and
            // corrected half-step velocities is stored in the state
            *molecule = _move_correct_half_step(molecule_last_step,
                molecule_half_step_to_correct);
        }
        // r(t + dt) has been calculated and corrected
        // now update forces in all atoms from f(t) to f(t + dt)
        // if calculate_observables = true, accumulators of observables that
        // should be updated in the force loop will be updated in this iteration
        _force_updater.update_forces(state, calculate_observables);
        // loop over molecules, perfrom verlet full-step and correction
        for (int im = 0; im < state.nm; ++im){
            // get pointer to position of molecule in memory
            // molecule contains corrected positions at full step,
            // corrected velocities at half step, and forces calculated at
            // full step
            molecule = &(state.molecules.at(im));
            // pass molecule to verlet integrator, which returns the molecule
            // with uncorrected velocities at full step.
            // store returned molecule in memory
            *molecule = _move_verlet_full_step(*molecule, false);
            // pass the molecule to RATTLE to iteratively correct velocities
            // at full step. Store the returned molecule in memory.
            *molecule = _move_correct_full_step(*molecule, calculate_observables);
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
