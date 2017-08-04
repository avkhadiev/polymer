// 2017 Artur Avkhadiev
/*! \file force_updater.h
*/
#ifndef POLYMER_FORCE_UPDATER_H
#define POLYMER_FORCE_UPDATER_H
#include "ljpotential.h"
#include "state.h"
#include "simple_state.h"
// Is only meant to be called inside the function update_forces

class ForceUpdater {
private:
    LJPotential _potential;
    // pointers to accumulators of observables
    bool _is_potential_energy_set;
    bool _is_neg_virial_acc_set;
    double *_potential_energy_acc;
    double *_neg_virial_acc;
    // Calculates force on atom i from atom j, fij. I
    // increases the atom_i->force.first and
    // decreases atom_j->force.first by fij.
    // Does not change the force time records (atom->force.second) of the atoms.
    void _update_forces_in_atomic_pair(Atom *atom_i, Atom *atom_j, bool calculate_observables);
    void _update_forces_in_atomic_pair(simple::Atom &atom_i,
        simple::Atom &atom_j,
        bool calculate_observables);
    void _update_forces_intramolecular(State& state, bool calculate_observables);
    void _update_forces_intramolecular(simple::AtomState& state, bool calculate_observables);
    void _update_forces_intermolecular(State& state, bool calculate_observables);
    void _update_forces_intermolecular(simple::AtomState& state, bool calculate_observables);
public:
    // getters
    const LJPotential& get_potential() const;
    bool is_potential_energy_acc_set() const;
    bool is_neg_virial_acc_set() const;
    // setters
    void set_potential(LJPotential& potential);
    void set_potential_energy_acc(double *potential_energy_acc);
    void set_neg_virial_acc(double *neg_virial_acc);
    // main functions
    void update_forces(State& state, bool calculate_observables = false);
    void update_forces(simple::AtomState &state,
        bool calculate_observables = false);
    ForceUpdater(LJPotential potential,
        double *potential_energy_acc = NULL,
        double *_neg_virial_acc = NULL);
    ~ForceUpdater();
};
#endif
