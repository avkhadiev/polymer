// 2017 Artur Avkhadiev
/*! \file force_updater.h
*/
#ifndef POLYMER_FORCE_UPDATER_H
#define POLYMER_FORCE_UPDATER_H
#include "ljpotential.h"
#include "simple_state.h"
#include "observable.h"
// Is only meant to be called inside the function update_forces

class ForceUpdater {
private:
    LJPotential _potential;
    // pointers to accumulators of observables
    struct obs_t {
        bool is_set;
        Observable *ptr;
    } _pe, _w;
    // Calculates force on atom i from atom j, fij. I
    // increases the atom_i->force.first and
    // decreases atom_j->force.first by fij.
    // Does not change the force time records (atom->force.second) of the atoms.
    void _update_forces_in_atomic_pair(simple::Atom &atom_i,
        simple::Atom &atom_j,
        bool calculate_observables);
    void _update_forces_intramolecular(simple::AtomState& state, bool calculate_observables);
    void _update_forces_intermolecular(simple::AtomState& state, bool calculate_observables);
public:
    // getters
    const LJPotential& get_potential() const;
    // setters
    void set_potential(LJPotential& potential);
    void set_pe(Observable *pe);
    void set_w(Observable *w);
    // main functions
    void update_forces(simple::AtomState &state,
        bool calculate_observables = false);
    ForceUpdater(LJPotential potential,
        Observable *pe = NULL,
        Observable *w = NULL);
    ~ForceUpdater();
};
#endif
