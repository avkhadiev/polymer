// 2017 Artur Avkhadiev
/*! \file force_updater.h
*/
#ifndef POLYMER_FORCE_UPDATER_H
#define POLYMER_FORCE_UPDATER_H
#include "potential.h"
#include "simple_state.h"
#include "observable.h"
// Is only meant to be called inside the function update_forces
class ForceUpdater {
private:
    Potential* _polymer_potential;
    Potential* _solvent_potential;
    Potential* _inter_potential;
    // pointers to accumulators of observables
    struct obs_t {
        bool is_set;
        Observable *ptr;
    } _pe, _w;
    // Calculates force on atom i from atom j, fij. I
    // increases the atom_i->force.first and
    // decreases atom_j->force.first by fij.
    // if use_cutoff, attempts to get cutoff information from the potential
    void _update_forces_in_atomic_pair(simple::Atom &atom_i,
        simple::Atom &atom_j,
        Potential *potential,
        bool calculate_observables);
    void _update_forces_intramolecular(simple::AtomState& state, bool calculate_observables);
    void _update_forces_intermolecular(simple::AtomState& state, bool calculate_observables);
public:
    // getters
    const Potential* get_polymer_potential() const {return _polymer_potential;};
    const Potential* get_solvent_potential() const {return _solvent_potential;};
    const Potential* get_inter_potential() const {return _inter_potential;};
    // setters
    void set_polymer_potential(Potential *p) {_polymer_potential = p;};
    void set_solvent_potential(Potential *p) {_solvent_potential = p;};
    void set_inter_potential(Potential *p) {_inter_potential = p;};
    void set_pe(Observable *pe);
    void set_w(Observable *w);
    void zero_observables();
    // main functions
    double calc_pot_energy(const simple::AtomState &state);
    void update_forces(simple::AtomState &state,
        bool calculate_observables = false);
    ForceUpdater(Potential* polymer_potential,
        Potential* solvent_potential,
        Potential* inter_potential,
        Observable *pe = NULL,
        Observable *w = NULL);
    ~ForceUpdater();
};
#endif
