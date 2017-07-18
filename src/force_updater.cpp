// 2017 Artur Avkhadiev
/*! \file force_updater.h
*/
#include <string>
#include <stdexcept>
#include "../include/force_updater.h"
#include "../include/vector.h"
#include "../include/ljpotential.h"
#include "../include/state.h"
ForceUpdater::ForceUpdater(LJPotential potential,
    double *potential_energy_acc,
    double *neg_virial_acc) :
    _potential (potential),
    _is_potential_energy_set (false),
    _is_neg_virial_acc_set (false) {
    // if non-null pointers to accumulators are given, they will be used to
    // calculate the respective observables
    if(potential_energy_acc != NULL){
        _is_potential_energy_set = true;
        _potential_energy_acc = potential_energy_acc;
    }
    if(neg_virial_acc != NULL){
        _is_neg_virial_acc_set = true;
        _neg_virial_acc = neg_virial_acc;
    }
}
ForceUpdater::~ForceUpdater(){}
const LJPotential& ForceUpdater::get_potential() const {
    return _potential;
}
bool ForceUpdater::is_potential_energy_acc_set() const {
    return _is_potential_energy_set;
}
bool ForceUpdater::is_neg_virial_acc_set() const {
    return _is_potential_energy_set;
}
void ForceUpdater::set_potential(LJPotential& potential){
    _potential = potential;
}
void ForceUpdater::set_potential_energy_acc(double *potential_energy_acc){
    if (potential_energy_acc != NULL){
        _is_potential_energy_set = true;
    }
    else {
        _is_potential_energy_set = false;
    }
    _potential_energy_acc = potential_energy_acc;
}
void ForceUpdater::set_neg_virial_acc(double *neg_virial_acc){
    if (neg_virial_acc != NULL){
        _is_neg_virial_acc_set = true;
    }
    else {
        _is_neg_virial_acc_set = false;
    }
    _neg_virial_acc = neg_virial_acc;
}
void ForceUpdater::_update_forces_in_atomic_pair(Atom *atom_i, Atom *atom_j, bool calculate_observables){
    Vector rij = subtract(atom_i->position.first, atom_j->position.first);
    double rijsq = normsq(rij);
    // if the distances between atoms is zero
    if (rijsq == 0) {
        std::string err_msg = "update_forces_in_atomic_pair: distance between atoms is 0.";
        throw std::invalid_argument(err_msg);
    }
    else {
        double inv_rijsq = 1/rijsq;
        // calculate the force of atom j on atom i
        Vector fij = multiply(rij, _potential.calculate_fstrength_over_r(inv_rijsq));
        // update forces
        atom_i->force.first += fij;
        atom_j->force.first -= fij;
        // if observables need to be calculated and respective pointers are set,
        // zero the accumulators
        if(calculate_observables){
            if(_is_potential_energy_set){
                double vij = _potential.calculate_pair_potential(inv_rijsq);
                *_potential_energy_acc += vij;
            }
            if(_is_neg_virial_acc_set){
                double wij = _potential.calculate_neg_pair_virial(inv_rijsq);
                *_neg_virial_acc += wij;
            }
        }
    }
};
void ForceUpdater::_update_forces_intramolecular(State &state, bool calculate_observables){
    Molecule *mol;
    Atom *atom_i;
    Atom *atom_j;
    // loop over all intramolecular interactions
    for(int im = 0; im < state.nm; ++im){
        mol = &(state.molecules.at(im));
        // loop over all atomic pairs in the molecule
        for(int ia = 0; ia < mol->na - 1; ++ia){
            atom_i = &(mol->atoms.at(ia));
            for(int ja = ia + 1; ja < mol->na; ++ja){
                atom_j = &(mol->atoms.at(ja));
                _update_forces_in_atomic_pair(atom_i, atom_j, calculate_observables);
            }
        }
    }
}
void ForceUpdater::_update_forces_intermolecular(State &state, bool calculate_observables){
    Molecule *mol_i;
    Molecule *mol_j;
    Atom *atom_i;
    Atom *atom_j;
    // loop over all intermolecular interactions
    for(int im = 0; im < state.nm - 1; ++im){
        mol_i = &(state.molecules.at(im));
        for(int jm = im + 1; jm < state.nm; ++jm){
            mol_j = &(state.molecules.at(jm));
            // loop over all atoms in molecule i
            for(int ia = 0; ia < mol_i->na; ++ia){
                atom_i = &(mol_i->atoms.at(ia));
                for(int ja = 0; ja < mol_j->na; ++ja){
                    atom_j = &(mol_j->atoms.at(ja));
                    _update_forces_in_atomic_pair(atom_i, atom_j, calculate_observables);
                }
            }
        }
    }
}
void ForceUpdater::update_forces(State &state, bool calculate_observables){
        Molecule *mol;
        Atom *atom;
        // set all forces to zero
        for(int im = 0; im < state.nm; ++im){
            mol = &(state.molecules.at(im));
            for(int ia = 0; ia < mol->na; ++ia){
                atom = &(mol->atoms.at(ia));
                // update time records on forces to coincide with the position times
                atom->force.first = vector(0.0, 0.0, 0.0);
            }
        }
        // if observables need to be calculated and respective pointers are set,
        // zero the accumulators
        if(calculate_observables){
            if(_is_potential_energy_set){
                *_potential_energy_acc = 0.0;
            }
            if(_is_neg_virial_acc_set){
                *_neg_virial_acc = 0.0;
            }
        }
        // update all inter- and intra- molecular interactions
        _update_forces_intermolecular(state, calculate_observables);
        _update_forces_intramolecular(state, calculate_observables);
        // divide the negative of the sum of pair virials by 3
        // (see definition of virial in Allen & Tildesley)
        if(calculate_observables){
            if(_is_neg_virial_acc_set){
                *_neg_virial_acc = *_neg_virial_acc / 3.0;
            }
        }
        // forces are updated; advance times attached to all forces
        for(int im = 0; im < state.nm; ++im){
            mol = &(state.molecules.at(im));
            for(int ia = 0; ia < mol->na; ++ia){
                atom = &(mol->atoms.at(ia));
                // update time records on forces to coincide with the position times
                atom->force.second = atom->position.second;
            }
        }
}
