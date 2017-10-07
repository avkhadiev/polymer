// 2017 Artur Avkhadiev
/*! \file force_updater.h
*/
#include <string>
#include <stdexcept>
#include "../include/force_updater.h"
#include "../include/vector.h"
#include "../include/ljpotential.h"
#include "../include/simple_state.h"
ForceUpdater::ForceUpdater(LJPotential potential,
    Observable *pe,
    Observable *w) :
    _potential (potential){
    _pe.ptr = pe;
    _w.ptr = w;
    if(_pe.ptr == NULL){
        _pe.is_set = false;
    }
    else{
        _pe.is_set = true;
    }
    if(_w.ptr == NULL){
        _w.is_set = false;
    }
    else{
        _w.is_set = true;
    }
}
ForceUpdater::~ForceUpdater(){}
const LJPotential& ForceUpdater::get_potential() const {
    return _potential;
}
void ForceUpdater::set_potential(LJPotential& potential){
    _potential = potential;
}
void ForceUpdater::set_pe(Observable *ptr){
    if (ptr != NULL){
        _pe.is_set = true;
    }
    else {
        _pe.is_set = false;
    }
    _pe.ptr = ptr;
}
void ForceUpdater::set_w(Observable *ptr){
    if (ptr != NULL){
        _w.is_set = true;
    }
    else {
        _w.is_set = false;
    }
    _w.ptr = ptr;
}
void ForceUpdater::zero_observables(){
    if (_pe.is_set) _pe.ptr->zero();
    if (_w.is_set) _w.ptr->zero();
}
void ForceUpdater::_update_forces_in_atomic_pair(simple::Atom &atom_i,
    simple::Atom &atom_j,
    bool calculate_observables){
        Vector rij = subtract(atom_i.position, atom_j.position);
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
            atom_i.force += fij;
            atom_j.force -= fij;
            // if observables need to be calculated and respective pointers are set,
            // zero the accumulators
            if(calculate_observables){
                if(_pe.is_set){
                    double vij = _potential.calculate_pair_potential(inv_rijsq);
                    _pe.ptr->update( _pe.ptr->value() + vij );
                }
                if(_w.is_set){
                    double wij = _potential.calculate_neg_pair_virial(inv_rijsq);
                    _w.ptr->update( _w.ptr->value() + wij );
                }
            }
        }
}
void ForceUpdater::_update_forces_intramolecular(simple::AtomState &state, bool calculate_observables){
    int na = state.polymer_nb() + 1;
    // loop over all intramolecular interactions
    for(int im = 0; im < state.nm(); ++im){
        // loop over all atomic pairs, except the ones connected with a bond
        for(int ia = 0; ia < na - 2; ++ia){
            simple::Atom& atom_i = state.polymers.at(im).atoms.at(ia);
            for(int ja = ia + 2; ja < na; ++ja){
                simple::Atom& atom_j = state.polymers.at(im).atoms.at(ja);
                _update_forces_in_atomic_pair(atom_i, atom_j,
                    calculate_observables);
            }
        }
    }
}
void ForceUpdater::_update_forces_intermolecular(simple::AtomState &state,
    bool calculate_observables){
    int na = state.polymer_nb() + 1;
    // loop over all intermolecular interactions
    for(int im = 0; im < state.nm(); ++im){
        for(int jm = im + 1; jm < state.nm(); ++jm){
            // loop over all atoms in molecule i
            for(int ia = 0; ia < na; ++ia){
                simple::Atom& atom_i = state.polymers.at(im).atoms.at(ia);
                for(int ja = 0; ja < na; ++ja){
                    simple::Atom& atom_j = state.polymers.at(im).atoms.at(ja);
                    _update_forces_in_atomic_pair(atom_i, atom_j, calculate_observables);
                }
            }
        }
    }
}
void ForceUpdater::update_forces(simple::AtomState &state, bool calculate_observables){
        int na = state.polymer_nb() + 1;
        // set all forces to zero
        for(int im = 0; im < state.nm(); ++im){
            for(int ia = 0; ia < na; ++ia){
                simple::Atom& atom = state.polymers.at(im).atoms.at(ia);
                atom.force = vector(0.0, 0.0, 0.0);
            }
        }
        // if observables need to be calculated and respective pointers are set,
        // zero the accumulators
        if(calculate_observables){
            if(_pe.is_set) _pe.ptr->zero();
            if(_w.is_set) _w.ptr->zero();
        }
        // update all inter- and intra- molecular interactions
        _update_forces_intermolecular(state, calculate_observables);
        _update_forces_intramolecular(state, calculate_observables);
        // divide the negative of the sum of pair virials by 3
        // (see definition of  in Allen & Tildesley)
        if(calculate_observables && _w.is_set) {
            _w.ptr->update( _w.ptr->value() / 3.0 );
        }
}
