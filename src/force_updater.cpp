// 2017 Artur Avkhadiev
/*! \file force_updater.h
*/
#include <string>
#include <stdexcept>
#include "../include/force_updater.h"
#include "../include/vector.h"
#include "../include/ljpotential.h"
#include "../include/simple_state.h"
ForceUpdater::ForceUpdater(Potential* polymer_potential,
    Potential *solvent_potential,
    Potential *inter_potential,
    Observable *pe,
    Observable *w) :
    _polymer_potential (polymer_potential),
    _solvent_potential (solvent_potential),
    _inter_potential (inter_potential){
        _pe.ptr = pe;
        _w.ptr = w;
        if ((polymer_potential == NULL)
            || (solvent_potential == NULL)
            || (inter_potential == NULL)){
            std::string err_msg = "force updater initialization: NULL pointer to one of the potentials given!";
            throw std::invalid_argument(err_msg);
        }
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
    simple::Atom &atom_j, Potential* potential, bool calculate_observables){
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
            Vector fij = multiply(rij, _polymer_potential->calculate_fstrength_over_r(inv_rijsq));
            // update forces
            atom_i.force += fij;
            atom_j.force -= fij;
            // if observables need to be calculated and respective pointers are set,
            // zero the accumulators
            if(calculate_observables){
                if(_pe.is_set){
                    double vij = _polymer_potential->calculate_pair_potential(inv_rijsq);
                    _pe.ptr->update( _pe.ptr->value() + vij );
                }
                if(_w.is_set){
                    double wij = _polymer_potential->calculate_neg_pair_virial(inv_rijsq);
                    _w.ptr->update( _w.ptr->value() + wij );
                }
            }
        }
}
void ForceUpdater::_update_forces_intramolecular(simple::AtomState &state, bool calculate_observables){
    // loop over all intramolecular interactions
    int na = state.polymer_nb() + 1;
    for(int im = 0; im < state.nm(); ++im){
        // loop over all atomic pairs in polymers,
        // except the ones connected with a bond
        for(int ia = 0; ia < na - 2; ++ia){
            simple::Atom& atom_i = state.polymers.at(im).atoms.at(ia);
            for(int ja = ia + 2; ja < na; ++ja){
                simple::Atom& atom_j = state.polymers.at(im).atoms.at(ja);
                _update_forces_in_atomic_pair(atom_i, atom_j, _polymer_potential, calculate_observables);
            }
        }
    }
    // loop over all solvent molecules
    for(int is = 0; is < state.nsolvents() - 1; ++is){
        simple::Atom& atom_i = state.solvents.at(is).atom;
        for(int js = 0; js < state.nsolvents(); ++js){
            simple::Atom& atom_j = state.solvents.at(js).atom;
            _update_forces_in_atomic_pair(atom_i, atom_j, _solvent_potential, calculate_observables);
        }
    }
}
void ForceUpdater::_update_forces_intermolecular(simple::AtomState &state,
    bool calculate_observables){
    int na = state.polymer_nb() + 1;
    // loop over all intermolecular interactions
    for(int im = 0; im < state.nm(); ++im){
        // polymer-polymer
        for(int jm = im + 1; jm < state.nm(); ++jm){
            for(int ia = 0; ia < na; ++ia){
                simple::Atom& atom_i = state.polymers.at(im).atoms.at(ia);
                for(int ja = 0; ja < na; ++ja){
                    simple::Atom& atom_j = state.polymers.at(im).atoms.at(ja);
                    _update_forces_in_atomic_pair(atom_i, atom_j, _polymer_potential, calculate_observables);
                }
            }
        }
        // polymer-solvent
        for(int ia = 0; ia < na; ++ia){
            simple::Atom& atom_i = state.polymers.at(im).atoms.at(ia);
            for(int js = 0; js < state.nsolvents(); ++js){
                simple::Atom& atom_j = state.solvents.at(js).atom;
                _update_forces_in_atomic_pair(atom_i, atom_j, _inter_potential, calculate_observables);
            }
        }
    }
}
void ForceUpdater::update_forces(simple::AtomState &state, bool calculate_observables){
        int na = state.polymer_nb() + 1;
        // set all forces to zero
        // for polymers
        for(int im = 0; im < state.nm(); ++im){
            for(int ia = 0; ia < na; ++ia){
                simple::Atom& atom = state.polymers.at(im).atoms.at(ia);
                atom.force = vector(0.0, 0.0, 0.0);
            }
        }
        // for solvents
        for (int is = 0; is < state.nsolvents(); ++is){
            state.solvents.at(is).set_f(vector(0.0, 0.0, 0.0));
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
}
