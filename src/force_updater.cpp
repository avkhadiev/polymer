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
    Potential *inter_potential) :
    _polymer_potential (polymer_potential),
    _solvent_potential (solvent_potential),
    _inter_potential (inter_potential){
        if ((polymer_potential == NULL)
            || (solvent_potential == NULL)
            || (inter_potential == NULL)){
            std::string err_msg = "force updater initialization: NULL pointer to one of the potentials given!";
            throw std::invalid_argument(err_msg);
        }
        if (DEBUG) fprintf(stdout, "%s\n", "ForceUpdater set up successfully!");
}
ForceUpdater::~ForceUpdater(){}
void ForceUpdater::zero_observables(){
    _polymer_potential->zero_observables();
    _solvent_potential->zero_observables();
    _inter_potential->zero_observables();
}
void ForceUpdater::_update_forces_in_atomic_pair(simple::Atom &atom_i,
    simple::Atom &atom_j, Potential* potential, bool calculate_observables){
        Vector ri = atom_i.position;
        Vector rj = atom_j.position;
        Vector fij = potential->fij(ri, rj, calculate_observables);
        // update forces
        atom_i.force += fij;
        atom_j.force -= fij;
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
}
void ForceUpdater::_update_forces_intermolecular(simple::AtomState &state,
    bool calculate_observables){
    // loop over all intermolecular interactions
    // loop over all solvent molecules
    for(int is = 0; is < state.nsolvents() - 1; ++is){
        simple::Atom& atom_i = state.solvents.at(is).atom;
        for(int js = is + 1; js < state.nsolvents(); ++js){
            simple::Atom& atom_j = state.solvents.at(js).atom;
            try
            {
                _update_forces_in_atomic_pair(atom_i, atom_j, _solvent_potential, calculate_observables);
            }
            catch (std::invalid_argument)
            {
                state.write_to_file("/Users/Arthur/stratt/polymer/test/solvents/", "log", true, true);
                fprintf(stderr, "%s\n", "wrote out crashing configuration to log");
                throw;
            }
        }
    }
    int na = state.polymer_nb() + 1;
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
void ForceUpdater::_zero_forces(simple::AtomState &state) {
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
        state.solvents.at(is).atom.force = vector(0.0, 0.0, 0.0);
    }
}
void ForceUpdater::update_forces(simple::AtomState &state, bool calculate_observables){
        _zero_forces(state);
        // if observables need to be calculated, zero the accumulators
        if(calculate_observables) zero_observables();
        // update all inter- and intra- molecular interactions
        _update_forces_intermolecular(state, calculate_observables);
        _update_forces_intramolecular(state, calculate_observables);
}
double ForceUpdater::polymer_potential_energy(const simple::AtomPolymer& polymer){
    double V = 0.0;
    int na = polymer.nb() + 1;
    Vector ri, rj;
    for(int ia = 0; ia < na - 2; ++ia){
        ri = polymer.atoms.at(ia).position;
        for(int ja = ia + 2; ja < na; ++ja){
            rj = polymer.atoms.at(ja).position;
            V += _polymer_potential->vij(ri, rj);
        }
    }
    return V;
}
double ForceUpdater::calc_pot_energy(const simple::AtomState &state){
    double V = 0.0;
    Vector ri, rj;
    // INTRAMOLECULAR
    // within a polymer
    // loop over all intramolecular interactions
    int na = state.polymer_nb() + 1;
    for(int im = 0; im < state.nm(); ++im){
        // loop over all atomic pairs in polymers,
        // except the ones connected with a bond
        for(int ia = 0; ia < na - 2; ++ia){
            ri = state.polymers.at(im).atoms.at(ia).position;
            for(int ja = ia + 2; ja < na; ++ja){
                rj = state.polymers.at(im).atoms.at(ja).position;
                V += _polymer_potential->vij(ri, rj);
            }
        }
    }
    // INTERMOLECULAR
    // solvent-solvent
    for(int is = 0; is < state.nsolvents() - 1; ++is){
        ri = state.solvents.at(is).atom.position;
        for(int js = is + 1; js < state.nsolvents(); ++js){
            rj = state.solvents.at(js).atom.position;
            V += _solvent_potential->vij(ri, rj);
        }
    }
    for(int im = 0; im < state.nm(); ++im){
        // polymer-polymer
        for(int jm = im + 1; jm < state.nm(); ++jm){
            for(int ia = 0; ia < na; ++ia){
                ri = state.polymers.at(im).atoms.at(ia).position;
                for(int ja = 0; ja < na; ++ja){
                    rj = state.polymers.at(im).atoms.at(ja).position;
                    V += _polymer_potential->vij(ri, rj);
                }
            }
        }
        // polymer-solvent
        for(int ia = 0; ia < na; ++ia){
            ri = state.polymers.at(im).atoms.at(ia).position;
            for(int js = 0; js < state.nsolvents(); ++js){
                rj = state.solvents.at(js).atom.position;
                V += _inter_potential->vij(ri, rj);
            }
        }
    }
    return V;
}
