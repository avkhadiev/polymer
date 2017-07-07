// 2017 Artur Avkhadiev
/*! \file bond.h
*/
#ifndef POLYMER_BOND_H
#define POLYMER_BOND_H
#include <vector>
#include <utility>                  /* std::pair, std::make_pair */
#include "atom.h"
struct bond_t {
    double fixed_length_sq;
    std::pair<double, double> current_length_sq;
    int atom1;
    int atom2;
    bool operator==(const bond_t& b) const
    {
        bool atoms_equal = ((atom1 == b.atom1 && atom2 == b.atom2)
            || (atom1 == b.atom2 && atom2 == b.atom1));
        return (fixed_length_sq == b.fixed_length_sq
            && current_length_sq == b.current_length_sq
            && atoms_equal);
    }
    bool operator!=(const bond_t& b) const
    {
        bool atoms_equal = ((atom1 == b.atom1 && atom2 == b.atom2)
            || (atom1 == b.atom2 && atom2 == b.atom1));
        return (fixed_length_sq != b.fixed_length_sq
            || current_length_sq != b.current_length_sq
            || !atoms_equal);
    }
};
typedef struct bond_t Bond;
/**
* Takes the bond and the relevant vector of atoms.
* If both atoms have their position vectors defined at the same time, updates
* the square of the bond length.
* If atoms have their positions defined at different times, raises an invalid
* argument exception.
*/
void update_bond_length(Bond *bond, std::vector<Atom> atoms);
Bond initialize_bond(int atom1, int atom2, double fixed_length);
/**
* Takes a bond and ouptuts its std::string representation
*/
std::string bond_to_string(Bond bond, bool verbose = true);
::std::ostream& operator<<(::std::ostream& os, const Bond& bond);
/*
* Checks if the bond has different atom1 and atom2 and that these indices are
* not greater than the size of the atom vector.
*/
void check_bond(Bond bond, std::vector<Atom> atoms);
/*
* Performs check_bond(bond, atoms) for each bond in bonds
*/
void check_bonds(std::vector<Bond> bonds, std::vector<Atom> atoms);
#endif
