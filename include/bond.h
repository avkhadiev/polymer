// 2017 Artur Avkhadiev
/*! \file bond.h
*/
#ifndef POLYMER_BOND_H
#define POLYMER_BOND_H
#include <vector>
#include "atom.h"
struct bond_t {
    double fixed_length_sq;
    std::pair<double, double> current_length_sq;
    Atom *atom1;
    Atom *atom2;
    bool operator==(const bond_t& b) const
    {
        return (fixed_length_sq == b.fixed_length_sq
            && current_length_sq == b.current_length_sq
            && *atom1 == *(b.atom1)
            && *atom1 == *(b.atom1));
    }
    bool operator!=(const bond_t& b) const
    {
        return (fixed_length_sq != b.fixed_length_sq
            || current_length_sq != b.current_length_sq
            || *atom1 != *(b.atom1)
            || *atom1 != *(b.atom1));
    }
};
typedef struct bond_t Bond;
/**
* Takes the bond and and follows the references to the two atoms.
* If both atoms have their position vectors defined at the same time, updates
* the square of the bond length.
* If atoms have their positions defined at different times, raises an invalid
* argument exception.
*/
void update_bond_length(Bond *bond);
Bond initialize_bond(Atom *atom1, Atom *atom2, double fixed_length);
/**
* Takes a bond and ouptuts its std::string representation
*/
std::string bond_to_string(Bond bond);
::std::ostream& operator<<(::std::ostream& os, const Bond& bond);
/*
* Checks if the atom pointers in the given bond point to one of the atoms in
* the given vector, and that the pointers don't point to the same atom.
*/
void check_bond(Bond bond, std::vector<Atom> atoms);
/*
* Performs check_bond(bond, atoms) for each bond in bonds
*/
void check_bonds(std::vector<Bond> bonds, std::vector<Atom> atoms);
#endif
