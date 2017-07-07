// 2017 Artur Avkhadiev
/*! \file molecule.h
*/
#ifndef POLYMER_MOLECULE_H
#define POLYMER_MOLECULE_H
#include <vector>
#include <map>
#include <string>
#include "atom.h"
#include "bond.h"
typedef struct molecule_t {
    int na;                         /*>> number of atoms */
    int nb;                         /*>> number of fixed bonds */
    std::vector<Atom> atoms;        /*>> list of atoms */
    std::vector<Bond> bonds;        /*>> list of bonds */
} Molecule;
Molecule initialize_molecule(std::vector<Atom> atoms, std::vector<Bond> bonds);
/**
* Checks if a molecule is valid:
*   that number of contrained bonds is less than or equal to the number of atoms;
*   that all bonds in the molecule are valid.
*   Throws an invalid argument exception if the molecule is invalid;
*   Returns silently otherwise
*/
void check_molecule(Molecule m);
/**
* Takes a molecule and ensures all time records on positions and velocities coincide.
* If the second argument is specified, ensures these times also equal to the
* second argument
*/
bool is_time_consistent(Molecule molecule, double time = -1);
/**
* Takes a pointer to a molecule and a time and sets all records of positions and velocities
* of all the atoms in the to that time
*/
void set_time(Molecule *molecule, double time);
/**
* Takes a molecule and ouptuts its std::string representation
*/
std::string molecule_to_string(Molecule m, bool verbose = true);
::std::ostream& operator<<(::std::ostream& os, const Molecule& m);
/**
* Takes a non-verbose represenation of a molecule from molecule_to_string and
* returns the corresponding molecule
*/
Atom string_to_atoms(std::string nonverbose_str);
#endif
