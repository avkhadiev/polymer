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
    std::map<int, Atom> atoms;      /*>> ordered atoms */
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
* Takes a molecule and ouptuts its std::string representation
*/
std::string molecule_to_string(Molecule m);
::std::ostream& operator<<(::std::ostream& os, const Molecule& m);
#endif
