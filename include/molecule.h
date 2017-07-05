// 2017 Artur Avkhadiev
/*! \file molecule.h
*/
#ifndef POLYMER_MOLECULE_H
#define POLYMER_MOLECULE_H
#include <vector>
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
* Takes a molecule and ouptuts its std::string representation
*/
std::string molecule_to_string(Molecule m);
::std::ostream& operator<<(::std::ostream& os, const Molecule& m);
#endif
