// 2017 Artur Avkhadiev
/*! \file molecule.h
*/
#ifndef POLYMER_MOLECULE_H
#define POLYMER_MOLECULE_H
#include <vector>
#include "atom.h"
#include "bond.h"
typedef struct molecule_t {
    int na;                         /*>> number of atoms */
    int nb;                         /*>> number of fixed bonds */
    std::vector<Atom> atoms;       /*>> list of atoms */
    std::vector<Bond> bonds;       /*>> list of bonds */
} Molecule;
#endif
