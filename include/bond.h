// 2017 Artur Avkhadiev
/*! \file bond.h
*/
#ifndef POLYMER_BOND_H
#define POLYMER_BOND_H
#include "atom.h"
typedef struct bond_t {
    double fixed_length;
    Atom& atom1;
    Atom& atom2;
} Bond;
#endif
