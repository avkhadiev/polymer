// 2017 Artur Avkhadiev
/*! \file atom.h
*/
#ifndef POLYMER_ATOM_H
#define POLYMER_ATOM_H
#include "vector.h"
typedef struct atom_t {
    double mass;
    Vector position;
    Vector velocity;
    Vector acceleration;
} Atom;
#endif
