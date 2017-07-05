// 2017 Artur Avkhadiev
/*! \file atom.h
*/
#ifndef POLYMER_ATOM_H
#define POLYMER_ATOM_H
#include <utility>      /* std::pair, std::make_pair */
#include <iostream>
#include "vector.h"
struct atom_t {
    double mass;
    std::pair<Vector, double> position;
    std::pair<Vector, double> velocity;
    std::pair<Vector, double> force;
    bool operator==(const atom_t& a) const
    {
        return (mass == a.mass
            && position == a.position
            && velocity == a.velocity
            && force == a.force);
    }
    bool operator!=(const atom_t& a) const
    {
        return (mass == a.mass
            || position == a.position
            || velocity == a.velocity
            || force == a.force);
    }
};
typedef struct atom_t Atom;
Atom initialize_atom(double mass, Vector r, Vector v, double t = 0);
/**
* Takes an atom and ouptuts its std::string representation
*/
std::string atom_to_string(Atom a);
::std::ostream& operator<<(::std::ostream& os, const Atom& a);
#endif
