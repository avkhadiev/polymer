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
Atom initialize_atom(double mass, Vector r, Vector v, double t = -1.0);
/**
* Takes an atom and ensures all time records on positions and velocities coincide.
* If the second argument is specified, ensures these times also equal to the
* second argument
*/
bool is_time_consistent(Atom atom, double time = -1);
/**
* Takes a pointer to atom and a time and sets all records of positions and velocities
* in the atom to that time
*/
void set_time(Atom *atom, double time);
/**
* Takes an atom and ouptuts its std::string representation
* the representaion is verbose by default
*/
std::string atom_to_string(Atom a, bool verbose = true);
::std::ostream& operator<<(::std::ostream& os, const Atom& a);
/**
* Takes a non-verbose represenation of an atom from atom_to_string and
* returns the corresponding atom. If time is specified, sets the time of
* positions and velocities to t.
*/
Atom string_to_atom(std::string nonverbose_str, double t = -1.0);
#endif
