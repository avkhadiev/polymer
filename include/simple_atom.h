// 2017 Artur Avkhadiev
/*! \file simple_atom.h
* Unlike the regular atom, a simple atom not store any time records:
*   storing time records together with velocity, position and force vectors was
*   useful for implementing and ensuring the correctness of integrators.
*   stripping simple atoms of time records will simplify the code and speed
*   up its execution.
*   This version of simple atoms will not store mass records. It is unnecessary
*   for current purposes, all atoms are the same mass and belong to a polymer
*   molecule; therefore, it makes sense to store the mass record elsewhere.
*/
#ifndef POLYMER_SIMPLE_ATOM_H
#define POLYMER_SIMPLE_ATOM_H
#include <string>
#include <iterator>
#include <iostream>
#include "vector.h"
namespace simple {
    class Atom {
    public:
        Vector position;
        Vector velocity;
        Vector force;
        Atom& operator=(const Atom &other);
        bool operator==(const Atom &other) const;
        bool operator!=(const Atom &other) const;
        /**
        * Takes an atom and ouptuts its std::string representation
        * In non-verbose representation, forces are not displayed.
        */
        std::string to_string(bool verbose = true) const;
        Atom(Vector position = vector(0.0, 0.0, 0.0),
            Vector velocity = vector(0.0, 0.0, 0.0),
            Vector force = vector(0.0, 0.0, 0.0));
        ~Atom();
    };
    /**
    * Takes a non-verbose represenation of a simple::Atom from atom_to_string
    * and returns the corresponding atom.
    */
    Atom string_to_atom(std::string nonverbose_str);
    ::std::ostream& operator<<(::std::ostream& os, const simple::Atom& a);
} // namespace simple
#endif
