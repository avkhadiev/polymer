// 2017 Artur Avkhadiev
/*! \file simple_bond.h
*
*/
#ifndef POLYMER_SIMPLE_BOND_H
#define POLYMER_SIMPLE_BOND_H
#include <string>
#include <iostream>
#include "vector.h"
#include "simple_atom.h"
namespace simple {
    class Bond {
    public:
        Vector position;
        Vector velocity;
        Bond& operator=(const Bond &other);
        bool operator==(const Bond &other) const;
        bool operator!=(const Bond &other) const;
        /**
        * Takes an bond and ouptuts its std::string representation
        * Verbose representation is multi-line with a header
        * Non-verbose representation is single-line
        */
        std::string to_string(bool verbose = true) const;
        /**
        * Direct initialization via position and velocity vectors
        * Position vector will be normalized
        */
        Bond(Vector new_position = vector(1.0, 0.0, 0.0),
            Vector new_velocity = vector(0.0, 0.0, 0.0));
        /**
        * Position is normalized difference of atomic positions (end - start)
        * Velocity is difference of atomic velocities (end - start)
        */
        Bond(Atom start, Atom end);
        ~Bond();
    };
    /**
    * Takes a non-verbose represenation of a simple::Bond from bond_to_string
    * and returns the corresponding bond.
    */
    Bond string_to_bond(std::string nonverbose_str);
    ::std::ostream& operator<<(::std::ostream& os, const simple::Bond& b);
} // namespace simple
#endif
