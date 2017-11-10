// 2017 Artur Avkhadiev
/*! \file simple_solvent.h
* A simple::Solvent consists in a single atom, a position and a velocity vector.
* All solvents have the same mass.
*/
#ifndef POLYMER_SIMPLE_SOLVENT_H
#define POLYMER_SIMPLE_SOLVENT_H
#include <string>
#include <iostream>
#include "vector.h"
#include "simple_atom.h"
namespace simple {
    class Solvent {
    protected:
        // all solvents in the simulation have the same mass
        static double _m;
    public:
        Atom atom;
        // operators
        bool operator==(const Solvent &other) const;
        bool operator!=(const Solvent &other) const{return !(*this == other);};
        // getters
        static double m() {return _m;}
        Vector r() const {return atom.position;};
        Vector v() const {return atom.velocity;};
        Vector f() const {return atom.force;};
        // setters
        static void set_m(double m) {_m = m;}
        void set_r(Vector r) {atom.position = r;};
        void set_v(Vector v) {atom.velocity = v;};
        void set_f(Vector f) {atom.force = f;};
        // I/O
        virtual std::string to_string(bool verbose = true) const;
        // creation/annihilation
        Solvent(Vector position = vector(0.0, 0.0, 0.0),
                Vector velocity = vector(0.0, 0.0, 0.0));
        virtual ~Solvent();
    };
    /**
    * Takes a solvent and ouptuts its std::string representation
    */
    ::std::ostream& operator<<(::std::ostream& os,
        const simple::Solvent& polymer);
    /**
    * Reads the stream from the current position, advancing the position until after the molecule instance; converts input data to a molecule instance
    */
    simple::Solvent string_to_solvent(std::ifstream& input_stream);
} // namespace simple
#endif
