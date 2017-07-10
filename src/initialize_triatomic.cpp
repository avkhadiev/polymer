// 2017 Artur Avkhadiev
/*! \file initialize_triatomic.cpp
*/
#include <string>
#include <utility>              /* std::pair, std::make_pair */
#include <stdexcept>
#include <cmath>                /* pow */
#include <gtest/gtest.h>
#include "../include/vector.h"
#include "../include/atom.h"
#include "../include/bond.h"
#include "../include/molecule.h"
#include "../include/state.h"

int main(int argc, char **argv){
    std::string fname;
    if (argc == 1){
        fname = std::string("test_state");
    }
    else {
        fname = std::string(argv[1]);
    }
    static const int na = 3;
    static const int nb = 2;
    static const int nm = 2;
    double t;
    Vector r[na];
    Vector v[na];
    double mass[na];
    double fixed_length_sq;
    Atom a[na];
    std::vector<Atom> atoms;
    Bond b[nb];
    std::vector<Bond> bonds;
    Molecule m[nm];
    std::vector<Molecule> molecules;
    /**
    /* SET UP
    */
    t = 3.0;
    // initialize atoms
    // (-1, -1, -1), (0, 0, 0), and (1, 1, 1)
    for (int i = 0; i < na;  ++i) {
        double coord = i - 1.0;
        mass[i] = 1.0;
        r[i] = vector(coord, coord, coord);
        v[i] = vector(coord, coord, coord);
        a[i] = initialize_atom(mass[i], r[i], v[i], t);
        atoms.push_back(a[i]);
    }
    // initialize bonds
    fixed_length_sq = 3.0;
    for (int i = 0; i < nb; ++i) {
        try
        {
            b[i] = initialize_bond(i, i + 1, fixed_length_sq);
        }
        catch (std::invalid_argument &e)
        {
            throw;
        }
        bonds.push_back(b[i]);
    }
    // initialize molecules
    for (int i = 0; i < nm; ++i) {
        try
        {
            m[i] = initialize_molecule(atoms, bonds, t);
        }
        catch (std::invalid_argument &e)
        {
            throw;
        }
        molecules.push_back(m[i]);
    }
    /**
    * INITIALIZE STATE AND WRITE OUT TO FILE
    */
    State s = initialize_state(molecules, t);
    std::string fname_verbose = fname + "_verbose";
    std::string fname_nonverbose = fname + "_nonverbose";
    write_state_to_file(s, fname_verbose, true);
    write_state_to_file(s, fname_nonverbose, false);
}
