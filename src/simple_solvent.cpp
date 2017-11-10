// 2017 Artur Avkhadiev
/*! \file simple_solvent.cpp
*/
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "../include/vector.h"
#include "../include/default_macros.h"
#include "../include/simple_solvent.h"
namespace simple {
    double Solvent::_m = DEFAULT_M;
    Solvent::Solvent(Vector r, Vector v) :
        atom (Atom(r, v)){}
    Solvent::~Solvent(){}
    bool Solvent::operator==(const Solvent &other) const{
        return atom == other.atom;
    }
    std::string Solvent::to_string(bool verbose) const{
        std::string m_str_name = "SIMPLE SOLVENT:";
        std::string m_str_atoms = atom.to_string(verbose) + "\n";
        std::string m_str;
        if (verbose) {
            // verbose output
            m_str = m_str_name + "\n" + m_str_atoms;
        }
        else {
            // non-verbose output
            m_str = m_str_atoms;
        }
        return m_str;
    }
    ::std::ostream& operator<<(::std::ostream& os,
        const simple::Solvent& molecule){
        return os << molecule.to_string().c_str();
    }
    Solvent string_to_solvent(std::ifstream& input_stream){
        Solvent molecule;
        std::string line;
        std::getline(input_stream, line);
        molecule.atom = string_to_atom(line);
        return molecule;
    }
} // namespace simple
