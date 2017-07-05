// 2017 Artur Avkhadiev
/*! \file molecule.cpp
*/
#include <vector>
#include <string>
#include <stdexcept>
#include <../include/atom.h>
#include <../include/bond.h>
#include <../include/molecule.h>
Molecule initialize_molecule(std::vector<Atom> atoms, std::vector<Bond> bonds) {
    Molecule molecule;
    int na = atoms.size();
    int nb = bonds.size();
    // number of fixed bonds cannot exceed the number of atoms
    // (it is actually one less than the number of atoms,
    // for non-cyclic molecules)
    if (nb > na) {
        std::string err_msg = "initialize_molecule: number of bonds exceeds number of atoms ";
        throw std::invalid_argument( err_msg );
    }
    else{
        try{
            check_bonds(bonds, atoms);
            molecule.na = na;
            molecule.nb = nb;
            molecule.atoms = atoms;
            molecule.bonds = bonds;
        }
        catch(std::invalid_argument &e){
            throw;
        }
    }
    return molecule;
}
std::string molecule_to_string(Molecule m){
    std::string m_str_header = "MOLECULE:\n";
    std::string m_str_na = "NA: " + std::to_string(m.na) + "\n";
    std::string m_str_nb = "NB: " + std::to_string(m.nb) + "\n";
    std::string m_str_atoms = "";
    for (Atom& atom : m.atoms) {
        m_str_atoms += atom_to_string(atom) + "\n";
    }
    std::string m_str_bonds = "";
    for (Bond& bond : m.bonds) {
        m_str_bonds += bond_to_string(bond) + "\n";
    }
    std::string m_str = m_str_header + m_str_na + m_str_nb + m_str_atoms + m_str_bonds;
    return m_str;
}
::std::ostream& operator<<(::std::ostream& os, const Molecule& m) {
    return os << molecule_to_string(m).c_str();
}
