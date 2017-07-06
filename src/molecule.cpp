// 2017 Artur Avkhadiev
/*! \file molecule.cpp
*/
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>                        /* std::find */
#include <iterator>                         /* std::distance */
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
            // order the atoms
            std::map<int, Atom> ordered_atoms;
            for (int i = 0; i < atoms.size(); ++i) {
                ordered_atoms.emplace(i, atoms.at(i));
            }
            molecule.atoms = ordered_atoms;
            molecule.bonds = bonds;
        }
        catch(std::invalid_argument &e){
            throw;
        }
    }
    return molecule;
}
void check_molecule(Molecule m) {
    bool valid = true;
    std::string err_msg = "check_molecule: ";
    if(m.na != m.atoms.size()){
        valid = false;
        err_msg += "number of atoms is not equal to the one declared; ";
    }
    if(m.nb != m.bonds.size()){
        valid = false;
        err_msg += "number of bonds is not equal to the one declare; ";
    }
    try
    {
        // extract values from the map of atoms keyed on their number
        std::vector<Atom> atom_vector;
        for( std::map<int, Atom>::iterator it = m.atoms.begin();
            it != m.atoms.end();
            ++it ) {
            atom_vector.push_back( it->second );
    }
        check_bonds(m.bonds, atom_vector);
    }
    catch (std::invalid_argument &e)
    {
        valid = false;
        err_msg += "check_bonds threw an exception.";
    }
    if (!valid) {
        throw std::invalid_argument( err_msg );
    }
    return;
}
std::string molecule_to_string(Molecule m){
    try
    {
        check_molecule(m);
    }
    catch (std::invalid_argument &e)
    {
        throw;
    }
    std::string m_str_header = "MOLECULE:\n";
    std::string m_str_na = "NA: " + std::to_string(m.na) + "\n";
    std::string m_str_nb = "NB: " + std::to_string(m.nb) + "\n";
    std::string m_str_atoms = "";
    for (int i = 0; i < m.atoms.size(); ++i) {
        m_str_atoms = std::to_string(i)
            + " "
            +  atom_to_string(m.atoms.at(i))
            + "\n";
    }
    Bond b;
    std::string m_str_bonds = "";
    for (int i = 0; i < m.bonds.size(); ++i) {
        // find the indices of the two atoms that the bond is pointing to
        b = m.bonds.at(i);
        Atom atom1 = *(b.atom1);
        Atom atom2 = *(b.atom2);
        std::vector<Atom>::iterator it1, it2;
        it1 = std::find(m.atoms.begin(), m.atoms.end(), atom1);
        it2 = std::find(m.atoms.begin(), m.atoms.end(), atom2);
        if (it1 != m.atoms.end() && it2 != m.atoms.end()){
            // get the atom indices
            int index1 = 1 + std::distance(m.atoms.begin(), it1);
            int index2 = 1 + std::distance(m.atoms.begin(), it2);
            if (index1 != index2) {
            std::string atoms_connected = std::to_string(index1) + "-" + std::to_string(index2);
            // append the bond string
            m_str_bonds = std::to_string(i + 1)
                + " "
                + atoms_connected
                + " "
                + bond_to_string(b)
                + "\n";
            }
            else {
                // atom1 and atom2 are the same atom!
                std::string err_msg = "molecule_to_string: bond in molecule is not valid, contains same-atom pointers";
                throw std::invalid_argument( err_msg );
            }
        }
        else {
            // bond was pointing to atom that was not in the list;
            // throw exception
            std::string err_msg = "molecule_to_string: bond in molecule is not valid, points to atom not in molecule";
            throw std::invalid_argument( err_msg );
        }
    }
    std::string m_str = m_str_header + m_str_na + m_str_nb + m_str_atoms + m_str_bonds;
    return m_str;
}
::std::ostream& operator<<(::std::ostream& os, const Molecule& m) {
    return os << molecule_to_string(m).c_str();
}
