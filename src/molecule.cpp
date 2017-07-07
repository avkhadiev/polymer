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
            molecule.atoms = atoms;
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
        check_bonds(m.bonds, m.atoms);
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
bool is_time_consistent(Molecule molecule, double time) {
    bool consistent = true;
    // if atoms list is empty, we are done; otherwise:
    if (!molecule.atoms.empty()){
        // if time was not optional and is not equal to one of the time values
        double some_time_record = molecule.atoms.back().position.second;
        if (time != -1 && time != some_time_record){
                consistent = false;
        }
        // otherwise set the time argument as not optional and equal to the
        // time record
        else {
            time = some_time_record;
            for(Atom& atom : molecule.atoms) {
                consistent = consistent && is_time_consistent(atom, time);
            }
        }
    }
    return consistent;
}
void set_time(Molecule *molecule, double time) {
    for (int i = 0; i < molecule->atoms.size(); ++i) {
        set_time(&(molecule->atoms.at(i)), time);
    }
}
std::string molecule_to_string(Molecule m, bool verbose){
    try
    {
        check_molecule(m);
    }
    catch (std::invalid_argument &e)
    {
        throw;
    }
    std::string m_str_header = "MOLECULE:\n";
    std::string m_str_na_header = "NA = ";
    std::string m_str_na = std::to_string(m.na);
    std::string m_str_nb_header = "NB = ";
    std::string m_str_nb = std::to_string(m.nb);
    std::string m_str_atoms = "";
    for (int i = 0; i < m.atoms.size(); ++i) {
        if (verbose) {
            m_str_atoms += std::to_string(i)
                + " "
                +  atom_to_string(m.atoms.at(i), verbose)
                + "\n";
        }
        else {
            m_str_atoms += atom_to_string(m.atoms.at(i), verbose) + "\n";
        }
    }
    Bond b;
    std::string m_str_bonds = "";
    for (int i = 0; i < m.bonds.size(); ++i) {
        b = m.bonds.at(i);
        // find the indicthe bond string
        if (verbose) {
            m_str_bonds += std::to_string(i)
                + " "
                + bond_to_string(b, verbose)
                + "\n";
        }
        else {
            m_str_bonds += bond_to_string(b, verbose) + "\n";
        }
    }
    std::string m_str;
    if (verbose) {
        m_str = m_str_header
            + m_str_na_header + m_str_na + "\n"
            + m_str_nb_header + m_str_nb + "\n"
            + m_str_atoms + m_str_bonds;
    }
    else {
        m_str = m_str_na + "\n" + m_str_nb + "\n" + m_str_atoms + m_str_bonds;
    }
    return m_str;
}
::std::ostream& operator<<(::std::ostream& os, const Molecule& m) {
    return os << molecule_to_string(m).c_str();
}
