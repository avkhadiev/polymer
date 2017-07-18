// 2017 Artur Avkhadiev
/*! \file molecule.cpp
*/
#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <algorithm>                        /* std::find */
#include <iterator>                         /* std::distance */
#include <../include/atom.h>
#include <../include/bond.h>
#include <../include/molecule.h>
Molecule initialize_molecule(std::vector<Atom> atoms, std::vector<Bond> bonds, double t) {
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
            // check if all bonds have pairs of different atom indices,
            // and that no two bonds have the same pair
            check_bonds(bonds, atoms);
            molecule.na = na;
            molecule.nb = nb;
            molecule.atoms = atoms;
            molecule.bonds = bonds;
            // set the time of all atoms' positions and velocities to t
            set_time(&molecule, t);
            // update all bonds' current length squared according to the atoms'
            // positions at time t.
            Bond *bond_ptr;
            for(int i = 0; i < nb; ++i) {
                bond_ptr = &(molecule.bonds.at(i));
                try {
                    update_bond_length(bond_ptr, molecule.atoms);
                }
                catch (std::invalid_argument &e) {
                    throw;
                }
            }
        }
        catch(std::invalid_argument &e){
            throw;
        }
    }
    return molecule;
}
void check_molecule(Molecule m) {
    bool valid = true;
    std::string na_declared;
    std::string na_actual;
    std::string nb_declared;
    std::string nb_actual;
    std::string err_msg = "check_molecule: ";
    std::string new_msg;
    if(m.na != m.atoms.size()){
        valid = false;
        na_declared = std::to_string(m.na);
        na_actual = std::to_string(m.atoms.size());
        new_msg = "number of atoms ("
            + na_actual
            + ") is not equal to the one declared ("
            + na_declared
            + "); ";
        err_msg += new_msg;
    }
    if(m.nb != m.bonds.size()){
        valid = false;
        nb_declared = std::to_string(m.nb);
        nb_actual = std::to_string(m.bonds.size());
        new_msg = "number of bonds ("
            + nb_actual
            + ") is not equal to the one declared ("
            + nb_declared
            + "); ";
        err_msg += new_msg;
    }
    try
    {
        check_bonds(m.bonds, m.atoms);
    }
    catch (std::invalid_argument &e)
    {
        valid = false;
        std::string check_bonds_err(e.what());
        err_msg += check_bonds_err;
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
            fprintf(stderr, "atoms time records do not match the required time (%.3f)\n", time);
            consistent = false;
        }
        // otherwise set the time argument as not optional and equal to the
        // time record
        else {
            time = some_time_record;
            for(int i = 0; i < molecule.na; ++i) {
                consistent = consistent && is_time_consistent(molecule.atoms.at(i), time);
                if (!consistent) {
                    fprintf(stderr, "atom %d is not time-consistent\n", i);
                    break;
                }
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
        // verbose output
        if (verbose) {
            m_str_atoms += std::to_string(i)
                + " "
                +  atom_to_string(m.atoms.at(i), verbose)
                + "\n";
        }
        else {
        // non-verbose output
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
        // verbose output
        m_str = m_str_header
            + m_str_na_header + m_str_na + "\n"
            + m_str_nb_header + m_str_nb + "\n"
            + m_str_atoms + m_str_bonds;
    }
    else {
        // non-verbose output
        m_str = m_str_na + " " + m_str_nb + "\n" + m_str_atoms + m_str_bonds;
    }
    return m_str;
}
::std::ostream& operator<<(::std::ostream& os, const Molecule& m) {
    return os << molecule_to_string(m).c_str();
}
Molecule string_to_molecule(std::ifstream& input_stream){
    // data format in non-verbose lines:
    // na nb
    // ... na atoms strings...
    // ... nb bond strings...
    // get all lines in a vector of strings
    std::string line;
    std::getline(input_stream, line);
    // get number of atoms and number of bonds
    std::istringstream ss(line.c_str());
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> words(begin, end);
#ifdef DEBUG
    printf("%-45s%-30s%-15s\n",
        "string_to_molecule reads:",
        "number of atoms",
        "number of bonds");
#endif
    int na = atoi(words.at(0).c_str());
    int nb = atoi(words.at(1).c_str());
#ifdef DEBUG
    printf("%-45s%-30d%-15d\n",
        "data:",
        na,
        nb);
#endif
    // read all the atoms from the input stream
    std::vector<Atom> atoms;
    Atom next_atom;
    for(int nlines = 0; nlines < na; ++nlines) {
        std::getline(input_stream, line);
        next_atom = string_to_atom(line);
        atoms.push_back(next_atom);
    }
    // read all the bonds from the input stream
    std::vector<Bond> bonds;
    Bond next_bond;
    for(int nlines = 0; nlines < nb; ++nlines) {
        std::getline(input_stream, line);
        try {
            next_bond = string_to_bond(line);
        }
        catch (std::invalid_argument &e) {
            throw;
        }
        bonds.push_back(next_bond);
    }
    // initialize the molecule
    Molecule m;
    try
    {
        m = initialize_molecule(atoms, bonds);
        return m;
    }
    catch (std::invalid_argument &e)
    {
        throw;
    }
}
