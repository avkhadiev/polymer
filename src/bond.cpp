// 2017 Artur Avkhadiev
/*! \file bond.cpp
*/
#include <algorithm>                /* std::find */
#include <utility>                  /* std::pair, std::make_pair */
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>                    /* pow */
#include <../include/atom.h>
#include <../include/bond.h>
#include <../include/vector.h>
void update_bond_length(Bond *bond, std::vector<Atom> atoms) {
    std::pair<Vector, double> r1_t = atoms.at(bond->atom1).position;
    std::pair<Vector, double> r2_t = atoms.at(bond->atom2).position;
    // if two atoms have positions defined at the same time
    if (r1_t.second == r2_t.second) {
        Vector r1 = r1_t.first;
        Vector r2 = r2_t.first;
        Vector r12 = subtract(r1, r2);
        double l = normsq(r12);
        double t = r1_t.second;
        std::pair<double, double> current_length_sq = std::make_pair(l, t);
        bond->current_length_sq = current_length_sq;
    }
    else {
        // positions of atoms were defined at different times
        // throw exception, keep current length
        std::string err_msg = "update_bond_length: positions of a1 and a2 are defined at different times";
        throw std::invalid_argument( err_msg );
    }
    return;
}
Bond initialize_bond(int atom1, int atom2, double fixed_length_sq) {
    // check if two atoms are the same
    if ( atom1 == atom2 ){
        std::string err_msg = "check_bond: two atoms are the same";
        throw std::invalid_argument( err_msg );
    }
    else {
        Bond bond = {.atom1 = atom1, .atom2 = atom2,
            .fixed_length_sq = fixed_length_sq,
            .current_length_sq = std::make_pair(-1, -1)};
        return bond;
    }
}
std::string bond_to_string(Bond bond, bool verbose){
    std::string b_str_header = "BOND: \n";
    std::string b_atom1_str = std::to_string(bond.atom1);
    std::string b_atom2_str = std::to_string(bond.atom2);
    std::string b_fixed_sq_header = "fixedsq = ";
    std::string b_fixed_sq_value = std::to_string(bond.fixed_length_sq);
    std::string b_current_sq_header = "currentsq:";
    std::string b_current_sq_value_header = "lsq = ";
    std::string b_current_sq_value = std::to_string(bond.current_length_sq.first);
    std::string b_current_sq_time_header = "t = ";
    std::string b_current_sq_time = std::to_string(bond.current_length_sq.second);
    std::string b_str;
    if (verbose) {
        // verbose output
        b_str = b_atom1_str + "-" + b_atom2_str + " " + b_str_header
            + b_fixed_sq_header + b_fixed_sq_value + "\n"
            + b_current_sq_header + " "
                + b_current_sq_time_header + b_current_sq_time + " "
                + b_current_sq_value_header + b_current_sq_value;
    }
    else {
        // non-verbose output
        b_str = b_atom1_str + " " + b_atom2_str + + " " + b_fixed_sq_value;
    }
    return b_str;
}
::std::ostream& operator<<(::std::ostream& os, const Bond& bond) {
    return os << bond_to_string(bond).c_str();
};
Bond string_to_bond(std::string nonverbose_str){
    std::istringstream ss(nonverbose_str.c_str());
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> words(begin, end);
    // convert strings to data
    int atom1 = atoi(words.at(0).c_str());
    int atom2 = atoi(words.at(1).c_str());
    double fixed_length_sq = atof(words.at(2).c_str());
    try
    {
        Bond b = initialize_bond(atom1, atom2, fixed_length_sq);
        return b;
    }
    catch (std::invalid_argument &e)
    {
        throw;
    }
}
void check_bond(Bond bond, std::vector<Atom> atoms) {
    if ( bond.atom1 == bond.atom2 ){
        std::string err_msg = "check_bond: two atoms are the same";
        throw std::invalid_argument( err_msg );
    }
    else{
        if (bond.atom1 >= atoms.size() || bond.atom2 >= atoms.size()){
            std::string err_msg = "check_bond: at least one of the atoms is not valid";
            throw std::invalid_argument( err_msg );
        }
    }
    return;
}
void check_bonds(std::vector<Bond> bonds, std::vector<Atom> atoms) {
    for (Bond bond : bonds) {
        try
        {
            check_bond(bond, atoms);
        }
        catch (std::invalid_argument &e)
        {
            throw;
        }
        // check if any two bonds are pointing to the same two atoms
        for (int i = 0; i < bonds.size() - 1; ++i){
            for (int j = i + 1; j < bonds.size(); ++j) {
                // if the pairs of atoms in the bonds are identical
                if (((bonds.at(i).atom1 == bonds.at(j).atom1)
                        && (bonds.at(i).atom2 == bonds.at(j).atom2))
                    || ((bonds.at(i).atom1 == bonds.at(j).atom2)
                        && (bonds.at(i).atom2 == bonds.at(j).atom1))) {
                    std::string err_message = "check_bonds: two bonds are pointing to the same two atoms";
                    throw std::invalid_argument( err_message );
                }
            }
        }
    }
};
