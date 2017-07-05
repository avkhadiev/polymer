// 2017 Artur Avkhadiev
/*! \file bond.cpp
*/
#include <algorithm>                /* std::find */
#include <utility>                  /* std::pair, std::make_pair */
#include <stdexcept>
#include <iostream>
#include <string>
#include <cmath>                    /* pow */
#include <../include/atom.h>
#include <../include/bond.h>
#include <../include/vector.h>
void update_bond_length(Bond *bond) {
    std::pair<Vector, double> r1_t = bond->atom1->position;
    std::pair<Vector, double> r2_t = bond->atom2->position;
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
Bond initialize_bond(Atom *atom1, Atom *atom2, double fixed_length) {
    // check if two atoms are the same
    if ( *atom1 == *atom2 ){
        std::string err_msg = "check_bond: two atoms are the same";
        throw std::invalid_argument( err_msg );
    }
    else {
        double fixed_length_sq = pow(fixed_length, 2.0);
        Bond bond = {.atom1 = atom1, .atom2 = atom2,
            .fixed_length_sq = fixed_length_sq,
            .current_length_sq = std::make_pair(-1, -1)};
        try
        {
            update_bond_length(&bond);
            return bond;
        }
        catch (std::invalid_argument &e)
        {
            // positions of atoms were defined at different times
            // cannot initialzie, throw exception
            throw;
        }
    }
}
std::string bond_to_string(Bond bond){
    std::string b_str_header = "BOND: \n";
    std::string b_fixed_sq_str = "fixedsq: " + std::to_string(bond.fixed_length_sq) + "\n";
    std::string b_current_sq_value_str = std::to_string(bond.current_length_sq.first);
    std::string b_current_sq_time_str = std::to_string(bond.current_length_sq.second);
    std::string b_current_sq_str = "currentsq: t = "
        + b_current_sq_time_str
        + " lsq: " + b_current_sq_value_str;
    std::string b_str = b_str_header
        + b_fixed_sq_str
        + b_current_sq_str;
    return b_str;
}
::std::ostream& operator<<(::std::ostream& os, const Bond& bond) {
    return os << bond_to_string(bond).c_str();
};
void check_bond(Bond bond, std::vector<Atom> atoms) {
    Atom atom1 = *(bond.atom1);
    Atom atom2 = *(bond.atom2);
    if ( atom1 == atom2 ){
        std::string err_msg = "check_bond: two atoms are the same";
        throw std::invalid_argument( err_msg );
    }
    else{
        bool found1 = (std::find(atoms.begin(), atoms.end(), atom1) != atoms.end());
        bool found2 = (std::find(atoms.begin(), atoms.end(), atom2) != atoms.end());
        if (!found1 || !found2){
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
    }
};
