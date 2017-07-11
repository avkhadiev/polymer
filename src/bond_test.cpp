// 2017 Artur Avkhadiev
/*! \file bond_test.cpp
*/
#include <utility>                  /* std::pair, std::make_pair */
#include <cmath>                    /* pow */
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <gtest/gtest.h>
#include "../include/vector.h"
#include <../include/atom.h>
#include <../include/bond.h>

class BondTest : public ::testing::Test {
 protected:
     Vector r1;
     Vector r2;
     Vector r3;
     Vector r4;
     Vector v;
     Vector f;
     double m;
     double t;
     Atom a1;
     Atom a2;
     Atom a3;
     Atom a4;
     std::vector<Atom> atoms;
     double fixed_length_sq;
     Bond b_check;
     Bond b_expect;
     std::vector<Bond> bonds;
     virtual void SetUp() {
         r1.x = 1.0; r1.y = 1.0; r1.z = 1.0;
         r2.x = -1.0; r2.y = -1.0; r2.z = -1.0;
         r3.x = 0.0; r3.y = 0.0; r3.z = 0.0;
         r4.x = 2.0; r4.y = 2.0; r4.z = 2.0;
         v.x = 0.0; v.y = 0.0; v.z = 0.0;
         f.x = 0.0; f.y = 0.0; f.z = 0.0;
         m = 3.0;
         t = 3.0;
         a1 = initialize_atom(m, r1, v, t);
         a2 = initialize_atom(m, r2, v, t);
         a3 = initialize_atom(m, r3, v, t);
         a4 = initialize_atom(m, r4, v, t);
         atoms.push_back(a1);
         atoms.push_back(a2);
         atoms.push_back(a3);
         atoms.push_back(a4);
         fixed_length_sq = 12.0;
     }
  // virtual void TearDown() {}
};

TEST_F(BondTest, BondInitialization) {
    //  initialization with two identical atoms
    EXPECT_THROW(initialize_bond(1, 1, fixed_length_sq), std::invalid_argument);
    //  initialization with proper current_length computation
    b_check = initialize_bond(0, 1, fixed_length_sq);
    std::pair<double, double> current_length_sq = std::make_pair(-1, -1);
    Bond b_expect = {.atom1 = 0, .atom2 = 1,
        .fixed_length_sq = fixed_length_sq,
        .current_length_sq = current_length_sq};
    EXPECT_EQ(b_expect, b_check);
    //  initialization with positions of atoms defined at different times.
    set_time(&(atoms.at(0)), 2*0);
    EXPECT_THROW(update_bond_length(&b_check, atoms), std::invalid_argument);
}

TEST_F(BondTest, CheckBonds) {
    // two pointers are the same
    b_check = initialize_bond(0, 1, fixed_length_sq);
    b_check.atom2 = 0;
    EXPECT_THROW(check_bond(b_check, atoms), std::invalid_argument);
    b_check.atom2 = 1;
    // everything is fine
    EXPECT_NO_THROW(check_bond(b_check, atoms));
    // one of the indices is not in the list
    b_check.atom1 = 10;
    EXPECT_THROW(check_bond(b_check, atoms), std::invalid_argument);
    b_check.atom1 = 0;
    // identical bonds
    b_expect = initialize_bond(0, 1, fixed_length_sq);
    bonds.push_back(b_check);
    bonds.push_back(b_expect);
    EXPECT_THROW(check_bonds(bonds, atoms), std::invalid_argument);
}

TEST_F(BondTest, IOTest) {
    b_check = initialize_bond(0, 1, fixed_length_sq);
    b_expect = initialize_bond(0, 1, fixed_length_sq);
    std::string b_str = bond_to_string(b_check, false);
    EXPECT_EQ(b_expect, string_to_bond(b_str));
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
