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

// TODO separate into test cases, add test fixtures

TEST(BondTest, BondInitialization) {
    Vector r1 = {.x = 1.0, .y = 1.0, .z = 1.0};
    Vector r2 = {.x = -1.0, .y = -1.0, .z = -1.0};
    Vector v = {.x = 0.0, .y = 0.0, .z = 0.0};
    Vector f = {.x = 0.0, .y = 0.0, .z = 0.0};
    double m = 1.0;
    double t = 3.0;
    Atom a1 = initialize_atom(m, r1, v, t);
    Atom a2 = initialize_atom(m, r2, v, t);
    std::vector<Atom> atoms;
    atoms.push_back(a1);
    atoms.push_back(a2);
    double fixed_length = 1.0;
    //  initialization with two identical atoms
    EXPECT_THROW(initialize_bond(1, 1, fixed_length), std::invalid_argument);
    //  initialization with proper current_length computation
    Bond b_check = initialize_bond(0, 1, fixed_length);
    std::pair<double, double> current_length_sq = std::make_pair(-1, -1);
    Bond b_expect = {.atom1 = 0, .atom2 = 1,
        .fixed_length_sq = pow(fixed_length, 2.0),
        .current_length_sq = current_length_sq};
    EXPECT_EQ(b_expect, b_check);
    //  initialization with positions of atoms defined at different times.
    atoms.at(0).position.second = t * 2;
    ASSERT_EQ(atoms.at(0).position.second, t*2);
    EXPECT_THROW(update_bond_length(&b_check, atoms), std::invalid_argument);
}
TEST(BondTest, CheckBonds) {
    Vector r1 = {.x = 1.0, .y = 1.0, .z = 1.0};
    Vector r2 = {.x = -1.0, .y = -1.0, .z = -1.0};
    Vector r3 = {.x = 0.0, .y = 0.0, .z = 0.0};
    Vector r4 = {.x = 2.0, .y = 2.0, .z = 2.0};
    Vector v = {.x = 0.0, .y = 0.0, .z = 0.0};
    Vector f = {.x = 0.0, .y = 0.0, .z = 0.0};
    double m = 1.0;
    double t = 3.0;
    Atom a1 = initialize_atom(m, r1, v, t);
    Atom a2 = initialize_atom(m, r2, v, t);
    Atom a3 = initialize_atom(m, r3, v, t);
    Atom a4 = initialize_atom(m, r4, v, t);
    std::vector<Atom> atoms;
    atoms.push_back(a1);
    atoms.push_back(a2);
    atoms.push_back(a3);
    atoms.push_back(a4);
    double fixed_length = 1.0;
    // two pointers are the same
    Bond b = initialize_bond(0, 1, fixed_length);
    b.atom2 = 0;
    EXPECT_THROW(check_bond(b, atoms), std::invalid_argument);
    b.atom2 = 1;
    // everything is fine
    EXPECT_NO_THROW(check_bond(b, atoms));
    // one of the indices is not in the list
    b.atom1 = 10;
    EXPECT_THROW(check_bond(b, atoms), std::invalid_argument);
    b.atom1 = 0;
    Bond b2 = initialize_bond(0, 1, fixed_length);
    std::vector<Bond> bonds;
    bonds.push_back(b);
    bonds.push_back(b2);
    // identical bonds
    EXPECT_THROW(check_bonds(bonds, atoms), std::invalid_argument);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
