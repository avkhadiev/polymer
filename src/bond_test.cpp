// 2017 Artur Avkhadiev
/*! \file bond_test.cpp
*/
#include <utility>                  /* std::pair, std::make_pair */
#include <cmath>                    /* pow */
#include <iostream>
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
    double fixed_length = 1.0;
    //  initialization with two identical atoms
    EXPECT_THROW(initialize_bond(&a1, &a1, fixed_length), std::invalid_argument);
    //  initialization with proper current_length computation
    Bond b_check = initialize_bond(&a1, &a2, fixed_length);
    std::pair<double, double> current_length_sq = std::make_pair( normsq(subtract(r1, r2)), t);
    Bond b_expect = {.atom1 = &a1, .atom2 = &a2,
        .fixed_length_sq = pow(fixed_length, 2.0),
        .current_length_sq = current_length_sq};
    EXPECT_EQ(b_expect, b_check);
    //  initialization with positions of atoms defined at different times.
    a1.position.second = t * 2;
    ASSERT_EQ(b_check.atom1->position.second, t*2);
    EXPECT_THROW(update_bond_length(&b_check), std::invalid_argument);
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
    std::vector<Atom> atoms_1;
    std::vector<Atom> atoms_2;
    std::vector<Atom> atoms_3;
    atoms_1.push_back(a1); atoms_1.push_back(a2);
    atoms_2.push_back(a1); atoms_2.push_back(a3);
    atoms_3.push_back(a3); atoms_3.push_back(a4);
    double fixed_length = 1.0;
    // two pointers are the same
    Bond b = initialize_bond(&a1, &a2, fixed_length);
    b.atom2 = &a1;
    EXPECT_THROW(check_bond(b, atoms_1), std::invalid_argument);
    b.atom2 = &a2;
    // everything is fine
    EXPECT_NO_THROW(check_bond(b, atoms_1));
    // one of the pointers is not in the list
    EXPECT_THROW(check_bond(b, atoms_2), std::invalid_argument);
    EXPECT_THROW(check_bond(b, atoms_3), std::invalid_argument);
    Bond b2 = initialize_bond(&a3, &a4, fixed_length);
    std::vector<Bond> bonds;
    bonds.push_back(b);
    bonds.push_back(b2);
    EXPECT_THROW(check_bonds(bonds, atoms_3), std::invalid_argument);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
