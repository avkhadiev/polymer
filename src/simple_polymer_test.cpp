// 2017 Artur Avkhadiev
/*! \file simple_polymer_test.cpp
*/
#include <string>
#include <iostream>
#include <gtest/gtest.h>

#include "../include/vector.h"
#include "../include/simple_atom.h"
#include "../include/simple_bond.h"
#include "../include/simple_polymer.h"

class SimplePolymerTest : public ::testing::Test {
 protected:
     // BASE
     Vector rcm, vcm;
     // ATOM REPRESENTATION
     Vector r0, r1, r2;
     Vector v0, v1, v2;
     simple::Atom atom0, atom1, atom2;
     std::vector<simple::Atom> atoms;
     simple::AtomPolymer atom_polymer;
     // BOND REPRESENTATION
     Vector d1, d2;
     Vector ddot1, ddot2;
     simple::Bond bond1, bond2;
     std::vector<simple::Bond> bonds;
     simple::BondPolymer bond_polymer;
     virtual void SetUp() {
         /*
         * TRIATOMIC MOLECULE:
         *  number of bonds = 2
         *  mass of atoms = 1.0
         *  bond length = 3.0
         */
         simple::BondPolymer::set_nb(2);
         simple::BondPolymer::set_m(1.0);
         simple::BondPolymer::set_d(3.0);
         rcm = vcm = vector(0.0, 0.0, 0.0);
         /*
         *  RCM at origin, VCM zero
         * ATOMIC REPRESENTATION:
         *  (R = -3Y, V = -1Z) -> (R = 0Y, V = +2Z) -> (R = 3Y, V = -1Z)
         */
         r1 = vector(0.0, 0.0, 0.0);
         r0 = subtract(r1, vector(0.0, 3.0, 0.0));
         r2 = add(r1, vector(0.0, 3.0, 0.0));
         v0 = vector(0.0, 0.0, -1.0);
         v1 = vector(0.0, 0.0, 2.0);
         v2 = vector(0.0, 0.0, -1.0);
         // make sure bond velocities are 0 when projected onto the bond
         // atom polymer initialization does not check for that --- it's not
         // supposed to (in the simulation,
         //     initialization is done with bond polymers).
         atom0 = simple::Atom(r0, v0);
         atoms.push_back(atom0);
         atom1 = simple::Atom(r1, v1);
         atoms.push_back(atom1);
         atom2 = simple::Atom(r2, v2);
         atoms.push_back(atom2);
         atom_polymer = simple::AtomPolymer(atoms);
         /*
         * BOND REPRESENTATION:
         *  BOND1_R = BOND2_R = (0, 1, 0)
         *  BOND1_V = (0, 0, 1), BOND2_V = (0, 0, -1)
         */
         d1 = d2 = vector(0.0, 1.0, 0.0);
         ddot1 = vector(0, 0, 1.0);
         ddot2 = vector(0, 0, -1.0);
         bond1 = simple::Bond(d1, ddot1);
         bonds.push_back(bond1);
         bond2 = simple::Bond(d2, ddot2);
         bonds.push_back(bond2);
         bond_polymer = simple::BondPolymer(bonds, rcm, vcm);
     }
  // virtual void TearDown() {}
};

TEST_F(SimplePolymerTest, Initalization) {
    simple::AtomPolymer atom_polymer_check;
    simple::BondPolymer bond_polymer_check;
    simple::AtomPolymer atom_polymer_expect = atom_polymer;
    simple::BondPolymer bond_polymer_expect = bond_polymer;
    // add overall drift; it should be subtracted during initialization
    Vector v_extra = vector(60.0, 0, 0);
    v0 += v_extra;
    v1 += v_extra;
    v2 += v_extra;
    atoms.clear();
    atom0 = simple::Atom(r0, v0);
    atoms.push_back(atom0);
    atom1 = simple::Atom(r1, v1);
    atoms.push_back(atom1);
    atom2 = simple::Atom(r2, v2);
    atoms.push_back(atom2);
    atom_polymer_check = simple::AtomPolymer(atoms);
    // add small projection along bond; it should be subtracted during
    // initialization
    Vector ddot1_proj = multiply(d1, 0.01);
    Vector ddot2_proj = multiply(d2, 0.01);
    ddot1 += ddot1_proj;
    ddot2 += ddot2_proj;
    bonds.clear();
    bond1 = simple::Bond(d1, ddot1);
    bonds.push_back(bond1);
    bond2 = simple::Bond(d2, ddot2);
    bonds.push_back(bond2);
    bond_polymer_check = simple::BondPolymer(bonds, rcm, vcm);
    // check equivalence
    EXPECT_EQ(atom_polymer_check.to_string(true),
        atom_polymer_expect.to_string(true));
    EXPECT_EQ(bond_polymer_check.to_string(true),
        bond_polymer_expect.to_string(true));
    printf("%s\n", atom_polymer_check.to_string(true).c_str());
    printf("%s\n", bond_polymer_check.to_string(true).c_str());
}

TEST_F(SimplePolymerTest, BondToAtomConversion) {
    EXPECT_EQ(atom_polymer, simple::AtomPolymer(bond_polymer));
}

TEST_F(SimplePolymerTest, AtomToBondConversion) {
    EXPECT_EQ(bond_polymer, simple::BondPolymer(atom_polymer));
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
