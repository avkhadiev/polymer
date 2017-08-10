// 2017 Artur Avkhadiev
/*! \file simple_state_test.cpp
*/
#include <string>
#include <stdexcept>
#include "../include/vector.h"
#include "../include/simple_atom.h"
#include "../include/simple_bond.h"
#include "../include/simple_polymer.h"
#include "../include/simple_state.h"

#include <gtest/gtest.h>

class SimpleStateTest : public ::testing::Test {
 protected:
     // IO settings
     bool verbose;
     bool overwrite;
     int nwriteouts;
     std::string fname;
     std::string outdir;
     std::string indir;
     // Base State Variables
     int nm;
     int nb;
     double m;
     double d;
     Vector zero_vector;
     Vector cm_displacement;
     double t;
     // Atom State Variableis
     Vector atom_displacement, r0, r1, r2, v0, v1, v2;
     std::vector<simple::Atom> atoms;
     std::vector<simple::AtomPolymer> atom_polymers;
     simple::AtomState atom_state;
     // Bond State Variables
     Vector d1, d2, ddot1, ddot2;
     std::vector<simple::Bond> bonds;
     std::vector<simple::BondPolymer> bond_polymers;
     simple::BondState bond_state;
     virtual void SetUp() {
         verbose = true;
         overwrite = true;
         nwriteouts = 100;
         fname ="simple_s_check";
         outdir = "/Users/Arthur/stratt/polymer/test/";
         indir = "/Users/Arthur/stratt/polymer/test/";
         /*
         * STATE:
         *  number of bonds = 2
         *  mass of atoms = 1.0
         *  bond length = 3.0
         * TRIATOMIC MOLECULE:      (1)                     (2)
         *   RCM                    origin            origin + cm_displacement
         *   VCM                    (0.0, 0.0, 0.0)   (0.0, 0.0, 0.0)
         */
         nm = 2;
         nb = 2;
         m = 1.0;
         d = 3.0;
         t = 0.0;
         zero_vector = vector(0.0, 0.0, 0.0);
         cm_displacement = vector(0.0, 0.0, 25.0);
         Vector rcm[nm];
         Vector vcm[nm];
         simple::BaseState::set_nm(nm);
         simple::BondPolymer::set_nb(nb);
         simple::BondPolymer::set_m(m);
         simple::BondPolymer::set_d(d);
         for(int i = 0; i < nm; ++i){
             rcm[i] = add(zero_vector, multiply(cm_displacement, i));
             vcm[i] = zero_vector;
         }
         /* ATOMIC REPRESENTATION IN CM FRAME:
         *  (R = -3Y, V = -1Z) -> (R = 0Y, V = +2Z) -> (R = 3Y, V = -1Z)
         */
         Vector atom_displacement = vector(0.0, 3.0, 0.0);
         Vector r1 = zero_vector;
         Vector r0 = subtract(r1, atom_displacement);
         Vector r2 = add(r1, atom_displacement);
         v0 = vector(0.0, 0.0, -1.0);
         v1 = vector(0.0, 0.0, 2.0);
         v2 = vector(0.0, 0.0, -1.0);
         simple::Atom atom[nm][nb + 1];
         for(int i = 0; i < nm; ++i){
            atom[i][0] = simple::Atom(r0, v0);
            atom[i][1] = simple::Atom(r1, v1);
            atom[i][2] = simple::Atom(r2, v2);
            atoms.push_back(atom[i][0]);
            atoms.push_back(atom[i][1]);
            atoms.push_back(atom[i][2]);
            atom_polymers.push_back(simple::AtomPolymer(atoms, rcm[i], vcm[i]));
            atoms.clear();
         }
         atom_state = simple::AtomState(atom_polymers, t);
         /*
         * BOND REPRESENTATION IN CM FRAME:
         *  BOND1_R = BOND2_R = (0, 1, 0)
         *  BOND1_V = (0, 0, 1), BOND2_V = (0, 0, -1)
         */
         d1 = vector(0.0, 1.0, 0.0);
         d2 = vector(0.0, 1.0, 0.0);
         ddot1 = vector(0, 0, 1.0);
         ddot2 = vector(0, 0, -1.0);
         simple::Bond bond[nm][nb];
         for(int i = 0; i < nm; ++i){
            bond[i][0] = simple::Bond(d1, ddot1);
            bond[i][1] = simple::Bond(d2, ddot2);
            bonds.push_back(bond[i][0]);
            bonds.push_back(bond[i][1]);
            bond_polymers.push_back(simple::BondPolymer(bonds, rcm[i], vcm[i]));
            bonds.clear();
         }
         bond_state = simple::BondState(bond_polymers, t);
     }
  // virtual void TearDown() {}
};

TEST_F(SimpleStateTest, BondStateIO) {
    // prepare storage
    simple::BondState s_expect = bond_state;
    simple::BondState s_check;
    std::string s_str_expect = s_expect.to_string(verbose, overwrite);
    std::string s_str_check;
    // write out state to file
    // write out many times to ensure proper read out later
    for (int i = 0; i < nwriteouts; ++i){
        s_expect.write_to_file(outdir, fname, !verbose, overwrite);
    }
    // read state back from file
    std::vector<simple::BondState> states(nwriteouts);
    simple::read_states_from_file(outdir, fname, states);
    s_check = states.back();
    // check the what has been written vs. read
    EXPECT_EQ(s_check, s_expect);
    //EXPECT_EQ(s_check.to_string(verbose), s_expect.to_string(verbose));
}

TEST_F(SimpleStateTest, AtomStateIO) {
    // prepare storage
    simple::AtomState s_expect = bond_state;
    simple::AtomState s_check;
    std::string s_str_expect = s_expect.to_string(verbose, overwrite);
    std::string s_str_check;
    // write out state to file
    // write out many times to ensure proper read out later
    s_expect.write_to_file(outdir, fname, !verbose, overwrite);
    for (int i = 0; i < nwriteouts - 1; ++i){
        s_expect.write_to_file(outdir, fname, !verbose, !overwrite);
    }
    // read state back from file
    std::vector<simple::AtomState> states(nwriteouts);
    simple::read_states_from_file(outdir, fname, states);
    s_check = states.back();
    // check the what has been written vs. read
    EXPECT_EQ(s_check, s_expect);
    //EXPECT_EQ(s_check.to_string(verbose), s_expect.to_string(verbose));
}

TEST_F(SimpleStateTest, AtomStateToBondState) {
    simple::AtomState s_expect = bond_state;
    simple::AtomState s_check = simple::AtomState(bond_state);
    EXPECT_EQ(s_check, s_expect);
    //printf("%s\n", s_check.to_string().c_str());
    //printf("%s\n", s_expect.to_string().c_str());
}

TEST_F(SimpleStateTest, BondStateToAtomState) {
    simple::AtomState s_expect = atom_state;
    simple::AtomState s_check = simple::BondState(atom_state);
    EXPECT_EQ(s_check, s_expect);
    //printf("%s\n", s_check.to_string().c_str());
    //printf("%s\n", s_expect.to_string().c_str());
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
