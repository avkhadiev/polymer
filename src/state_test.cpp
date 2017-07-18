// 2017 Artur Avkhadiev
/*! \file state_test.cpp
*/
#include <string>
#include <stdexcept>
#include "../include/vector.h"
#include "../include/atom.h"
#include "../include/bond.h"
#include "../include/molecule.h"
#include "../include/state.h"

#include <gtest/gtest.h>

class StateTest : public ::testing::Test {
 protected:
     static const int na = 3;
     static const int nb = 2;
     static const int nm = 2;
     double t;
     Vector r[nm][na];
     Vector v[nm][na];
     double mass[nm][na];
     double fixed_length_sq;
     Atom a[nm][na];
     std::vector<Atom> atoms;
     Bond b[nb];
     std::vector<Bond> bonds;
     Molecule m[nm];
     std::vector<Molecule> molecules;
     bool verbose;
     bool overwrite;
     std::string fname;
     std::string outdir;
     std::string indir;
     std::vector<State> states;
     virtual void SetUp() {
         t = 3.0;
         // initialize atoms
         // (-1, -1, -1), (0, 0, 0), and (1, 1, 1)
         for (int j = 0; j < nm; ++j) {
             for (int i = 0; i < na;  ++i) {
                 double coord = i - 1.0 + 5.0 * j;
                 mass[j][i] = 1.0 + j;
                 r[j][i] = vector(coord, coord, coord);
                 v[j][i] = vector(coord, coord, coord);
                 a[j][i] = initialize_atom(mass[j][i], r[j][i], v[j][i], t);
                 atoms.push_back(a[j][i]);
             }
         }
         // initialize bonds
         fixed_length_sq = 3.0;
         for (int i = 0; i < nb; ++i) {
             try
             {
                 b[i] = initialize_bond(i, i + 1, fixed_length_sq);
             }
             catch (std::invalid_argument &e)
             {
                 throw;
             }
             bonds.push_back(b[i]);
         }
         // initialize molecules
         for (int j = 0; j < nm; ++j) {
             try
             {
                 std::vector<Atom>::const_iterator first = atoms.begin() + j * na;
                 std::vector<Atom>::const_iterator last = atoms.begin() + j*na + na;
                 std::vector<Atom> atoms_for_molecule (first, last);
                 m[j] = initialize_molecule(atoms_for_molecule, bonds, t);
             }
             catch (std::invalid_argument &e)
             {
                 throw;
             }
             molecules.push_back(m[j]);
         }
         // set up arguments for I/O;
         verbose = true;
         overwrite = true;
         fname = "s_check";
         outdir = "/Users/Arthur/stratt/polymer/test/";
         indir = outdir;
     }
  // virtual void TearDown() {}
};

TEST_F(StateTest, CheckState) {
    ASSERT_NO_THROW(initialize_state(molecules, t));
    State s = initialize_state(molecules, t);
    s.molecules.pop_back();
    // number of molecules is not equal to the one declared
    EXPECT_THROW(check_state(s), std::invalid_argument);
}

TEST_F(StateTest, SetTime) {
    State s_check = initialize_state(molecules, t);
    set_time(&s_check, 2*t);
    EXPECT_EQ(true, is_time_consistent(s_check));
    // can't simply change s.time: molecules and atoms also contain time records
    s_check.time = t;
    EXPECT_EQ(false, is_time_consistent(s_check));
}


TEST_F(StateTest, IOTest) {
    State s_check = initialize_state(molecules, t);
    State s_expect = s_check;
    // write out s_check to file twice
    write_state_to_file(s_check, outdir, fname, !verbose, overwrite);
    // read s_check from file
    try
    {
        states = read_states_from_file(indir, fname);
    }
    catch (std::invalid_argument &e)
    {
        throw;
    }
    // extract the state read
    s_check = states.back();
    // check the what has been written vs. read
    std::string s_str_check = state_to_string(s_check, false);
    std::string s_str_expect = state_to_string(s_expect, false);
    EXPECT_EQ(s_str_expect, s_str_check);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
