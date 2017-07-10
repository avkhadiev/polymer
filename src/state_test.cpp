// 2017 Artur Avkhadiev
/*! \file state_test.cpp
*/
#include <string>
#include <stdexcept>
#include <gtest/gtest.h>
#include "../include/vector.h"
#include "../include/atom.h"
#include "../include/bond.h"
#include "../include/molecule.h"
#include "../include/state.h"

class StateTest : public ::testing::Test {
 protected:
     static const int na = 3;
     static const int nb = 2;
     static const int nm = 2;
     double t;
     Vector r[na];
     Vector v[na];
     double mass[na];
     double fixed_length_sq;
     Atom a[na];
     std::vector<Atom> atoms;
     Bond b[nb];
     std::vector<Bond> bonds;
     Molecule m[nm];
     std::vector<Molecule> molecules;
     virtual void SetUp() {
         t = 3.0;
         // initialize atoms
         // (-1, -1, -1), (0, 0, 0), and (1, 1, 1)
         for (int i = 0; i < na;  ++i) {
             double coord = i - 1.0;
             mass[i] = 1.0;
             r[i] = vector(coord, coord, coord);
             v[i] = vector(coord, coord, coord);
             a[i] = initialize_atom(mass[i], r[i], v[i], t);
             atoms.push_back(a[i]);
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
         for (int i = 0; i < nm; ++i) {
             try
             {
                 m[i] = initialize_molecule(atoms, bonds, t);
             }
             catch (std::invalid_argument &e)
             {
                 throw;
             }
             molecules.push_back(m[i]);
         }
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
    std::cout << "**************************" << std::endl;
    std::cout << "Verbose output of state:" << std::endl;
    std::cout << "**************************" << std::endl;
    std::cout << state_to_string(s_check, true) << std::endl;
    std::cout << "**************************" << std::endl;
    std::cout << "Non-verbose output of state:" << std::endl;
    std::cout << "**************************" << std::endl;
    std::cout << state_to_string(s_check, false) << std::endl;
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
