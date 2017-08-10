// 2017 Artur Avkhadiev
/*! \file simple_bond_test.cpp
*/
//#include "atom.cpp"
#include <string>
#include <iostream>
#include <gtest/gtest.h>
#include "../include/vector.h"
#include "../include/simple_atom.h"
#include "../include/simple_bond.h"

class SimpleBondTest : public ::testing::Test {
 protected:
     Vector r1, r2;
     Vector v1, v2;
     simple::Atom atom1, atom2;
     std::string a_str;
     virtual void SetUp() {
         r1 = vector(1.0, 0.0, 0.0);
         r2 = -r1;
         v1 = vector(0.0, 0.0, 1.0);
         v2 = -v1;
         atom1 = simple::Atom(r1, v1);
         atom2 = simple::Atom(r2, v2);
     }
  // virtual void TearDown() {}
};

TEST_F(SimpleBondTest, Initalization) {
    simple::Bond bond_expect = simple::Bond(subtract(r2, r1),
        divide(subtract(v2, v1), norm(subtract(r2, r1))));
    simple::Bond bond_check = simple::Bond(atom1, atom2);
    EXPECT_EQ(bond_check, bond_expect);
}

TEST_F(SimpleBondTest, IO) {
    simple::Bond b_expect = simple::Bond(atom1, atom2);
    std::string b_str = b_expect.to_string(false);      // non-verbose output
    simple::Bond b_check = simple::string_to_bond(b_str);
    EXPECT_EQ(b_check, b_expect);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
