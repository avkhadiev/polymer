// 2017 Artur Avkhadiev
/*! \file simple_atom_test.cpp
*/
//#include "atom.cpp"
#include <string>
#include <iostream>
#include <gtest/gtest.h>
#include "../include/vector.h"
#include "../include/simple_atom.h"

class SimpleAtomTest : public ::testing::Test {
 protected:
     Vector r;
     Vector v;
     Vector f;
     simple::Atom a_expect;
     simple::Atom a_check;
     std::string a_str;
     virtual void SetUp() {
         r = vector(1.0, 2.0, 3.0);
         v = vector(4.0, 5.0, 6.0);
         f = vector(0.0, 0.0, 0.0);
         a_expect = simple::Atom(r, v, f);
         a_check = a_expect;
     }
  // virtual void TearDown() {}
};

TEST_F(SimpleAtomTest, IOTest) {
    std::string a_str = a_expect.to_string(false);      // non-verbose output
    a_check = simple::string_to_atom(a_str);
    EXPECT_EQ(a_check, a_expect);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
