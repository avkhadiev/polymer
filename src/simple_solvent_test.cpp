// 2017 Artur Avkhadiev
/*! \file simple_polymer_test.cpp
*/
#include <string>
#include <iostream>
#include <gtest/gtest.h>
#include "../include/vector.h"
#include "../include/simple_atom.h"
#include "../include/simple_solvent.h"
class SimpleSolventTest : public ::testing::Test {
 protected:
     Vector r, v;
     simple::Atom atom;
     simple::Solvent molecule;
     virtual void SetUp() {
         r = v = vector(0.0, 0.0, 0.0);
         atom = simple::Atom(r, v);
         molecule = simple::Solvent(r, v);
     }
  // virtual void TearDown() {}
};

TEST_F(SimpleSolventTest, Initalization) {
    EXPECT_EQ(molecule.atom, atom);
}

TEST_F(SimpleSolventTest, IO) {
    bool verbose = false;
    EXPECT_EQ(molecule.to_string(verbose), atom.to_string(verbose) + "\n");
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
