// 2017 Artur Avkhadiev
/*! \file ljpotential_test.cpp
*/
//#include "atom.cpp"
#include <string>
#include <utility>              /* std::pair, std::make_pair */
#include <stdexcept>
#include <gtest/gtest.h>
#include "../include/vector.h"
#include "../include/ljpotential.h"

class LJPotentialTest : public ::testing::Test {
 protected:
     Vector r1;
     Vector r2;
     double inv_rijsq;
     LJPotential ljpotential_default;
     std::string sim_name;
     std::string outdir;
     std::string indir;
     virtual void SetUp() {
         r1 = vector(2.0, 2.0, 1.0);
         r2 = vector(1.0, 2.0, 1.0);
         inv_rijsq = 1.0/normsq(subtract(r1, r2));
         ljpotential_default = LJPotential();
        sim_name = "testsim";
        outdir = "/Users/Arthur/stratt/polymer/test/";
        indir = outdir;
     }
  // virtual void TearDown() {}
};

TEST_F(LJPotentialTest, SetupDefault) {
    EXPECT_EQ(1.0, ljpotential_default.get_epsilon());
    EXPECT_EQ(1.0, ljpotential_default.get_sigma());
}

TEST_F(LJPotentialTest, CalculatePairPotential) {
    double pair_potential_expect = 0.0;
    EXPECT_EQ(pair_potential_expect,
        ljpotential_default.calculate_pair_potential(inv_rijsq));
}

TEST_F(LJPotentialTest, CalculatePairVirial) {
    double neg_pair_virial_expect = 24.0;
    EXPECT_EQ(neg_pair_virial_expect,
        ljpotential_default.calculate_neg_pair_virial(inv_rijsq));
}

TEST_F(LJPotentialTest, CalculateFStrengthOverR) {
    double fstrength_over_r_expect = 24.0;
    EXPECT_EQ(fstrength_over_r_expect,
        ljpotential_default.calculate_fstrength_over_r(inv_rijsq));
    // compare force calculation with -w(rij)/rijsq * rijsq
}
int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
