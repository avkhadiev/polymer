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
     double epsilon;
     double sigma;
     double calculate_prefactors;
     LJPotential ljpotential_default;
     LJPotential ljpotential_custom;
     std::string sim_name;
     std::string outdir;
     std::string indir;
     virtual void SetUp() {
         r1 = vector(2.0, 2.0, 1.0);
         r2 = vector(1.0, 2.0, 1.0);
         inv_rijsq = 1.0/normsq(subtract(r1, r2));
         epsilon = 2.0;
         sigma = 3.0;
         calculate_prefactors = true;
         ljpotential_default = LJPotential();
         ljpotential_custom = LJPotential(epsilon,
            sigma,
            calculate_prefactors);
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

TEST_F(LJPotentialTest, SetupCustom) {
    EXPECT_EQ(epsilon, ljpotential_custom.get_epsilon());
    EXPECT_EQ(sigma, ljpotential_custom.get_sigma());
    EXPECT_EQ(calculate_prefactors,
        ljpotential_custom.get_calculate_prefactors());
    // test setters
    epsilon += 1.0;
    sigma += 1.0;
    ljpotential_custom.set_epsilon(epsilon);
    ljpotential_custom.set_sigma(sigma);
    ljpotential_custom.set_calculate_prefactors(!calculate_prefactors);
    EXPECT_EQ(epsilon, ljpotential_custom.get_epsilon());
    EXPECT_EQ(sigma, ljpotential_custom.get_sigma());
    EXPECT_EQ(!calculate_prefactors,
        ljpotential_custom.get_calculate_prefactors());
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
    // calculate with no prefactors
    ljpotential_default.set_calculate_prefactors(false);
    neg_pair_virial_expect = 1.0;
    EXPECT_EQ(neg_pair_virial_expect,
        ljpotential_default.calculate_neg_pair_virial(inv_rijsq));
}

TEST_F(LJPotentialTest, CalculateFStrengthOverR) {
    double fstrength_over_r_expect = 24.0;
    EXPECT_EQ(fstrength_over_r_expect,
        ljpotential_default.calculate_fstrength_over_r(inv_rijsq));
    // compare force calculation with -w(rij)/rijsq * rijsq
    // calculate with no prefactors
    ljpotential_default.set_calculate_prefactors(false);
    fstrength_over_r_expect = 1.0;
    EXPECT_EQ(fstrength_over_r_expect,
        ljpotential_default.calculate_fstrength_over_r(inv_rijsq));
}

TEST_F(LJPotentialTest, IOTest) {
    LJPotential ljpotential_check = LJPotential(epsilon, sigma);
    ljpotential_check.writeout_parameters_to_file(outdir, sim_name);
    LJPotential ljpotential_expect = LJPotential();
    ljpotential_expect.read_parameters_from_file(indir, sim_name);
    EXPECT_EQ(ljpotential_expect.get_epsilon(),
        ljpotential_check.get_epsilon());
    EXPECT_EQ(ljpotential_expect.get_sigma(),
            ljpotential_check.get_sigma());
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
