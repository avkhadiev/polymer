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

TEST(LJPotentialTest, SetupDefault) {
    LJPotential ljpotential_default = LJPotential();
    EXPECT_EQ(1.0, ljpotential_default.get_epsilon());
    EXPECT_EQ(1.0, ljpotential_default.get_sigma());
}

TEST(LJPotentialTest, SetupCustom) {
    double epsilon = 2.0;
    double sigma = 3.0;
    LJPotential ljpotential_custom = LJPotential( epsilon, sigma );
    EXPECT_EQ(epsilon, ljpotential_custom.get_epsilon());
    EXPECT_EQ(sigma, ljpotential_custom.get_sigma());
}

TEST(LJPotentialTest, CalculateForce) {
    LJPotential ljpotential_default = LJPotential();
    Vector r1 = {.x = 2.0, .y = 2.0, .z = 1.0};
    Vector r2 = {.x = 1.0, .y = 2.0, .z = 1.0};
    double force_strength = 24.0;
    Vector r12 = subtract(r1, r2);
    Vector r21 = subtract(r2, r1);
    EXPECT_EQ(multiply(r12, force_strength), ljpotential_default.calculate_force(r1, r2));
    EXPECT_EQ(multiply(r21, force_strength), ljpotential_default.calculate_force(r2, r1));
    EXPECT_THROW(ljpotential_default.calculate_force(r1, r1), std::invalid_argument);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
