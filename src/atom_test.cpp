// 2017 Artur Avkhadiev
/*! \file atom_test.cpp
*/
//#include "atom.cpp"
#include <string>
#include <utility>              /* std::pair, std::make_pair */
#include <stdexcept>
#include <gtest/gtest.h>
#include "../include/vector.h"
#include "../include/atom.h"

TEST(AtomTest, Initialization) {
    Vector r = {.x = 1.0, .y = 2.0, .z = 3.0};
    Vector v = {.x = 4.0, .y = 4.0, .z = 4.0};
    Vector f = {.x = 0.0, .y = 0.0, .z = 0.0};
    double t = 3;
    std::pair<Vector, double> position = std::make_pair(r, t);
    std::pair<Vector, double> velocity = std::make_pair(v, t);
    std::pair<Vector, double> force = std::make_pair(f, t);
    double m = 3;
    Atom a_expect = {.mass = m,
        .position = position,
        .velocity = velocity,
        .force = force};
    Atom a_check = initialize_atom(m, r, v, t);
    EXPECT_EQ(a_expect, a_check);
    a_expect.mass = m * 2;
    EXPECT_NE(a_expect, a_check);
    EXPECT_EQ(true, is_time_consistent(a_expect));
    EXPECT_EQ(true, is_time_consistent(a_expect, t));
    EXPECT_EQ(false, is_time_consistent(a_expect, t*3));
    a_expect.mass = m;
    a_expect.position.second = 3 * t;
    EXPECT_NE(a_expect, a_check);
    EXPECT_EQ(false, is_time_consistent(a_expect));
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
