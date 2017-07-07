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

class AtomTest : public ::testing::Test {
 protected:
     Vector r;
     Vector v;
     Vector f;
     double t;
     std::pair<Vector, double> position;
     std::pair<Vector, double> velocity;
     std::pair<Vector, double> force;
     double m;
     Atom a_expect;
     Atom a_check;
     std::string a_str;
     virtual void SetUp() {
         r.x = 1.0; r.y = 2.0; r.z = 3.0;
         v.x = 4.0; v.y = 5.0; v.z = 6.0;
         f.x = 0.0; f.y = 0.0; f.z = 0.0;
         t = 3.0;
         position = std::make_pair(r, t);
         velocity = std::make_pair(v, t);
         force = std::make_pair(f, -1.0);
         m = 3;
         a_expect.mass = m;
         a_expect.position = position;
         a_expect.velocity = velocity;
         a_expect.force = force;
         a_check = initialize_atom(m, r, v, t);
     }
  // virtual void TearDown() {}
};

TEST_F(AtomTest, Initialization) {
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

TEST_F(AtomTest, IOTest) {
    std::string a_str = atom_to_string(a_expect, false);
    a_check = string_to_atom(a_str);
    set_time(&a_check, t);
    EXPECT_EQ(a_expect, a_check);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
