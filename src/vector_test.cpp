// 2017 Artur Avkhadiev
/*! \file vector_test.cpp
    \brief tests for vector operations
*/
#include "vector.cpp"
#include <gtest/gtest.h>
TEST(VectorTest, VectorToVectorAddition) {
    Vector v1 = {.x = 1.0, .y = 1.0, .z = 1.0};
    Vector v2 = {.x = -1.0, .y = -1.0, .z = -1.0};
    Vector v = {.x = 0.0, .y = 0.0, .z = 0.0};
    EXPECT_EQ(v, add(v1, v2));
    EXPECT_NE(v, add(v1, v));
}

TEST(VectorTest, VectorToScalarAddition) {
    Vector v1 = {.x = 1.0, .y = 1.0, .z = 1.0};
    double s1 = -1.0;
    double s2 = 1.0;
    Vector v = {.x = 0.0, .y = 0.0, .z = 0.0};
    EXPECT_EQ(v, add(v1, s1));
    EXPECT_NE(v, add(v1, s2));
}

TEST(VectorTest, VectorFromVectorSubtraction) {
    Vector v1 = {.x = 1.0, .y = 1.0, .z = 1.0};
    Vector v2 = {.x = 0.0, .y = 0.0, .z = 0.0};
    EXPECT_EQ(v1, subtract(v1, v2));
    EXPECT_NE(v1, subtract(v2, v1));
}


TEST(VectorTest, DotProduct) {
    Vector v1 = {.x = 1.0, .y = 0.0, .z = 0.0};
    Vector v2 = {.x = 0.0, .y = -1.0, .z = 0.0};
    double s1 = 0.0;
    double s2 = 1.0;
    EXPECT_EQ(s1, dot(v1, v2));
    EXPECT_EQ(s2, dot(v1, v1));
    EXPECT_EQ(s2, dot(v2, v2));
}

TEST(VectorTest, Exponentiaion) {
    Vector v1 = {.x = 2.0, .y = 3.0, .z = 4.0};
    double s = 2.0;
    Vector v = {.x = 4.0, .y = 9.0, .z =  16.0};
    EXPECT_EQ(v, pow(v1, s));
}

TEST(VectorTest, CrossProduct) {
    Vector v1 = {.x = 1.0, .y = 0.0, .z = 0.0};
    Vector v2 = {.x = 0.0, .y = -1.0, .z = 0.0};
    Vector v3 = {.x = 0.0, .y = 0.0, .z = -1.0};
    Vector v4 = {.x = 0.0, .y = 0.0, .z =  1.0};
    Vector v5 = {.x = 0.0, .y = 0.0, .z =  0.0};
    EXPECT_EQ(v3, cross(v1, v2));
    EXPECT_EQ(v4, cross(v2, v1));
    EXPECT_EQ(v5, cross(v1, v1));
    EXPECT_EQ(v5, cross(v2, v2));
}

TEST(VectorTest, NormTest) {
    Vector v1 = {.x = 0.0, .y = 3.0, .z = 4.0};
    double s1 = 25.0;
    double s2 = 5.0;
    EXPECT_EQ(s1, normsq(v1));
    EXPECT_EQ(s2, norm(v1));
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
