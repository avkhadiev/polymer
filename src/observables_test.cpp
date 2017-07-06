// 2017 Artur Avkhadiev
/*! \file observable_test.cpp
*/
#include <map>
#include <utility>      /* std::pair, std::make_pair */
#include <vector>
#include <string>
#include <gtest/gtest.h>
#include "../include/vector.h"
#include "../include/observables.h"

class ObservableTest : public ::testing::Test {
 protected:
     std::string name;
     std::string units;
     std::string observable_string;
     ScalarObservable so;
     VectorObservable vo;
     virtual void SetUp() {
         name = "mock name";
         units = "mock units";
         so = declare_scalar_observable(name, units);
         vo = declare_vector_observable(name, units);
         observable_string = name + ", " + units;
     }
  // virtual void TearDown() {}
};

TEST_F(ObservableTest, ScalarInitilization) {
    EXPECT_EQ(name, so.name);
    EXPECT_EQ(units, so.units);
    EXPECT_EQ(true, so.value_time.empty());
}
TEST_F(ObservableTest, VectorInitilization) {
    EXPECT_EQ(name, vo.name);
    EXPECT_EQ(units, vo.units);
    EXPECT_EQ(true, vo.value_time.empty());
}
TEST_F(ObservableTest, StringOutput) {
    EXPECT_EQ(observable_string, scalar_observable_to_string(so));
    EXPECT_EQ(observable_string, vector_observable_to_string(vo));
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
