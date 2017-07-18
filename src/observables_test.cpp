// 2017 Artur Avkhadiev
/*! \file observable_test.cpp
*/
#include <vector>
#include <string>
#include <map>
#include <gtest/gtest.h>
#include "../include/state.h"
#include "../include/parsing.h"
#include "../include/vector.h"
#include "../include/observable_container.h"

class ObservableTest : public ::testing::Test {
 protected:
     std::string scalar_name;
     std::string vector_name;
     std::string scalar_axis_name;
     std::string vector_axis_name;
     std::string units;
     std::string so_string;
     std::string vo_string;
     double sacc;
     Vector vacc;
     double mock_time;
     double mock_scalar_value;
     Vector mock_vector_value;
     std::vector<std::pair<double, double> > value_time_scalar;
     std::vector<std::pair<Vector, double> > value_time_vector;
     ScalarObservable so;
     VectorObservable vo;
     std::string sim_name;
     std::string outdir;
     std::string indir;
     bool overwrite;
     virtual void SetUp() {
         vector_name = "mock vector";
         vector_axis_name = vector_name;
         scalar_name = "mock scalar";
         scalar_axis_name = scalar_name;
         units = "mock units";
         sacc = 0.0;
         vacc = vector(0.0, 0.0, 0.0);
         mock_time = 1.0;
         mock_scalar_value = 2.0;
         mock_vector_value = vector(3.0, 4.0, 5.0);
         value_time_scalar.push_back(std::pair<double, double>(mock_scalar_value, mock_time));
         value_time_vector.push_back(std::pair<Vector, double>(mock_vector_value, mock_time));
         so = declare_scalar_observable(scalar_name, units);
         vo = declare_vector_observable(vector_name, units);
         so_string = scalar_axis_name + ", " + units;
         vo_string = vector_axis_name + ", " + units;
         sim_name = "testsim";
         outdir = "/Users/Arthur/stratt/polymer/test/";
         indir = outdir;
         overwrite = true;
     }
  // virtual void TearDown() {}
};

TEST_F(ObservableTest, ScalarInitilization) {
    EXPECT_EQ(scalar_name, so.name);
    EXPECT_EQ(units, so.units);
    EXPECT_EQ(scalar_axis_name, so.axis_name);
    EXPECT_EQ(sacc, so.accumulator);
    EXPECT_EQ(true, so.value_time.empty());
}
TEST_F(ObservableTest, VectorInitilization) {
    EXPECT_EQ(vector_name, vo.name);
    EXPECT_EQ(units, vo.units);
    EXPECT_EQ(vector_axis_name, vo.axis_name);
    EXPECT_EQ(vacc, vo.accumulator);
    EXPECT_EQ(true, vo.value_time.empty());
}
TEST_F(ObservableTest, StringOutput) {
    EXPECT_EQ(so_string, scalar_observable_to_string(so));
    EXPECT_EQ(vo_string, vector_observable_to_string(vo));
}

TEST_F(ObservableTest, ClearTest){
    // add stuff to the observable value_time vector
    so.value_time.push_back(std::make_pair(2.0, 2.0));
    vo.value_time.push_back(std::make_pair(vector(2.0, 2.0, 2.0), 2.0));
    ASSERT_EQ(false, so.value_time.empty());
    ASSERT_EQ(false, vo.value_time.empty());
    // clear value_time vectors in observables
    clear_observable_records(&so);
    clear_observable_records(&vo);
    // check whether the records are emptied
    ASSERT_EQ(true, so.value_time.empty());
    ASSERT_EQ(true, vo.value_time.empty());
}

TEST_F(ObservableTest, IOTest) {
    so.value_time = value_time_scalar;
    vo.value_time = value_time_vector;
    write_scalar_observable_to_file(so, outdir, sim_name, overwrite);
    write_vector_observable_to_file(vo, outdir, sim_name, overwrite);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
