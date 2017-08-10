// 2017 Artur Avkhadiev
/*! \file ljpotential_observable_container_test.cpp
*/
#include <string>
#include <utility>              /* std::pair, std::make_pair */
#include <vector>
#include <stdexcept>
#include <gtest/gtest.h>
#include "../include/observable_container.h"

class ObservableContainerTest : public ::testing::Test {
 protected:
    ScalarObservable mock_scalar;
    VectorObservable mock_vector;
    std::string sim_name;
    std::string outdir;
    std::string indir;
    bool overwrite;
    ObservableContainer test_container;
    virtual void SetUp() {
        mock_scalar = declare_scalar_observable("Mock Scalar",
         "scalar unit",
         "MS");
        mock_vector = declare_vector_observable("Mock Vector",
         "vector unit",
         "MV");
        sim_name = "testsim";
        outdir = "/Users/Arthur/stratt/polymer/test/";
        indir = outdir;
        overwrite = true;
        test_container = ObservableContainer();
        test_container.add_scalar(mock_scalar);
        test_container.add_vector(mock_vector);
     }
  // virtual void TearDown() {}
};

TEST_F(ObservableContainerTest, Initialization){
    // test scalars
    ScalarObservable& ms =
        test_container.get_scalar(mock_scalar.name);
    EXPECT_EQ(mock_scalar.name, ms.name);
    EXPECT_EQ(mock_scalar.units, ms.units);
    EXPECT_EQ(mock_scalar.axis_name, ms.axis_name);
    EXPECT_THROW(test_container.get_scalar(mock_vector.name), std::invalid_argument);
    // test vectors
    VectorObservable mv =
        test_container.get_vector(mock_vector.name);
    EXPECT_EQ(mock_vector.name, mv.name);
    EXPECT_EQ(mock_vector.units, mv.units);
    EXPECT_EQ(mock_vector.axis_name, mv.axis_name);
    EXPECT_THROW(test_container.get_vector(mock_scalar.name), std::invalid_argument);
}

TEST_F(ObservableContainerTest, ScalarManipulation){
    double *s_acc = &(test_container.get_scalar(mock_scalar.name).accumulator);
    // recording observables via accumulators
    ScalarObservable& ms =
        test_container.get_scalar(mock_scalar.name);
    EXPECT_EQ(true, ms.value_time.empty());
    *s_acc = 3.0;
    double t = 1.0;
    test_container.update(mock_scalar.name, t);
    EXPECT_EQ(false, ms.value_time.empty());
    EXPECT_EQ(*s_acc, ms.value_time.at(0).first);
    EXPECT_EQ(t, ms.value_time.at(0).second);
}

TEST_F(ObservableContainerTest, VectorManipulation){
    Vector *v_acc = &(test_container.get_vector(mock_vector.name).accumulator);
    VectorObservable &mv =
        test_container.get_vector(mock_vector.name);
    EXPECT_EQ(true, mv.value_time.empty());
    *v_acc = vector(3.0, 3.0, 3.0);
    double t = 1.0;
    test_container.update(mock_vector.name, t);
    EXPECT_EQ(false, mv.value_time.empty());
    EXPECT_EQ(*v_acc, mv.value_time.at(0).first);
    EXPECT_EQ(t, mv.value_time.at(0).second);
}

TEST_F(ObservableContainerTest, Writeout){
    double *s_acc = &(test_container.get_scalar(mock_scalar.name).accumulator);
    Vector *v_acc = &(test_container.get_vector(mock_vector.name).accumulator);
    double t = 0.0;
    for (int i = 0; i < 5; ++i){
        t += 10.0;
        *s_acc += 1.0;
        *v_acc = add(*v_acc, vector(1.0, 2.0, 3.0));
        test_container.update(t);
    }
    test_container.writeout(outdir, sim_name, overwrite);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
