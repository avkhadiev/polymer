// 2017 Artur Avkhadiev
/*! \file observable_test.cpp
*/
#include <vector>
#include <string>
#include <map>
#include <gtest/gtest.h>
#include "../include/observable.h"

class ObservableTest : public ::testing::Test {
 protected:
     std::string long_name;
     std::string short_name;
     std::string latex_name;
     std::string units;
     observable::name_t name;
     observable::update_time_t update_time;
     observable::calculate_avg_t calculate_instructions;
     bool print_inst_val;
     bool e_format;
     bool verbose;
     virtual void SetUp() {
        long_name = "Mock Observable";
        short_name = "mock obs";
        latex_name = "\\mathsrc{O}_{\\mathrm{mock}}";
        units = "1";
        name = {.full = long_name,
                .abridged = short_name,
                .latex = latex_name,
                .units = units};
        update_time = observable::MAIN_LOOP;
        calculate_instructions = {.mean = true, .meansq = true};
        print_inst_val = true;
        e_format = true;
        verbose = true;
     }
  // virtual void TearDown() {}
};

TEST_F(ObservableTest, Initilization) {
    Observable mo = Observable(name,
        update_time,
        calculate_instructions,
        print_inst_val,
        e_format);
    EXPECT_EQ(mo.long_name(), long_name);
    EXPECT_EQ(mo.short_name(), parse_string(short_name));
    EXPECT_EQ(mo.tex_name(), latex_name);
    EXPECT_EQ(mo.units(), units);
    EXPECT_EQ(mo.value, 0.0);
    EXPECT_EQ(mo.block_number(), 0);
}
TEST_F(ObservableTest, StringOutput) {
    Observable mo = Observable(name,
        update_time,
        calculate_instructions,
        print_inst_val,
        e_format);
    mo.value = 1.0/3.0;
    fprintf(stdout, "%s\n", mo.to_string(verbose).c_str());
    fprintf(stdout, "%s\n", mo.to_string(!verbose).c_str());
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
