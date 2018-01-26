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
        std::string sim_name, outdir;
    // long names, short names, tex names, and units
    std::string ln1, ln2, sn1, sn2, tn1, tn2, units;
    observable::name_t name1, name2;
    observable::update_time_t update_time;
    observable::calculate_avg_t calculate_instructions;
    bool print_inst_val;
    bool e_format;
    bool verbose;
    virtual void SetUp() {
        sim_name = "test_sim";
        outdir =  "/Users/Arthur/stratt/polymer/test/";
        ln1 = "Mock Observable 1";
        ln2 = "Mock Observable 2";
        sn1 = "mock obs 1";
        sn2 = "mock obs 2";
        tn1 = "\\mathsrc{O}_{\\mathrm{mock} 1}";
        tn2 = "\\mathsrc{O}_{\\mathrm{mock} 2}";
        units = "1";
        name1 = {.full = ln1,
                .abridged = sn1,
                .latex = tn1,
                .units = units};
        name2 = {.full = ln2,
                .abridged = sn2,
                .latex = tn2,
                .units = units};
        update_time = observable::MAIN_LOOP;
        calculate_instructions = {.mean = true, .meansq = true};
        print_inst_val = true;
        e_format = true;
        verbose = true;
     }
  // virtual void TearDown() {}
};

TEST_F(ObservableContainerTest, Initialization){
    Observable mo1 = Observable(name1,
        update_time,
        calculate_instructions,
        print_inst_val,
        e_format);
    Observable mo2 = Observable(name2,
        update_time,
        calculate_instructions,
        print_inst_val,
        e_format);
    mo1.value = 1.0 / 3.0;
    mo2.value = 2.0 / 3.0;
    std::vector<Observable*> observables = {&mo1, &mo2};
    ObservableContainer test_container = ObservableContainer(observables);
    EXPECT_EQ(test_container.nobservables(), observables.size());
    std::cout << test_container.status_string(verbose);
}

TEST_F(ObservableContainerTest, Writeout){
    Observable mo1 = Observable(name1,
        update_time,
        calculate_instructions,
        print_inst_val,
        e_format);
    Observable mo2 = Observable(name2,
        update_time,
        calculate_instructions,
        print_inst_val,
        e_format);
    mo1.value = 1.0 / 3.0;
    mo2.value = 2.0 / 3.0;
    std::vector<Observable*> observables = {&mo1, &mo2};
    ObservableContainer test_container = ObservableContainer(observables);
    std::cout << "Expect values 1/3, 2/3" << std::endl;
    std::cout << test_container.status_string(!verbose);
    test_container.prepare_datafiles(outdir, sim_name);
    test_container.write_data(outdir, sim_name);
    mo1.value = 2.0 / 3.0;
    mo2.value = 1.0 / 3.0;
    std::cout << "Changed values to 2/3, 1/3" << std::endl;
    std::cout << test_container.status_string(!verbose);
    test_container.record_observables();
    test_container.write_data(outdir, sim_name);
    std::cout << "Observales recorded, so values should be now 0" << std::endl;
    std::cout << test_container.status_string(!verbose);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
