// 2017 Artur Avkhadiev
/*! \file simple_simulation_test.cpp
*/
#include <vector>
#include <string>
#include <map>
#include "../include/state.h"
#include "../include/observable_container.h"
#include "../include/ljpotential.h"
#include "../include/force_updater.h"
#include "../include/integrator.h"
#include "../include/rattle_integrator.h"
#include "../include/simple_simulation.h"
#include <gtest/gtest.h>

class SimpleSimulationTest : public ::testing::Test {
protected:
    std::string sim_name;
    std::string cndir;
    std::string tpdir;
    std::string dtdir;
    double dt;
    size_t icalc;
    size_t iprint;
    size_t isave;
    size_t idata;
    size_t itape;
    virtual void SetUp() {
        sim_name = "test_sim";
        cndir = tpdir = dtdir =  "/Users/Arthur/stratt/polymer/test/";
        dt = 0.01;
        icalc = 1;
        iprint = 0;
        isave = 0;
        idata = 0;
        itape = 0;
    }
    // virtual void TearDown() {}
};

TEST_F(SimpleSimulationTest, Initialization){
    simple::Simulation sim = simple::Simulation(sim_name, cndir, tpdir, dtdir,
                                dt, icalc, iprint, isave, idata, itape);
    sim.evolve(0.1);           /**> time in LJ units*/
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
