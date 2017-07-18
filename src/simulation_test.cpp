// 2017 Artur Avkhadiev
/*! \file observable_test.cpp
*/
#include <vector>
#include <string>
#include <map>
#include "../include/state.h"
#include "../include/observable_container.h"
#include "../include/ljpotential.h"
#include "../include/force_updater.h"
#include "../include/integrator.h"
#include "../include/verlet_integrator.h"
#include "../include/simulation.h"
#include <gtest/gtest.h>
class SimulationTest : public ::testing::Test {
protected:
    // mock scalar / vector observable is already defined
    // and put in test container
    // sim_name is already defined
    // outdir and indir are specified
    ScalarObservable mock_scalar;
    VectorObservable mock_vector;
    ObservableContainer test_container;
    LJPotential potential;
    std::string sim_name;
    double timestep;
    double measurestep;
    virtual void SetUp() {
        // set up observable container
        mock_scalar = declare_scalar_observable("Mock Scalar",
         "scalar unit",
         "MS");
        mock_vector = declare_vector_observable("Mock Vector",
         "vector unit",
         "MV");
        test_container = ObservableContainer();
        test_container.add_scalar_observable(&mock_scalar);
        test_container.add_vector_observable(&mock_vector);
        // set up the potential
        potential = LJPotential();
        sim_name = "test_sim";
        timestep = 3.0;
        measurestep = 2.0;
    }
    // virtual void TearDown() {}
};
TEST_F(SimulationTest, Initialization){
    // set up the simulation
    // set up the integrator, give pointer to store kinetic energy
    VerletIntegrator verlet_integrator
        = VerletIntegrator(ForceUpdater(LJPotential()), test_container.get_scalar_observable_accumulator("Mock Scalar"));
    Integrator& integrator = verlet_integrator;
    ObservableContainer& container = test_container;
    Simulation test_sim = Simulation(sim_name,
        integrator,
        container);
    test_sim.set_timestep(timestep);
    test_sim.set_measurestep(measurestep);
    // test basic getters
    EXPECT_EQ(sim_name, test_sim.get_name());
    EXPECT_EQ(timestep, test_sim.get_timestep());
    EXPECT_EQ(measurestep, test_sim.get_measurestep());
    // get list of observables
    ObservableContainer& container_check =
        test_sim.get_observables();
    ScalarObservable *ms = container_check.get_scalar_observable("Mock Scalar");
    EXPECT_EQ("Mock Scalar", ms->name);
    EXPECT_EQ("scalar unit", ms->units);
    EXPECT_EQ("MS", ms->axis_name);
}
int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
