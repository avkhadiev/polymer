// 2017 Artur Avkhadiev
/*! \file lj_verlet_test.h
*/
#ifndef LJ_VERLET_TEST_H
#define LJ_VERLET_TEST_H
#include <map>
#include "observables.h"
#include "observable_container.h"
#include "simulation.h"
extern int ncycles;
extern double timestep;
extern double measurestep;
extern const std::string default_outdir;
// number of observations to store in memory before writeout
extern const int observations_before_writeout;
// observable names
extern const std::string kinetic_energy_str;
extern const std::string potential_energy_str;
extern const std::string energy_str;
extern const std::string neg_virial_str;
// observable axes labels
extern const std::string kinetic_energy_axis_str;
extern const std::string potential_energy_axis_str;
extern const std::string energy_axis_str;
extern const std::string neg_virial_axis_str;
// observable units
extern const std::string units_energy;
// container classs
class LJVerletTestContainer :
    public ObservableContainer {
private:
    // scalar observables
    ScalarObservable _kinetic_energy;
    ScalarObservable _potential_energy;
    ScalarObservable _energy;
    ScalarObservable _neg_virial;
public:
    // a constructor and a destructor
    LJVerletTestContainer();
    ~LJVerletTestContainer();
};
// state for the simulation
extern double atomic_mass;
// simulation class
class LJVerletTestSimulation :
    public Simulation {
public:
    // specifies initial distance between the two atoms in units of sigma
    void initialize_state(double initial_distance);
    // inherits EVOLVE from the base class
    virtual void calculate_remaining_observables();
    // initializes base class via
    // Simulation(name, integrator, observables)
    LJVerletTestSimulation(std::string name,
        std::string outdir,
        Integrator &integrator,
        ObservableContainer &observables);
    ~LJVerletTestSimulation();
};
#endif
