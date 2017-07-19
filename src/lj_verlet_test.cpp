// 2017 Artur Avkhadiev
/*! \file lj_verlet.cpp
* simple simulation test to test the Verlet integrator and the force loop with
* LJ potential.
* Two atoms start on one line with the given initial separation, with zero velocity => total energy is known.
* The simulation evolves the system and tracks the change in potential energy, kinetic energy, and the virial. It outputs corresponding CSV files to the test/
 foler
*/
#include <iostream>
#include <cmath>                            /* pow */
#include <string>
#include <stdexcept>
#include "../include/vector.h"
#include "../include/ljpotential.h"
#include "../include/force_updater.h"
#include "../include/observable_container.h"
#include "../include/integrator.h"
#include "../include/verlet_integrator.h"
#include "../include/lj_verlet_test.h"
int ncycles = 100000;
// in LJ units
double timestep = pow(10, -5.0);
double measurestep = pow(10, -4.0);
// default output directory
const std::string default_outdir = "/Users/Arthur/stratt/polymer/test/";
// number of observations to store in memory before writeout
const int observations_before_writeout = 100;
// observable names
const std::string kinetic_energy_str = "Kinetic Energy";
const std::string potential_energy_str = "Potential Energy";
const std::string energy_str = "Energy";
const std::string neg_virial_str = "Virial";
// observable axes labels
const std::string kinetic_energy_axis_str = "K";
const std::string potential_energy_axis_str = "V_{LJ}";
const std::string energy_axis_str = "E";
const std::string neg_virial_axis_str = "-W_{LJ}";
// observable units
const std::string units_energy = "\\epsilon";
// observable container class
LJVerletTestContainer::LJVerletTestContainer()
{
    _kinetic_energy = declare_scalar_observable(kinetic_energy_str,
        units_energy,
        kinetic_energy_axis_str);
    _potential_energy = declare_scalar_observable(potential_energy_str,
        units_energy,
        potential_energy_axis_str);
    _energy = declare_scalar_observable(energy_str,
        units_energy,
        energy_axis_str);
    _neg_virial = declare_scalar_observable(neg_virial_str,
        units_energy,
        neg_virial_axis_str);
    add_scalar_observable(&_kinetic_energy);
    add_scalar_observable(&_potential_energy);
    add_scalar_observable(&_energy);
    add_scalar_observable(&_neg_virial);
};
LJVerletTestContainer::~LJVerletTestContainer(){
};
// mass unit is 1 atomic mass (all atoms are identical)
double atomic_mass = 1.0;
// simulation class
LJVerletTestSimulation::LJVerletTestSimulation(std::string name,
    std::string outdir,
    Integrator &integrator,
    ObservableContainer &observables) :
    Simulation(name, outdir, integrator, observables, observations_before_writeout)
{
    _timestep = timestep;
    _measurestep = measurestep;
}
LJVerletTestSimulation::~LJVerletTestSimulation(){
}
void LJVerletTestSimulation::initialize_state(double initial_distance){
    // one atom is centered at (0, 0, 0), another one is a given distance away
    // along the x-axis
    Vector pos1 = vector(0.0, 0.0, 0.0);
    Vector pos2 = add(pos1, vector(initial_distance, 0.0, 0.0));
    Vector vel = vector(0.0, 0.0, 0.0);
    Atom atom1 = initialize_atom(atomic_mass, pos1, vel);
    Atom atom2 = initialize_atom(atomic_mass, pos2, vel);
    std::vector<Bond> no_bonds = {};
    std::vector<Atom> atoms1;
    std::vector<Atom> atoms2;
    atoms1.push_back(atom1);
    atoms2.push_back(atom2);
    Molecule molecule1 = initialize_molecule(atoms1, no_bonds);
    Molecule molecule2 = initialize_molecule(atoms2, no_bonds);
    std::vector<Molecule> molecules;
    molecules.push_back(molecule1);
    molecules.push_back(molecule2);
    _state = ::initialize_state(molecules, 0.0);
    // evaluate forces on the two atoms
    bool calculate_observables = true;
    _integrator.get_force_updater().update_forces(_state,
        calculate_observables);
}
void LJVerletTestSimulation::calculate_remaining_observables(){
    // update accumulator as a sum of potential and kinetic energies
    double k = *(_observables.get_scalar_observable_accumulator(kinetic_energy_str));
    double v = *(_observables.get_scalar_observable_accumulator(potential_energy_str));
    double *e = _observables.get_scalar_observable_accumulator(energy_str);
    *e = k + v;
}
// takes in arguments: initial separation of the atom
int main(int argc, char **argv){
    std::string fname;
    if (argc < 2){
        fprintf(stderr,
            "%s\n",
            "Usage: ./lj_verlet_test initial_separation optional <outdir>");
    }
    else {
        // declare the title of the study
        std::string sim_name = "lj_verlet";
        // save initial distance
        double initial_distance = atof(argv[1]);
        std::string outdir;
        // second (optional) argument is the output directory
        if (argc == 3){
            outdir = std::string(argv[2]);
        }
        else {
            outdir = default_outdir;
        }
        // set up the observable container
        printf("%s\n", "setting up the container...");
        LJVerletTestContainer container = LJVerletTestContainer();
        ObservableContainer& observables = container;
        printf("%s\n", "container set.");
        // set up the force loop, tell it where are the accumulators for
        // potential energy and virial (they will be updated if specified)
        printf("%s\n", "setting up the force loop and integrator...");
        ForceUpdater force_updater = ForceUpdater(LJPotential(),
            observables.get_scalar_observable_accumulator(potential_energy_str),
            observables.get_scalar_observable_accumulator(neg_virial_str));
        // set up the integrator, tell it where is the accumulator for the
        // kinetic energy
        VerletIntegrator verlet = VerletIntegrator(force_updater,
            observables.get_scalar_observable_accumulator(kinetic_energy_str));
        Integrator& integrator = verlet;
        printf("%s\n", "force loop and integrator are set.");
        // set up the simulation
        printf("%s\n", "setting up the simulation...");
        LJVerletTestSimulation simulation
            = LJVerletTestSimulation(sim_name, outdir, integrator, observables);
        printf("%s\n", "simulation is set.");
        // initialize the state
        printf("%s\n", "initializing the state...");
        simulation.initialize_state(initial_distance);
        // evolve
        try
        {
            simulation.evolve(ncycles);
        }
        catch (std::invalid_argument &e)
        {
            fprintf(stdout, "ERROR occured. See error log.\n");
        }
    }
}
