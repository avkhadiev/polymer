// 2017 Artur Avkhadiev
/*! \file lj_verlet.cpp
* simple simulation test to test the Verlet integrator and the force loop with
* LJ potential.
* A diatomic molecule starts with its center of mass at the origin, aligned
* with the x-axis.
* The simulation evolves the system and tracks the change in potential energy, kinetic energy, and the virial. It outputs corresponding CSV files to the test/
 folder
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
#include "../include/rattle_integrator.h"
#include "../include/lj_rattle_test.h"
int ncycles = 100000;
// in LJ units
double timestep = pow(10, -3.0);
double measurestep = pow(10, -2.0);
// default output directory
const std::string default_outdir = "/Users/Arthur/stratt/polymer/test/";
// number of observations to store in memory before writeout
const int observations_before_writeout = 10000;
// integrator settings
const double tol = pow(10, -3.0);
const double rvtol = pow(10, -3.0);
const double tiny = pow(10, -7.0);
const int maxiter = pow(10, 5.0);
// observable names
const std::string kinetic_energy_str = "Kinetic Energy";
const std::string potential_energy_str = "Potential Energy";
const std::string energy_str = "Energy";
const std::string neg_virial_str = "Negative Virial";
const std::string neg_constraint_virial_str = "Negative Contraint Virial";
// observable axes labels
const std::string kinetic_energy_axis_str = "K";
const std::string potential_energy_axis_str = "V_{LJ}";
const std::string energy_axis_str = "E";
const std::string neg_virial_axis_str = "-W_{LJ}";
const std::string neg_constraint_virial_axis_str = "-W_{C}";
// observable units
const std::string units_energy = "\\epsilon";
// observable container class
LJRattleTestContainer::LJRattleTestContainer()
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
    _neg_constraint_virial
        = declare_scalar_observable(neg_constraint_virial_str,
            units_energy,
            neg_virial_axis_str);
    add_scalar_observable(&_kinetic_energy);
    add_scalar_observable(&_potential_energy);
    add_scalar_observable(&_energy);
    add_scalar_observable(&_neg_virial);
    add_scalar_observable(&_neg_constraint_virial);
};
LJRattleTestContainer::~LJRattleTestContainer(){
};
// mass unit is 1 atomic mass (all atoms are identical)
double atomic_mass = 1.0;
// simulation class
LJRattleTestSimulation::LJRattleTestSimulation(std::string name,
    std::string outdir,
    Integrator &integrator,
    ObservableContainer &observables) :
    Simulation(name, outdir, integrator, observables, observations_before_writeout)
{
    _timestep = timestep;
    _measurestep = measurestep;
}
LJRattleTestSimulation::~LJRattleTestSimulation(){
}
void LJRattleTestSimulation::initialize_state(double bond_length){
    // one atom is centered at (0, 0, 0), another one is a given distance away
    // along the x-axis
    double half_bond = bond_length / 2.0;
    Vector pos1 = vector(half_bond, 0.0, 0.0);
    Vector pos2 = vector(-half_bond, 0.0, 0.0);
    Vector vel1 = vector(0.0, 0.0, 0.0);
    Vector vel2 = vector(0.0, 0.0, 0.0);
    Atom atom1 = initialize_atom(atomic_mass, pos1, vel1);
    Atom atom2 = initialize_atom(atomic_mass, pos2, vel2);
    std::vector<Atom> atoms;
    atoms.push_back(atom1);
    atoms.push_back(atom2);
    Bond bond = initialize_bond(0, 1, pow(bond_length, 2.0));
    std::vector<Bond> bonds;
    bonds.push_back(bond);
    Molecule molecule = initialize_molecule(atoms, bonds);
    std::vector<Molecule> molecules;
    molecules.push_back(molecule);
    _state = ::initialize_state(molecules, 0.0);
    // evaluate forces on the two atoms
    bool calculate_observables = true;
    _integrator.get_force_updater().update_forces(_state,
        calculate_observables);
}
void LJRattleTestSimulation::calculate_remaining_observables(){
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
            "Usage: ./lj_rattle_test initial_separation optional <outdir>");
    }
    else {
        // declare the title of the study
        std::string sim_name = "lj_rattle";
        // save initial distance
        double bond_length = atof(argv[1]);
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
        LJRattleTestContainer container = LJRattleTestContainer();
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
        double *kinetic_energy_acc = NULL;
        double *neg_constraint_virial_acc = NULL;
        RattleIntegrator rattle = RattleIntegrator(force_updater,
            tol,
            rvtol,
            tiny,
            maxiter,
            kinetic_energy_acc,
            neg_constraint_virial_acc);
        Integrator& integrator = rattle;
        printf("%s\n", "force loop and integrator are set.");
        // set up the simulation
        printf("%s\n", "setting up the simulation...");
        LJRattleTestSimulation simulation
            = LJRattleTestSimulation(sim_name, outdir, integrator, observables);
        printf("%s\n", "simulation is set.");
        // initialize the state
        printf("%s\n", "initializing the state...");
        simulation.initialize_state(bond_length);
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
