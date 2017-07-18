// 2017 Artur Avkhadiev
/*! \file simulation.cpp
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
Simulation::Simulation(std::string name,
    Integrator& integrator,
    ObservableContainer &observables) :
    _name (name),
    _integrator (integrator),
    _observables (observables)
{
    // default state contains no atoms or molecules
    std::vector<Molecule> default_molecules = {};
    _state = {.time = 0.0, .nm = 0, .molecules = default_molecules};
    // all counters are zero
    _timestep = 0.0;
    _measurestep  = 0.0;
    _cycle = 0;
};
Simulation::~Simulation(){};
// Getters
std::string Simulation::get_name(){
    return _name;
}
const LJPotential& Simulation::get_potential() const {
    return _integrator.get_force_updater().get_potential();
}
const Integrator& Simulation::get_integrator() const {
    return _integrator;
}
ObservableContainer& Simulation::get_observables() const {
    return _observables;
}
const State& Simulation::get_state() const {
    return _state;
}
int Simulation::get_cycle() {
    return _cycle;
};
double Simulation::get_time() {
    return _state.time;
}
double Simulation::get_timestep() {
    return _timestep;
}
double Simulation::get_measurestep() {
    return _measurestep;
}
// Setters
void Simulation::set_name(std::string name){
    _name = name;
};
void Simulation::set_potential(LJPotential &potential){
    _integrator.get_force_updater().set_potential(potential);
}
void Simulation::set_integrator(Integrator &integrator){
    _integrator = integrator;
}
void Simulation::set_observables(ObservableContainer &observables){
    _observables = observables;
}
void Simulation::set_time(double time){
    // boy should I have used namespaces! :(
    ::set_time(&_state, time);
}
void Simulation::set_timestep(double timestep){
    _timestep = timestep;
}
void Simulation::set_measurestep(double measurestep){
    _measurestep = measurestep;
};
void Simulation::writeout_observables_to_file(std::vector<std::string> names, std::string outdir, bool overwrite)
{
    _observables.writeout_observables_to_file(names,
        outdir, _name, overwrite);
}
void Simulation::evolve(int ncycles){
}
