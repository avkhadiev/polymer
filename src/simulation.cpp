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
    std::string outdir,
    Integrator& integrator,
    ObservableContainer &observables,
    int observations_before_writeout) :
    _name (name),
    _integrator (integrator),
    _observables (observables),
    _outdir (outdir),
    _observations_before_writeout (observations_before_writeout)
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
std::string Simulation::get_outdir(){
    return _outdir;
}
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
void Simulation::set_outdir(std::string outdir){
    _outdir = outdir;
}
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
    // make sure steps are not 0
    if (_timestep == 0){
        std::string err_msg = "evolve: timestep is equal to 0";
        throw std::invalid_argument(err_msg);
    }
    if (_measurestep == 0){
        std::string err_msg = "evolve: measurestep is equal to 0";
        throw std::invalid_argument(err_msg);
    }
    // find how many timesteps are in the measure step (when to write out)
    int writeoutcycle = int((_measurestep / _timestep) + 0.5);
    printf("%s %d\n", "writeout cycle is", writeoutcycle);
    int nobservations = 0;
    bool overwrite_observables = true;
    bool overwrite_state = true;
    bool verbose_state = true;
    // iterate through the cycles, evolving the system and writing out
    // the observables
    bool calculate_observables;
    for(int icycle = 0; icycle < ncycles; ++icycle){
        if (icycle % writeoutcycle != 0){
            // integrate without computing observables
            calculate_observables = false;
            try
            {
                _integrator.move(_timestep, _state, calculate_observables);
            }
            catch (std::invalid_argument &e)
            {
                fprintf(stderr, "%s\n", e.what());
                fprintf(stderr, "%s\n",
                    state_to_string(_state, verbose_state).c_str());
                throw;
            }
        }
        else {
            // integrate with computing observables
            calculate_observables = true;
            write_state_to_file(_state, _outdir,
                _name,
                verbose_state,
                overwrite_state);
            // only overwrite state output file for the first time
            overwrite_state = false;
            try
            {
                _integrator.move(_timestep, _state, calculate_observables);
            }
            catch (std::invalid_argument &e)
            {
                fprintf(stderr, "%s\n", e.what());
                fprintf(stderr, "%s\n",
                    state_to_string(_state, verbose_state).c_str());
                throw;
            }
            // update accumulators of other observables that were not computed
            // in the force loop / integration step
            calculate_remaining_observables();
            // record all accumulators
            _observables.update_observables_through_accumulators({}, _state.time);
            // writeout if necessary
            nobservations += 1;
            if (nobservations > _observations_before_writeout) {
                writeout_observables_to_file({},
                    _outdir,
                    overwrite_observables);
                // do not overwrite after the first time
                overwrite_observables = false;
                // reset the counter
                nobservations = 0;
                // empty observable_container
                _observables.clear_observables_records();
            }
        }
        _cycle += 1;
    }
    // writeout the remaining observations
    writeout_observables_to_file({},
        _outdir,
        overwrite_observables);
    fprintf(stdout, "%s %s\n",
        "wrote out csv data to",
        _outdir.c_str());
}
void Simulation::calculate_remaining_observables(){
}
