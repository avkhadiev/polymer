// 2017 Artur Avkhadiev
/*! \file simulation.h
*/
#ifndef POLYMER_SIMULATION_H
#define POLYMER_SIMULATION_H
#include <vector>
#include <map>
#include "state.h"
#include "observable_container.h"
#include "ljpotential.h"
#include "integrator.h"
class Simulation {
    private:
        LJPotential& _potential;
        Integrator& _integrator;
        ObservableContainer& _observables;
        State &_state;
        typedef struct indicators_t {
            int step;
        } _Indicators;
        typedef struct parameters_t {
            double timestep;
            double measurestep;
        } _Parameters;
    public:
        // Getters
        LJPotential& get_potential();
        Integrator& get_integrator();
        ObservableContainer& get_observables();
        State& get_state();
        double get_time();
        double get_step();
        double get_timestep();
        double get_measurestep();
        // Setters
        void set_potential(LJPotential &potential);
        void set_integrator(Integrator &integrator);
        void set_observables(ObservableContainer &observables);
        void set_time(double time);
        void set_step(int step);
        void set_timestep(double timestep);
        void set_measurestep(double measurestep);
        // constructors and a destructor
        Simulation();
        Simulation(LJPotential& potential,
            ObservableContainer &observables,
            double timestep,
            double measurestep);
        ~Simulation();
};
#endif
