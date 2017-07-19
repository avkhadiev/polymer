// 2017 Artur Avkhadiev
/*! \file simulation.h
*/
#ifndef POLYMER_SIMULATION_H
#define POLYMER_SIMULATION_H
#include <vector>
#include <string>
#include <map>
#include "state.h"
#include "observable_container.h"
#include "ljpotential.h"
#include "force_updater.h"
#include "integrator.h"
class Simulation {
    protected:
        std::string _name;
        Integrator& _integrator;
        ObservableContainer& _observables;
        State _state;
        std::string _outdir;
        double _timestep;
        double _measurestep;
        int _cycle;
        int _observations_before_writeout;
    public:
        // getters
        std::string get_outdir();
        std::string get_name();
        const LJPotential& get_potential() const;
        const Integrator& get_integrator() const;
        ObservableContainer& get_observables() const;
        const State& get_state() const;
        int get_cycle();
        double get_time();
        double get_timestep();
        double get_measurestep();
        // setters
        void set_outdir(std::string outdir);
        void set_name(std::string name);
        void set_potential(LJPotential &potential);
        void set_integrator(Integrator &integrator);
        void set_observables(ObservableContainer &observables);
        void set_time(double time);
        void set_timestep(double timestep);
        void set_measurestep(double measurestep);
        // writeout
        // writes out observables to oudir/sim_name_<observable_name>.dat
        // if the vector of names is empty, outputs all observables.
        virtual void writeout_observables_to_file(std::vector<std::string> names,
            std::string outdir,
            bool overwrite = false);
        // main member functions
        virtual void evolve(int ncycles);
        // calculate observables after the integration step.
        // note that some observables are computed in the force loop or
        // during integration
        virtual void calculate_remaining_observables();
        // constructors and a destructor
        Simulation(std::string name,
            std::string outdir,
            Integrator& integrator,
            ObservableContainer &observables,
            int observations_before_writeout = 0);
        ~Simulation();
};
#endif
