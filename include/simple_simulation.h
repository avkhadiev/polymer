// 2017 Artur Avkhadiev
/*! \file simple_simulation.h
*/
#ifndef POLYMER_SIMPLE_SIMULATION_H
#define POLYMER_SIMPLE_SIMULATION_H
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include "parsing.h"
#include "simple_state.h"
#include "observable_container.h"
#include "ljpotential.h"
#include "force_updater.h"
#include "rattle_integrator.h"
#include "config_handler.h"
namespace simple{
    class Simulation {
    protected:
        std::string _name;      /**> run title, makes up file names */
        std::string _cndir;     /**> directory of config file       */
        std::string _tpdir;     /**> directory of `tape` file       */
        std::string _dtdir;     /**> directory for observables      */
        std::string _infile;    /**> input config file address      */
        std::string _cnfname;   /**> determined by _name            */
        std::string _tpfname;   /**> determined by _name            */
        ConfigHandler &_cfg;
        Integrator &_int;
        ObservableContainer &_obs;
        /**
        * contains
        * running averages, instantaneous values, and their evolution in time
        */
        double _dt;            /**> timestep in LJ units */
        size_t _step;          /**> number of steps made */
        size_t _calcstep;      /**> number of steps for calculating averages  */
        size_t _icalc;         /**> calculate observables every icalc steps   */
        size_t _iprint;        /**> writeout status stdout every iprint steps */
        size_t _isave;         /**> writeout to configfile every isave steps  */
        size_t _idata;         /**> writeout observables every idata steps    */
        size_t _itape;         /**> writeout to tapefile every itape steps    */
        /**
        * calculates any observables not already computed in the process of
        * numerical EOM integration
        */
        virtual void _calculate();
        /**
        * Readin utilities
        */
        bool _is_input_given;             /**> should input config be read? */
        void _read_config();              /**> read state + accumulators    */
        /**
        * Writeout utilities from within evolve()
        */
        void _write_data();             /**> writeout data                */
        void _write_status();           /**> report status to stdout      */
        void _write_config();           /**> overwrite state + accumulators */
        // append state to tape file
        void _prepare_tpstream(std::ofstream& tpstream);
        void _write_tape(std::ofstream& tpstream);
    public:
        bool is_input_given() {return _is_input_given;};
        void evolve(double runtime);
        // constructors and a destructor
        Simulation(std::string name,
            std::string cndir,
            std::string tpdir,
            std::string dtdir,
            std::string infile,
            ConfigHandler& config_handler,
            Integrator& integrator,
            ObservableContainer& container,
            double dt = 0.001,
            size_t icalc = 0,
            size_t iprint = 0,
            size_t isave = 0,
            size_t idata = 0,
            size_t itape = 0);
        ~Simulation();
    };
} // namespace simple
#endif
