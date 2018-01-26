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
#include "force_updater.h"
#include "rattle_integrator.h"
#include "config_handler.h"
#include "settings_parser.h"
namespace simple{
    class Simulation {
    protected:
        std::string _name;      /**> run title, makes up file names */
        std::string _cndir;     /**> directory of config file       */
        std::string _tpdir;     /**> directiry for ``tape'' file    */
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
        bool _should_write_data;
        double _dt;            /**> timestep in LJ units */
        size_t _step;          /**> number of steps made */
        size_t _calcstep;      /**> number of steps for calculating averages  */
        size_t _blockstep;     /**> number of blocks calculated */
        size_t _icalc;         /**> calculate observables every icalc steps   */
        size_t _iblock;        /**> average block every iblock icalc steps    */
        size_t _iprint;        /**> writeout status every iprint icalc steps  */
        size_t _isave;         /**> write config every isave icalc steps      */
        size_t _itape;         /**> writeout to tape file every itape steps   */
        /**
        * if should_write_data: setup datafile streams --- instantaneous (& average, if neccessary)
        * calls start_run on each observable
        */
        void _prepare_observables();
        /**
        * calculates any instantaneous observables
        * not already computed in the process of
        * numerical EOM integration
        */
        void _calculate_observables(size_t iblock);
        /**
        * records all instantanoues observables in the observable container
        */
        void _record_observables();
        /**
        * Readin utilities
        */
        bool _is_input_given;             /**> should input config be read? */
        void _read_config();              /**> read state configuration */
        /**
        * Writeout utilities from within evolve()
        */
        void _prepare_outstream(std::ofstream& stream, std::string fout, bool truncate);
        void _prepare_tpstream(std::ofstream &stream){
            std::string fout = _tpdir + _tpfname;
            bool truncate = false;
            _prepare_outstream(stream, fout, truncate);
            fprintf(stdout, "%s\n", "Tape file opened for writing.");
        };
        void _prepare_cfstream(std::ofstream &stream){
            std::string fout = _cndir + _cnfname;
            bool truncate = true;
            _prepare_outstream(stream, fout, truncate);
            fprintf(stdout, "%s\n", "Config file opened for writing.");
        };
        void _write_data();             /**> writeout data                */
        void _write_run_summary();      /**> writeut run averages         */
        void _write_status();           /**> report status to stdout      */
        void _write_config(std::ofstream& stream);
        void _write_trajectory(std::ofstream& stream);
    public:
        void set_should_write_data(bool should) {_should_write_data = should;};
        bool should_write_data() const {return _should_write_data;};
        bool is_input_given() {return _is_input_given;};
        //virtual void initialize();
        //void cool(double temp_change, double cool_time);
        void relax(double relax_time);
        void evolve(double runtime);
        void write_run_summary();           /**> writeut run averages      */
        // constructors and a destructor
        // TODO just pass the settings parser instead
        Simulation(std::string name,
            std::string cndir,
            std::string tpdir,
            std::string dtdir,
            std::string infile,
            ConfigHandler& config_handler,
            Integrator& integrator,
            ObservableContainer& container,
            bool should_write_data = true,
            double dt = 0.001,
            size_t icalc = 10,          // calculate every 10 steps
            size_t iblock = 100,        // average every 100 calcsteps
            size_t iprint = 1000,       // print status every 1000 steps
            size_t isave = 1000,        // save config + obs every 1000 steps
            size_t itape = 100);        // save config every 100 steps
        ~Simulation();
    };
} // namespace simple
#endif
