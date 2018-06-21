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
#include "geodesic_path_computer.h"
namespace simple{
    class Simulation {
    protected:
        std::string _name;      /**> run title, makes up file names */
        std::string _cndir;     /**> directory of config file       */
        std::string _tpdir;     /**> directiry for ``tape'' file    */
        std::string _dtdir;     /**> directory for observables      */
        std::string _cnfname;   /**> determined by _name            */
        std::string _tpfname;   /**> determined by _name            */
        ObservableContainer &_obs;
        /**
        * contains
        * running averages, instantaneous values, and their evolution in time
        */
        bool _should_write_data;
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
        void _calculate_observables(size_t iblock, simple::AtomState& state);
        /**
        * records all instantanoues observables in the observable container
        */
        void _record_observables();
        /**
        * Readin utilities
        */
        void _read_config();      /**> read state configuration */
        /**
        * Writeout utilities
        */
        virtual void _prepare_outstream(std::ofstream& stream, std::string fout, bool truncate);
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
        virtual void _write_status() = 0;  /**> report status to stdout      */
        virtual void _write_config(std::ofstream& stream) = 0;
    public:
        void set_should_write_data(bool should) {_should_write_data = should;};
        bool should_write_data() const {return _should_write_data;};
        void write_run_summary();           /**> writeut run averages      */
        // constructors and a destructor
        // TODO just pass the settings parser instead
        Simulation(std::string name,
            std::string cndir,
            std::string tpdir,
            std::string dtdir,
            ObservableContainer& container,
            bool should_write_data = true,
            size_t icalc = 10,          // calculate every 10 steps
            size_t iblock = 100,        // average every 100 calcsteps
            size_t iprint = 1000,       // print status every 1000 steps
            size_t isave = 1000,        // save config + obs every 1000 steps
            size_t itape = 100);        // save config every 100 steps
        ~Simulation();
    };
    class MDSimulation : public Simulation {
    protected:
        std::string _infile;              /**> input config file address    */
        ConfigHandler &_cfg;
        Integrator &_int;
        double _dt;                       /**> timestep in LJ units         */
        /**
        * Readin utilities
        */
        bool _is_input_given;             /**> should input config be read? */
        void _read_config();              /**> read state configuration     */
        virtual void _prepare_outstream(std::ofstream& stream, std::string fout, bool truncate);
        virtual void _write_status();     /**> report status to stdout      */
        virtual void _write_config(std::ofstream& stream);
        // TODO
        void _write_trajectory(std::ofstream& stream);
    public:
        bool is_input_given() {return _is_input_given;};
        //virtual void initialize();
        //void cool(double temp_change, double cool_time);
        void relax(double relax_time);
        void evolve(double runtime);
        MDSimulation(std::string name,
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
        ~MDSimulation();
    };
    class GeodesicSimulation : public Simulation {
    protected:
        geodesic::Record _initial;        /**> initial record            */
        geodesic::Record _final;          /**> final record              */
        geodesic::Path _path;             /**> geodesic path             */
        std::string _pfname;
        geodesic::PathComputer* _comp;    /**> path computer             */

        double _el;                       /**> Landscape energy          */
        double _dtau;                     /**> step in progress variable */
        size_t _maxiter;
        size_t _max_escape_iter;
        /* create path file and trajectory file, output header line */
        virtual void _prepare_outstream(std::ofstream& stream,
            std::string fout, bool truncate);
        virtual void _write_status();     /**> report status to stdout      */
        /* append information to files already created */
        virtual void _write_config(std::ofstream& stream);
        // TODO
        void _write_trajectory(std::ofstream& stream);
        void _prepare_pfstream(std::ofstream &stream){
            std::string fout = _cndir + _pfname;
            bool truncate = true;
            _prepare_outstream(stream, fout, truncate);
            fprintf(stdout, "%s\n", "Path file opened for writing.");
        };
    public:
        observable::PE _pe;
        // temporary observables to check algorithm implementation
        geodesic::OmegaProj omega_proj1;
        geodesic::OmegaProj omega_proj2;
        geodesic::Path path() const {return _path;};
        double landscape_energy() const {return _el;};
        void set_landscape_energy(double el) {_el = el;};
        void compute_path();
        GeodesicSimulation(std::string name,
            std::string initial,
            std::string final,
            std::string cndir,
            std::string tpdir,
            std::string dtdir,
            ObservableContainer& container,
            geodesic::PathComputer* comp,
            double landscape_energy,
            bool should_write_data = true,
            double dtau = 0.001,
            size_t icalc = 10,          // calculate every 10 steps
            size_t iblock = 100,        // average every 100 calcsteps
            size_t iprint = 1000,       // print status every 1000 steps
            size_t isave = 1000,        // save config + obs every 1000 steps
            size_t itape = 100,         // save config every 100 steps
            size_t maxiter = 100000000,        // max propagation steps
            size_t max_escape_iter = 10000);   // max escape steps
        ~GeodesicSimulation();
    };
    class SLERPGeodesicSimulation : public GeodesicSimulation {
    public:
        SLERPGeodesicSimulation(std::string name,
            std::string initial,
            std::string final,
            std::string cndir,
            std::string tpdir,
            std::string dtdir,
            ObservableContainer& container,
            geodesic::SLERP* comp,
            double landscape_energy,
            bool should_write_data = true,
            double dtau = 0.001,
            size_t icalc = 10,          // calculate every 10 steps
            size_t iblock = 100,        // average every 100 calcsteps
            size_t iprint = 1000,       // print status every 1000 steps
            size_t isave = 1000,        // save config + obs every 1000 steps
            size_t itape = 100,         // save config every 100 steps
            size_t maxiter = 100000000,        // max propagation steps
            size_t max_escape_iter = 10000);   // max escape steps
        ~SLERPGeodesicSimulation();
    };
    class ShortStepGeodesicSimulation : public GeodesicSimulation {
    public:
        ShortStepGeodesicSimulation(std::string name,
            std::string initial,
            std::string final,
            std::string cndir,
            std::string tpdir,
            std::string dtdir,
            ObservableContainer& container,
            geodesic::ShortStep* comp,
            double landscape_energy,
            bool should_write_data = true,
            double dtau = 0.001,
            size_t icalc = 10,          // calculate every 10 steps
            size_t iblock = 100,        // average every 100 calcsteps
            size_t iprint = 1000,       // print status every 1000 steps
            size_t isave = 1000,        // save config + obs every 1000 steps
            size_t itape = 100,         // save config every 100 steps
            size_t _maxiter = 500000,   // max propagation steps
            size_t _max_escape_iter = 10000);   // max escape steps
        ~ShortStepGeodesicSimulation();
    };
    class ShoveGeodesicSimulation : public GeodesicSimulation {
    public:
        ShoveGeodesicSimulation(std::string name,
            std::string initial,
            std::string final,
            std::string cndir,
            std::string tpdir,
            std::string dtdir,
            ObservableContainer& container,
            geodesic::SHOVE* comp,
            double landscape_energy,
            bool should_write_data = true,
            double sigma = 0.005,       // max step in configuration space
            size_t icalc = 10,          // calculate every 10 steps
            size_t iblock = 100,        // average every 100 calcsteps
            size_t iprint = 1000,       // print status every 1000 steps
            size_t isave = 1000,        // save config + obs every 1000 steps
            size_t itape = 100,         // save config every 100 steps
            size_t _maxiter = 500000,   // max propagation steps
            size_t _max_escape_iter = 10000);   // max escape steps
        ~ShoveGeodesicSimulation();
    };
    class PlerpGeodesicSimulation : public GeodesicSimulation {
    public:
        PlerpGeodesicSimulation(std::string name,
            std::string initial,
            std::string final,
            std::string cndir,
            std::string tpdir,
            std::string dtdir,
            ObservableContainer& container,
            geodesic::PLERP* comp,
            double landscape_energy,
            bool should_write_data = true,
            double sigma = 0.005,       // max step in configuration space
            size_t icalc = 10,          // calculate every 10 steps
            size_t iblock = 100,        // average every 100 calcsteps
            size_t iprint = 1000,       // print status every 1000 steps
            size_t isave = 1000,        // save config + obs every 1000 steps
            size_t itape = 100,         // save config every 100 steps
            size_t _maxiter = 500000,   // max propagation steps
            size_t _max_escape_iter = 10000);   // max escape steps
        ~PlerpGeodesicSimulation();
    };
} // namespace simple
#endif
