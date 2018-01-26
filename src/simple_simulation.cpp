// 2017 Artur Avkhadiev
/*! \file simple_simulation.cpp
*/
#include <vector>
#include <sys/stat.h>                       /**> checks for file existence */
#include <math.h>
#include <string>
#include <cassert>
#include <stdexcept>
#include <map>
#include <sstream>
#include <fstream>
#include "../include/default_macros.h"
#include "../include/simple_simulation.h"
namespace simple{
    Simulation::Simulation(std::string name,
        std::string cndir,
        std::string tpdir,
        std::string dtdir,
        std::string infile,
        ConfigHandler& config_handler,
        Integrator& integrator,
        ObservableContainer& container,
        bool should_write_data,
        double dt,
        size_t icalc,
        size_t iblock,
        size_t iprint,
        size_t isave,
        size_t itape) :
        _name (parse_string(name)),
        _cndir(cndir),
        _tpdir(tpdir),
        _dtdir(dtdir),
        _infile(infile),
        _cnfname(_name + "_cn.txt"),
        _tpfname(_name + "_tp.txt"),
        _cfg(config_handler),
        _int(integrator),
        _obs(container),
        _should_write_data(should_write_data),
        _dt(dt),
        _step(0),
        _calcstep(0),
        _icalc(icalc),
        _iblock(iblock),
        _iprint(iprint),
        _isave(isave),
        _itape(itape)
    {
        if (_dt == 0.0){
            fprintf(stderr, "%s. %s: %f\n",
                "Cannot have zero timestep",
                "Setting default value", 0.001);
            _dt = 0.001;
        }
        //*********************************************************************
        if (_infile == ""){
            _is_input_given = false;
        }
        else {
            _is_input_given = true;
        }
        //*********************************************************************
        _prepare_observables();
        //*********************************************************************
        fprintf(stdout, "%s\n", "Simulation is set up:");
        fprintf(stdout, "%s: %s\n", "Name", _name.c_str());
        fprintf(stdout, "%s\n", _cfg.get_info_str().c_str());
        ForceUpdater& fupd = integrator.get_force_updater();
        if (_cfg.nmolecules() > 0){
            fprintf(stdout, "%s:\n%s", "Polymer-Polymer potential",
                fupd.get_polymer_potential()->get_str().c_str());
        }
        if (_cfg.nsolvents() > 0){
            fprintf(stdout, "%s:\n%s", "Solvent-Solvent potential",
                fupd.get_solvent_potential()->get_str().c_str());
        }
        if ((_cfg.nsolvents() > 0) and (_cfg.nmolecules() > 0)){
            fprintf(stdout, "%s:\n%s", "Polymer-Solvent potential",
                fupd.get_inter_potential()->get_str().c_str());
        }
        if (_is_input_given){
            fprintf(stdout, "%s: %s\n",
                "Input Config File",
                _infile.c_str());
            _read_config();                 // may change state
        }
        if (VERBOSE) {
            fprintf(stdout, "%s: %s%s\n",
                "Config File",
                _cndir.c_str(),
                _cnfname.c_str());
            fprintf(stdout, "%s: %s%s\n",
                "Tape File",
                _tpdir.c_str(),
                _tpfname.c_str());
            fprintf(stdout, "%s: %d\n", "Should write data (0/1)?", _should_write_data);
            fprintf(stdout, "%s: %f\n", "Timestep", _dt);
            fprintf(stdout, "%s: %lu time steps\n", "Calculation Interval", _icalc);
            fprintf(stdout, "%s: %lu calc steps\n", "Block and Data Write Interval", _iblock);
            fprintf(stdout, "%s: %lu time steps\n", "Status Interval", _iprint);
            fprintf(stdout, "%s: %lu time steps\n", "Config Save Interval", _isave);
            fprintf(stdout, "%s: %lu time steps\n", "Tape Interval", _itape);
        }
    }
    Simulation::~Simulation(){}
    void Simulation::_prepare_observables(){
        if (_icalc == 0) {
            // if no instantaneous values are to be calculated...
            // there will be no block averaging either
            _iblock = 0;
        }
        if (_iblock == 0) {
            // if there is no block averaging...
            // iblock is still used to regulate output
            // turn off any averaging
            _obs.set_average_data(false);
        }
        if (_icalc != 0) _obs.run_begin(_dtdir, _name, _should_write_data);
    }
    void Simulation::_write_status(){
        fprintf(stdout, "%s\n", "STATUS:");
        fprintf(stdout, "%s: %f\n", "Time", _cfg.atom_state().time());
        fprintf(stdout, "%s", _obs.status_string(VERBOSE).c_str());
    }
    void Simulation::_read_config(){
        std::string fin = _infile;
        std::ifstream readout;
        readout.open(fin, std::ifstream::in);
        if (!readout.is_open()) {
            std::string err_msg = "read_config: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
            perror("open");
        }
        else{
            // get state
            _cfg.read_bond_state(readout);
        }
        readout.close();
    }
    void Simulation::_write_config(std::ofstream& writeout){
        // open file in truncate mode,
        // write out state,
        // write out accumulators
        bool verbose = false;
        bool output_header = false;
        if (writeout.is_open()) {
            writeout << _cfg.bond_state().to_string(verbose, output_header);
        }
        else {
            perror("open");
        }
    }
    void Simulation::_prepare_outstream(std::ofstream &stream,
        std::string fout, bool truncate){
        bool verbose = false;
        // check for file existence
        struct stat buffer;
        if (stat(fout.c_str(), &buffer) == 0) {
            if (truncate){
                // if requested, overwrite the information
                stream.open(fout, std::ofstream::out | std::ofstream::trunc);
            }
            else {
                // otherwise append the information
                stream.open(fout, std::ofstream::out | std::ofstream::app);
            }
        }
        // if doesn't exist, create new one
        else{
            truncate = true;                    //**> now will truncate */
            fprintf(stderr, "%s\n",
                "stream file doesn't exist, creating new one.");
            stream.open(fout, std::ofstream::out | std::ofstream::trunc);
        }
        // now file has to be opened
        if (!stream.is_open()) {
            // if file still could not be opened
            std::string err_msg = "write_state_to_file: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
            perror("open");
        }
        else{
            // output header information if truncation was done
            if (truncate){
                stream << _cfg.atom_state().header_str(verbose) << std::endl;
            }
        }
    }
    void Simulation::_write_data(){
        _obs.write_data(_dtdir, _name);
    }
    void Simulation::write_run_summary(){
        _obs.write_run_summary(_dtdir, _name);
    }
    void Simulation::_calculate_observables(size_t iblock){
        // if this is the first step after
        // averaging the previous block, start a new block.
        if (_calcstep % iblock == 1){
            ++_blockstep;
            _obs.block_begin();
        }
        // calculate whatever instantaneous values are to be computed
        // in the main loop
        _obs.calculate_observables(_cfg.atom_state());
        // update block average accumulators with
        // the new instantaneous values
        _obs.update_block();
        if (_calcstep % iblock == 0){
            // if this is the end of the block,
            // perform all the necessary averaging and record
            _obs.block_end();
        }
    }
    void Simulation::_record_observables(){
        _obs.record_observables();
    }
    void Simulation::relax(double relax_time){
        // do not write to trajectory or configuration files during relaxation
        size_t isave_og = _isave;
        size_t itape_og = _itape;
        _isave = 0;
        _itape = 0;
        bool should_avg_og = _obs.should_average();
        _obs.set_average_data(false);
        fprintf(stdout, "%s %f\n", "RELAXATION TIME: ", relax_time);
        Simulation::evolve(relax_time);
        _obs.zero_accumulators();
        _obs.set_average_data(should_avg_og);
        _isave = isave_og;
        _itape = itape_og;
    }
    void Simulation::evolve(double runtime){
        if (runtime < 0){
            fprintf(stderr, "%s (%f), %s\n",
                "runtime is < 0", runtime, "will use absolute value.");
            runtime = abs(runtime);
        }
        size_t icalc = _icalc;
        size_t iblock = _iblock;
        size_t iprint = _iprint;
        size_t isave = _isave;
        size_t itape = _itape;
        // _i<interval> == 0 => don't perform the action, but evolve() will
        // use modular arithmetic: if _step % i<interval> == 0
        size_t nsteps = _step + (size_t)(trunc(runtime / _dt));
        if (icalc == 0) icalc = nsteps + 1;
        if (iblock == 0) {
            // if there is no block averaging...
            // iblock is still used to regulate output
            // iblock = 100 means that data will be grouped in blocks of 100
            // calc steps and will be we output in chunks of this size, or at
            // the end of the run, whatever happens first
            // make iblock = number of calcsteps per run + 1
            // turn off any averaging
            iblock = 100;
        }
        if (iprint == 0) iprint = nsteps + 1;
        if (isave == 0) isave = nsteps + 1;
        if (itape == 0) itape = nsteps + 1;
        // write out initial state data as necessary
        std::ofstream tpstream;
        std::ofstream cfstream;
        if (_itape != 0){
            _prepare_tpstream(tpstream);
            _write_config(tpstream);
        }
        if (_isave != 0) _prepare_cfstream(cfstream);
        bool calc;
        _calcstep = 0;
        _blockstep = 0;
        while(_step < nsteps){
            _step = _step + 1;
            _step % icalc == 0? calc = true : calc = false;
            _int.move(_dt, _cfg.atom_state(), calc);
            if (calc) {
                ++_calcstep;
                // manage observable calculation and averaging
                _calculate_observables(iblock);
                // printing updates and saving a configuration will only happen
                // if observables are being calculated
                if (_step % iprint == 0) _write_status();
                if (_step % isave == 0) _write_config(cfstream);
            }
            if (_step % itape == 0) _write_config(tpstream);
            if (calc){
                _record_observables();
                if ((_should_write_data) and (_calcstep % iblock == 0)) {
                    _write_data();
                }
            }
        }
        if (_itape != 0){
            // write to files if this state has not been already written,
            // to save the ultimate configuration after the run
            if (_step % itape != 0) _write_config(tpstream);
            tpstream.close();
            fprintf(stdout, "%s\n", "Tape file closed.");
        }
        if (_isave != 0){
            cfstream.close();
            fprintf(stdout, "%s\n", "Config file closed.");
        }
        _obs.run_end();
    }
} // namespace simple
