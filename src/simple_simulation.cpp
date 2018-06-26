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
        ObservableContainer& container,
        bool should_write_data,
        size_t icalc,
        size_t iblock,
        size_t iprint,
        size_t isave,
        size_t itape) :
        _name (parse_string(name)),
        _cndir(cndir),
        _tpdir(tpdir),
        _dtdir(dtdir),
        _cnfname(_name + "_cn.txt"),
        _tpfname(_name + "_tp.txt"),
        _obs(container),
        _should_write_data(should_write_data),
        _step(0),
        _calcstep(0),
        _icalc(icalc),
        _iblock(iblock),
        _iprint(iprint),
        _isave(isave),
        _itape(itape)
    {
        _prepare_observables();
    }
    Simulation::~Simulation(){}
    MDSimulation::MDSimulation(std::string name,
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
        Simulation(name, cndir, tpdir, dtdir, container,
                    should_write_data, icalc, iblock, iprint, isave, itape),
        _infile(infile),
        _cfg(config_handler),
        _int(integrator),
        _dt(dt)
    {
        _cnfname = _name + "_cn.txt";
        _tpfname = _name + "_tp.txt";
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
        if (VERBOSE) {
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
    MDSimulation::~MDSimulation(){}
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
    void MDSimulation::_write_status(){
        fprintf(stdout, "%s\n", "STATUS:");
        fprintf(stdout, "%s: %f\n", "Time", _cfg.atom_state().time());
        fprintf(stdout, "%s", _obs.status_string(VERBOSE).c_str());
    }
    void MDSimulation::_read_config(){
        std::string fin = _infile;
        std::ifstream readout;
        readout.open(fin.c_str(), std::ifstream::in);
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
    void MDSimulation::_write_config(std::ofstream& writeout){
        // open file in truncate mode,
        // write out state,
        // write out accumulators
        bool verbose = false;
        bool output_header = false;
        if (writeout.is_open()) {
            writeout << _cfg.atom_state().to_string(verbose, output_header);
        }
        else {
            perror("open");
        }
    }
    void Simulation::_prepare_outstream(std::ofstream &stream,
        std::string fout, bool truncate){
        // check for file existence
        struct stat buffer;
        if (stat(fout.c_str(), &buffer) == 0) {
            if (truncate){
                // if requested, overwrite the information
                stream.open(fout.c_str(), std::ofstream::out | std::ofstream::trunc);
            }
            else {
                // otherwise append the information
                stream.open(fout.c_str(), std::ofstream::out | std::ofstream::app);
            }
        }
        // if doesn't exist, create new one
        else{
            truncate = true;                    //**> now will truncate */
            // fprintf(stderr, "%s\n",
            //    "stream file doesn't exist, creating new one.");
            stream.open(fout.c_str(), std::ofstream::out | std::ofstream::trunc);
        }
    }
    void MDSimulation::_prepare_outstream(std::ofstream &stream,
        std::string fout, bool truncate){
        bool verbose = false;
        Simulation::_prepare_outstream(stream, fout, truncate);
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
    void Simulation::_calculate_observables(size_t iblock,
                                            simple::AtomState& state){
        // if this is the first step after
        // averaging the previous block, start a new block.
        if (_calcstep % iblock == 1){
            ++_blockstep;
            _obs.block_begin();
        }
        // calculate whatever instantaneous values are to be computed
        // in the main loop
        _obs.calculate_observables(state);
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
    void MDSimulation::relax(double relax_time){
        // do not write to trajectory or configuration files during relaxation
        size_t isave_og = _isave;
        size_t itape_og = _itape;
        _isave = 0;
        _itape = 0;
        bool should_avg_og = _obs.should_average();
        _obs.set_average_data(false);
        if (VERBOSE) {
            fprintf(stdout, "%s %f\n", "RELAXATION TIME: ", relax_time);
        }
        evolve(relax_time);
        _obs.zero_accumulators();
        _obs.set_average_data(should_avg_og);
        _isave = isave_og;
        _itape = itape_og;
    }
    void MDSimulation::evolve(double runtime){
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
                _calculate_observables(iblock, _cfg.atom_state());
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
            //fprintf(stdout, "%s\n", "Tape file closed.");
        }
        if (_isave != 0){
            cfstream.close();
            //fprintf(stdout, "%s\n", "Config file closed.");
        }
        _obs.run_end();
    }
    GeodesicSimulation::GeodesicSimulation(std::string name,
        std::string ini_path,
        std::string fin_path,
        std::string cndir,
        std::string tpdir,
        std::string dtdir,
        ObservableContainer& container,
        geodesic::PathComputer* comp,
        double landscape_energy,
        bool should_write_data,
        double dtau,
        size_t icalc,               // calculate every 10 steps
        size_t iblock,              // average every 100 calcsteps
        size_t iprint,              // print status every 1000 steps
        size_t isave,               // save config + obs every 1000 steps
        size_t itape,               // save config every 100 steps
        size_t maxiter,             // max propagation steps
        size_t max_escape_iter) :   // max escape steps
        Simulation(name, cndir, tpdir, dtdir, container,
                    should_write_data, icalc, iblock, iprint, isave, itape),
        _initial(ini_path),
        _final(fin_path),
        _path(ini_path, fin_path),
        _pfname(_name + "_path.txt"),
        _comp(comp),
        _el(landscape_energy),
        _dtau(dtau),
        _maxiter(maxiter),
        _max_escape_iter(max_escape_iter),
        //  mean  err  print e_format
        _pe(false, false, false, false)
    {
        // add potential energy and omega projections as observables
        _obs.add_observable(&_pe);
        _obs.add_observable(_path.get_length());
        _prepare_observables();
        /*********************************************************************/
        if (_dtau <= 0.0){
            fprintf(stderr, "%s. %s: %f\n",
                "Cannot have zero or negative step in progress variable tau",
                "Setting default value", 0.001);
            _dtau = 0.001;
        }
        // set up the upper bound on the number of steps ---
        // this is just a fuse to halt paths that most likely won't converge
        _maxiter = (size_t)(ceil(_path.euc_sep.value * 100 / _dtau));
    }
    GeodesicSimulation::~GeodesicSimulation(){}
    void GeodesicSimulation::_prepare_outstream(std::ofstream &stream,
        std::string fout, bool truncate){
        Simulation::_prepare_outstream(stream, fout, truncate);
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
                bool output_header = false;
                stream << _path.header_str() << std::endl;
                // header output as a part of the path already
                _initial.write(stream, output_header);
            }
        }
    }
    void GeodesicSimulation::_write_status(){
        fprintf(stdout, "%s\n", "STATUS:");
        fprintf(stdout, "%s", _obs.status_string(VERBOSE).c_str());
    }
    void GeodesicSimulation::_write_config(std::ofstream& writeout){
        // open file in truncate mode,
        // write out state,
        // write out accumulators
        bool output_header = false;
        if (writeout.is_open()) {
            _path.current_tail().write(writeout, output_header);
        }
        else {
            perror("open");
        }
    }
    void GeodesicSimulation::compute_path(){
        size_t icalc = _icalc;
        size_t iblock = _iblock;
        size_t iprint = _iprint;
        size_t isave = _isave;
        size_t itape = _itape;
        //size_t itape = _itape;
        // _i<interval> == 0 => don't perform the action, but compute_path()
        // will use modular arithmetic: if _step % i<interval> == 0
        if (icalc == 0) icalc = _maxiter + 1;
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
        if (iprint == 0) iprint = _maxiter + 1;
        if (isave == 0) isave = _maxiter + 1;
        if (itape == 0) itape = _maxiter + 1;
        bool calc;
        bool escaped;
        size_t escape_iter;
        double delta_pe = pow(10.0, -8.0);
        _calcstep = 0;
        _blockstep = 0;
        geodesic::Record rec = _path.current_tail();
        // output initial record to path file with header with overwrite
        std::ofstream pfstream;
        bool output_header = true;
        if (_isave != 0) {
            _prepare_pfstream(pfstream);
            //rec.write(pfstream, output_header);
        }
        while( _comp->move(rec, _path.final(), _dtau) && (_step < _maxiter) ){
            _comp->update_PE(rec);
            if (rec.pe() > _el + delta_pe){
                fprintf(stdout, "%s\n", "Forbidden area encountered!");
                // TODO
                escaped = false;
                escape_iter = 0;
                while(!escaped && escape_iter < _max_escape_iter){
                    ++escape_iter;
                    //_comp->escape(rec, _path.final(), param);
                    // FIXME`
                    escaped = true;
                    //_update_record_from_links(rec);
                    //_comp->update_PE(rec);
                    //if (rec.pe() <= _el) escaped = true;
                }
                fprintf(stdout, "%s\n", "Forbidden area escaped!");
            }
            _path.append(rec);
            // take care of observables and bookkeeping
            _step = _step + 1;
            //if (_step % iprint == 0) {
            //    fprintf(stderr, "%s\n", "********************************");
            //    fprintf(stderr, "%s %zu\n", "STEP", _step);
            //    fprintf(stderr, "%s\n", "********************************");
            //}
            _step % icalc == 0? calc = true : calc = false;
            if (calc) {
                ++_calcstep;
                _pe.value = rec.pe();                 /* update observable */
                // temporary observables for checking algorithm implementation:
                // manage observable calculation and averaging
                _calculate_observables(iblock,   _path.current_tail().atom_state());
                // printing updates and saving a configuration will only happen
                // if observables are being calculated
                if (_step % iprint == 0) _write_status();
                if (_step % isave == 0) _write_config(pfstream);
            }
            // TODO add tape file?
            //if (_step % itape == 0) _write_config(pfstream);
            if (calc){
                _record_observables();
                if ((_should_write_data) and (_calcstep % iblock == 0)) {
                    _write_data();
                }
            }
            rec = _path.current_tail();
        }
        bool converged;
        if (_step >= _maxiter){
            converged = false;
        }
        else {
            converged = true;
        }
        fprintf(stdout, "%d\n", (int)converged);
        // write last status once run is complete
        //_write_status();
        // compute how close each atom is to its final position
        simple::AtomPolymer cur = rec.atom_state().polymers.at(0);
        simple::AtomPolymer fin = _path.final().atom_state().polymers.at(0);
        //fprintf(stdout, "%zu\n", _step);
        double diffsq = 0.0;
        for (size_t i = 0; i < BasePolymer::nb() + 1; ++i){
            diffsq += normsq(subtract(fin.atoms.at(i).position,
                                 cur.atoms.at(i).position));
        }
        double diff = sqrt(diffsq);
        fprintf(stdout, "%s: %5.7f\n",
            "Config space distance to the endpoint",
            diff);
        if (converged){
            //fprintf(stdout, "%s\n",
            //    "Will add that distance to the path length");
            _path.length.value += diff;
        }
        fprintf(stdout, "%5.7f\n", _path.length.value);
        // output last data results if they haven't been output just now
        _record_observables();
        if ((_should_write_data) and (_calcstep % iblock != 0)) {
            _write_data();
        }
        if (_isave != 0){
            // output endpoint record before closing
            output_header = false;
            _path.final().write(pfstream, output_header);
            pfstream.close();
        }
        _obs.run_end();
    }
    SLERPGeodesicSimulation::SLERPGeodesicSimulation(std::string name,
        std::string initial,
        std::string final,
        std::string cndir,
        std::string tpdir,
        std::string dtdir,
        ObservableContainer& container,
        geodesic::SLERP* comp,
        double landscape_energy,
        bool should_write_data,
        double dtau,
        size_t icalc,               // calculate every 10 steps
        size_t iblock,              // average every 100 calcsteps
        size_t iprint,              // print status every 1000 steps
        size_t isave,               // save config + obs every 1000 steps
        size_t itape,               // save config every 100 steps
        size_t maxiter,             // max propagation steps
        size_t max_escape_iter) :
        GeodesicSimulation(name, initial, final, cndir, tpdir, dtdir,
            container, comp, landscape_energy, should_write_data,
            dtau, icalc, iblock, iprint, isave, itape,
            maxiter, max_escape_iter){}
    SLERPGeodesicSimulation::~SLERPGeodesicSimulation(){}
    ShortStepGeodesicSimulation::ShortStepGeodesicSimulation(std::string name,
        std::string initial,
        std::string final,
        std::string cndir,
        std::string tpdir,
        std::string dtdir,
        ObservableContainer& container,
        geodesic::ShortStep* comp,
        double landscape_energy,
        bool should_write_data,
        double dtau,
        size_t icalc,               // calculate every 10 steps
        size_t iblock,              // average every 100 calcsteps
        size_t iprint,              // print status every 1000 steps
        size_t isave,               // save config + obs every 1000 steps
        size_t itape,               // save config every 100 steps
        size_t maxiter,             // max propagation steps
        size_t max_escape_iter) :
        GeodesicSimulation(name, initial, final, cndir, tpdir, dtdir,
            container, comp, landscape_energy, should_write_data,
            dtau, icalc, iblock, iprint, isave, itape,
            maxiter, max_escape_iter){}
    ShortStepGeodesicSimulation::~ShortStepGeodesicSimulation(){}
    ShoveGeodesicSimulation::ShoveGeodesicSimulation(std::string name,
        std::string initial,
        std::string final,
        std::string cndir,
        std::string tpdir,
        std::string dtdir,
        ObservableContainer& container,
        geodesic::SHOVE* comp,
        double landscape_energy,
        bool should_write_data,
        double sigma,               // max step in configuration space
        size_t icalc,               // calculate every 10 steps
        size_t iblock,              // average every 100 calcsteps
        size_t iprint,              // print status every 1000 steps
        size_t isave,               // save config + obs every 1000 steps
        size_t itape,               // save config every 100 steps
        size_t maxiter,             // max propagation steps
        size_t max_escape_iter) :   // max escape steps
    GeodesicSimulation(name, initial, final, cndir, tpdir, dtdir,
        container, comp, landscape_energy, should_write_data,
        sigma, icalc, iblock, iprint, isave, itape,
        maxiter, max_escape_iter){}
    ShoveGeodesicSimulation::~ShoveGeodesicSimulation(){}
    PlerpGeodesicSimulation::PlerpGeodesicSimulation(std::string name,
        std::string initial,
        std::string final,
        std::string cndir,
        std::string tpdir,
        std::string dtdir,
        ObservableContainer& container,
        geodesic::PLERP* comp,
        double landscape_energy,
        bool should_write_data,
        double sigma,               // max step in configuration space
        size_t icalc,               // calculate every 10 steps
        size_t iblock,              // average every 100 calcsteps
        size_t iprint,              // print status every 1000 steps
        size_t isave,               // save config + obs every 1000 steps
        size_t itape,               // save config every 100 steps
        size_t maxiter,             // max propagation steps
        size_t max_escape_iter) :   // max escape steps
    GeodesicSimulation(name, initial, final, cndir, tpdir, dtdir,
        container, comp, landscape_energy, should_write_data,
        sigma, icalc, iblock, iprint, isave, itape,
        maxiter, max_escape_iter){}
    PlerpGeodesicSimulation::~PlerpGeodesicSimulation(){}
} // namespace simple
