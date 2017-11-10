// 2017 Artur Avkhadiev
/*! \file simple_simulation.cpp
*/
#include <vector>
#include <sys/stat.h>                       /**> checks for file existence */
#include <string>
#include <stdexcept>
#include <map>
#include <sstream>
#include <fstream>
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
        double dt,
        size_t icalc,
        size_t iprint,
        size_t isave,
        size_t idata,
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
        _dt(dt),
        _step(0),
        _calcstep(0),
        _icalc(icalc),
        _iprint(iprint),
        _isave(isave),
        _idata(idata),
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
        fprintf(stdout, "%s\n", "Simulation is set up:");
        fprintf(stdout, "%s: %s\n", "Name", _name.c_str());
        fprintf(stdout, "%s: %d\n",
            "Number of Molecules",
            _cfg.nmolecules());
        fprintf(stdout, "%s: %d\n",
            "Number of Bonds in a Molecule",
            _cfg.polymer_nb());
        if (_is_input_given){
            fprintf(stdout, "%s: %s\n",
                "Input Config File",
                _infile.c_str());
            _read_config();                 // may change state & observables
        }
        fprintf(stdout, "%s: %s%s\n",
            "Config File",
            _cndir.c_str(),
            _cnfname.c_str());
        fprintf(stdout, "%s: %s%s\n",
            "Tape File",
            _tpdir.c_str(),
            _tpfname.c_str());
        fprintf(stdout, "%s: %f\n", "Timestep", _dt);
        fprintf(stdout, "%s: %lu steps\n", "Calculation Interval", _icalc);
        fprintf(stdout, "%s: %lu steps\n", "Status Print Interval", _iprint);
        fprintf(stdout, "%s: %lu steps\n", "Observable Print Interval", _idata);
        fprintf(stdout, "%s: %lu steps\n", "Config Update Interval", _isave);
        fprintf(stdout, "%s: %lu steps\n", "Tape Update Interval", _itape);
    }
    Simulation::~Simulation(){}
    void Simulation::_write_status(){
        fprintf(stdout, "%s\n", "STATUS:");
        fprintf(stdout, "%s: %f\n", "Time", _cfg.atom_state().time());
        fprintf(stdout, "%s", _obs.status_string().c_str());
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
            // get steps
            std::string line;
            getline(readout, line);
            sscanf(line.c_str(), "%zu", &_step);
            getline(readout, line);
            sscanf(line.c_str(), "%zu", &_calcstep);
            // get observables
            _obs.read_status(readout);
        }
        readout.close();
    }
    void Simulation::_write_config(){
        // open file in truncate mode,
        // write out state,
        // write out accumulators
        bool verbose = false;
        bool output_header = true;
        std::string fout;
        fout = _cndir + _cnfname;
        std::ofstream writeout;
        // open writeout for output operations in truncate mode
        writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
        // now file has to be opened
        if (writeout.is_open()) {
            std::string s_str;
            writeout << _cfg.bond_state().to_string(verbose, output_header);
            writeout << _step << std::endl;
            writeout << _calcstep << std::endl;
            writeout << _obs.config_string();
        }
        else {
            // if file still could not be opened
            std::string err_msg = "write_state_to_file: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
            perror("open");
        }
        writeout.close();
    }
    void Simulation::_prepare_tpstream(std::ofstream &tpstream){
        bool verbose = false;
        std::string fout = _tpdir + _tpfname;
        // check for file existence
        // if doesn't exist, create new one and output header information
        struct stat buffer;
        if (stat(fout.c_str(), &buffer) == 0) {
            tpstream.open(fout, std::ofstream::out | std::ofstream::app);
        }
        else{
            fprintf(stderr, "%s\n",
                "tape file doesn't exist, creating new one.");
            tpstream.open(fout, std::ofstream::out | std::ofstream::trunc);
            tpstream << _cfg.atom_state().header_str(verbose) << std::endl;
        }
        // now file has to be opened
        if (!tpstream.is_open()) {
            // if file still could not be opened
            std::string err_msg = "write_state_to_file: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
            perror("open");
        }
        else{
            fprintf(stdout, "%s\n", "Tape file opened for writing.");
        }
    }
    void Simulation::_write_tape(std::ofstream& tpstream){
        bool verbose = false;
        bool output_header = false;
        // now file has to be opened
        if (tpstream.is_open()) {
            tpstream << _cfg.bond_state().to_string(verbose,
                output_header);
        }
        else {
            perror("open");
        }
    }
    void Simulation::_write_data(){
        bool overwrite = false;
        _obs.write_data(_tpdir, _name, overwrite);
    }
    void Simulation::_calculate(){
        // calc step has to be updated manually
        _obs.update(_cfg.atom_state(), _calcstep);
    }
    void Simulation::evolve(double runtime){
        size_t nsteps = _step + (size_t)(runtime / _dt);
        size_t icalc = _icalc;
        size_t iprint = _iprint;
        size_t isave = _isave;
        size_t idata = _idata;
        size_t itape = _itape;
        // _i<interval> == 0 => don't perform the action, but evolve() will
        // use modular arithmetic: if _step % i<interval> == 0
        if (icalc == 0) icalc = nsteps + 1;
        if (iprint == 0) iprint = nsteps + 1;
        if (isave == 0) isave = nsteps + 1;
        if (idata == 0) idata = nsteps + 1;
        if (itape == 0) itape = nsteps + 1;
        // can't output data more often that new data is calculated
        if (idata < icalc) idata = icalc;
        bool calc;
        // write out initial state data as necessary
        std::ofstream tpstream;
        if (_itape != 0){
            _prepare_tpstream(tpstream);
            _write_tape(tpstream);
        }
        if (_idata != 0){
            bool overwrite = true;
            _obs.write_data(_dtdir, _name, overwrite);
        }
        // if (_iprint != 0) _write_status();
        while(_step < nsteps){
            _step = _step + 1;
            _step % icalc == 0? calc = true : calc = false;
            _int.move(_dt, _cfg.atom_state(), calc);
            if (calc) {
                ++_calcstep;
                _calculate();
            }
            if (_step % iprint == 0) _write_status();
            if (_step % itape == 0) _write_tape(tpstream);
            if (_step % idata == 0) _write_data();
            if (_step % isave == 0) _write_config();
        }
        // write to tape if this state has not been already written
        // and tape interval is non-zero
        if (_step % itape != 0 && itape != nsteps + 1){
            _write_tape(tpstream);
            tpstream.close();
            fprintf(stdout, "%s\n", "Tape file closed.");
        }
        if (idata != nsteps + 1) _write_data();
    }
} // namespace simple
