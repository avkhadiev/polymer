// 2017 Artur Avkhadiev
/*! \file simple_simulation.cpp
*/
#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <sstream>
#include <fstream>
#include "../include/simple_simulation.h"
//namespace parameters{
//    int NM = 1;
//    int NB = 1;
//    double M = 1.0;
//    double D = 3.0;
//}
//namespace obs{
//    ScalarObservable inst_K
//        = declare_scalar_observable("Instanteneous Kinetic Energy",
//            "\\epsilon",
//            "K");
//    ScalarObservable avg_K
//        = declare_scalar_observable("Average Kinetic Energy",
//            "\\epsilon",
//            "\\bar{K}");
//    ScalarObservable inst_V
//        = declare_scalar_observable("Instantaneous Potential Energy",
//            "\\epsilon",
//            "V_{LJ}");
//    ScalarObservable avg_V
//        = declare_scalar_observable("Average Potential Energy",
//            "\\epsilon",
//            "\\bar{V}_{LJ}");
//    ScalarObservable W
//        = declare_scalar_observable("Negative Virial",
//            "\\epsilon",
//            "-W_{LJ}");
//    ScalarObservable WC
//        = declare_scalar_observable("Negative Constraint Virial",
//            "\\epsilon",
//            "-W_{C}");
//} // namespace observables
//namespace rattle {
//    double tol = pow(10, -3);
//    double tiny = pow(10, -7);
//    int maxiter = pow(10, 3);
//}
namespace simple{
    Simulation::Simulation(std::string name,
        std::string cndir,
        std::string tpdir,
        std::string dtdir,
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
        //*********************************************************************
        fprintf(stdout, "%s\n", "Simulation is set up:");
        fprintf(stdout, "%s: %s\n", "Name", _name.c_str());
        fprintf(stdout, "%s: %d\n",
            "Number of Molecules",
            _cfg.nmolecules());
        fprintf(stdout, "%s: %d\n",
            "Number of Bonds in a Molecule",
            _cfg.polymer_nb());
        fprintf(stdout, "%s: %s/%s\n",
            "Config File",
            _cndir.c_str(),
            _cnfname.c_str());
        fprintf(stdout, "%s: %s/%s\n",
            "Tape File",
            _tpdir.c_str(),
            _tpfname.c_str());
        fprintf(stdout, "%s: %f\n", "Timestep", _dt);
        fprintf(stdout, "%s: %lu steps\n", "Calculation Interval", _icalc);
        fprintf(stdout, "%s: %lu steps\n", "Status Print Interval", _iprint);
        fprintf(stdout, "%s: %lu steps\n", "Observable Print Interval", _idata);
        fprintf(stdout, "%s: %lu steps\n", "Config Update Interval", _isave);
        fprintf(stdout, "%s: %lu steps\n", "Tape Update Interval", _itape);
        _write_status();
    }
    Simulation::~Simulation(){}
    void Simulation::_write_status(){
        fprintf(stdout, "%s\n", "STATUS:");
        fprintf(stdout, "%s: %f\n", "Time", _cfg.atom_state().time());
        fprintf(stdout, "%s\n", _cfg.atom_state().to_string(true, false).c_str());
        fprintf(stdout, "%s\n", _write_accumulators().c_str());
    }
    std::string Simulation::_write_accumulators(){
        return "";
    }
    void Simulation::_read_accumulators(std::ifstream& input_stream){
    }
    void Simulation::_read_config(){
        std::string fin;              /*>> stores path to input file */
        fin = _cndir + _cnfname;      /*>> stores path to input directory */
        std::ifstream readout;
        readout.open(fin, std::ifstream::in);
        if (!readout.is_open()) {
            std::string err_msg = "read_config: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
            perror("open");
        }
        else{
            _cfg.read_atom_state(readout);
            _read_accumulators(readout);
        }
        readout.close();
    }
    void Simulation::_write_config(){
        // open file in truncate mode,
        // write out state,
        // write out accumulators
        bool verbose = true;
        bool output_header = true;
        std::string fout;
        fout = _cndir + _cnfname;
        std::ofstream writeout;
        // open writeout for output operations in truncate mode
        writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
        // now file has to be opened
        if (writeout.is_open()) {
            std::string s_str;
            writeout << _cfg.atom_state().to_string(verbose,
                output_header);
            writeout << _write_accumulators();
        }
        else {
            // if file still could not be opened
            std::string err_msg = "write_state_to_file: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
            perror("open");
        }
        writeout.close();
    }
    void Simulation::_write_tape(){
        bool verbose = true;
        bool output_header = false;
        std::string fout;
        fout = _tpdir + _tpfname;
        std::ofstream writeout;
        // open writeout for output operations in append mode
        writeout.open(fout, std::ofstream::out | std::ofstream::app);
        if (!writeout.is_open()) {
            // if file does not exist, open in truncate mode
            // header is written only in the beginning of tape file
            output_header = true;
            writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
        }
        // now file has to be opened
        if (writeout.is_open()) {
            writeout << _cfg.atom_state().to_string(verbose,
                output_header);
        }
        else {
            // if file still could not be opened
            std::string err_msg = "write_state_to_file: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
            perror("open");
        }
        writeout.close();
    }
    void Simulation::_write_data(){
        bool overwrite = false;
        _obs.writeout(_tpdir, _name, overwrite);
        // flush all records to save memory and not output them next time
        _obs.clear();
    }
    void Simulation::_calculate(){
        _calcstep += 1;
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
        if (idata < icalc) idata = icalc;
        // can't output data more often that new data is calculated
        // shouldn't print to tape less often that configuration is output
        if (idata < icalc) idata = icalc;
        if (isave < itape) isave = itape;
        bool calc;
        while(_step < nsteps){
            _step = _step + 1;
            _step % icalc == 0? calc = true : calc = false;
            _int.move(_dt, _cfg.atom_state(), calc);
            if (calc) _calculate();
            if (_step % iprint == 0) _write_status();
            if (_step % itape == 0) _write_tape();
            if (_step % idata == 0) _write_data();
            if (_step % isave == 0) _write_config();
            if (_step % itape == 0) _write_tape();
        }
        _write_config();
        _write_data();
    }
} // namespace simple
