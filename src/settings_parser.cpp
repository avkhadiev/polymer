// 2017 Artur Avkhadiev
/*! \file settings_parser.cpp
*/
#include <sys/stat.h>                       /**> checks for file existence */
#include <limits>                           /**> controls output precision */
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include "../include/settings_parser.h"
namespace default_settings {
    const std::string potential_header =
        "# potential: epp ess eps ss sp rc_ss rc_sp";
    const std::string state_config_header
        = "# state configuration: mp ms np nb nc d/sp";
    const std::string init_header = "# initialization: rho_s kbT/e pol_conformation pol_energy_per_atom";
    const std::string integration_header = "# integration: runtime dt tol";
    const std::string geodesic_header
        = "# geodesic: geodir_input geodir_output geodir_data el dtau save_md_path";
    const std::string io_header = "i/o: cndir datadir tpdir icalc iprint isave iblock itape should_write_data";
    // POTENTIAL SETTINGS
    const double epp = 1.0;
    const double ess = 1.0;
    const double eps = 1.0;
    const double sp = 1.0;
    const double ss = 1.0;
    const double rc_ss = 5.0;
    const double rc_ps = rc_ss;
    // STATE CONFIGURATION SETTINGS
    const double mp = 1.0;
    const double ms = 1.0;
    const int np = 1;
    const int nb = 10;
    const int nc = 4;
    const double d = 0.86;
    // INITIALIZATION SETTINGS
    const double rho_s = 0.7;
    const double temperature = 1.0;
    const int polymer_planar_conformation = 0;
    // INTEGRATION SETTINGS
    const double runtime = 10.0;
    const double dt = pow(10, -4.0);
    const double tol = pow(dt, 2.0);
    // const double tiny = pow(10, -8.0);
    // const int maxiter = pow(10, 4.0);
    // GEODESIC SETTINGS
    const std::string geodir_input
        = "/Users/Arthur/stratt/polymer/test/geodesics/";
    const std::string geodir_output
        = "/Users/Arthur/stratt/polymer/test/geodesics/";
    const std::string geodir_data
        = "/Users/Arthur/stratt/polymer/test/geodesics/";
    const double el = 0.0;
    const double dtau = 0.001;
    const int save_md_path = 0;
    // I/O SETTINGS
    const std::string cndir = "/Users/Arthur/stratt/polymer/test/";
    const std::string dtdir = "/Users/Arthur/stratt/polymer/test/";
    const std::string tpdir = "/Users/Arthur/stratt/polymer/test/";
    const bool should_write_data = true;
    const size_t icalc = 10;
    const size_t iprint = 10;
    const size_t isave = 100;
    const size_t iblock = 100;
    const size_t itape = 10;
} // namespace default_settings
SettingsParser::SettingsParser() :
    // HEADER STRINGS FOR OUTPUT
    potential_header(default_settings::potential_header),
    state_config_header(default_settings::state_config_header),
    init_header(default_settings::init_header),
    integration_header(default_settings::integration_header),
    geodesic_header(default_settings::geodesic_header),
    io_header(default_settings::io_header),
    // POTENTIAL SETTINGS
    epp(default_settings::epp),
    ess(default_settings::ess),
    eps(default_settings::eps),
    sp(default_settings::sp),
    ss(default_settings::ss),
    rc_ss(default_settings::rc_ss),
    rc_ps(default_settings::rc_ps),
    // STATE CONFIGURATION SETTINGS
    mp(default_settings::mp),
    ms(default_settings::ms),
    np(default_settings::np),
    nb(default_settings::nb),
    nc(default_settings::nc),
    d(default_settings::d),
    // INITIALIZATION SETTINGS
    rho_s(default_settings::rho_s),
    temperature(default_settings::temperature),
    polymer_planar_conformation(default_settings::polymer_planar_conformation),
    // INTEGRATION SETTINGS
    runtime(default_settings::runtime),
    dt(default_settings::dt),
    tol(default_settings::tol),
    // GEODESOC SETTINGS
    geodir_input(default_settings::geodir_input),
    geodir_output(default_settings::geodir_output),
    geodir_data(default_settings::geodir_data),
    el(default_settings::el),
    dtau(default_settings::dtau),
    save_md_path(default_settings::save_md_path),
    // I/O SETTINGS
    cndir(default_settings::cndir),
    dtdir(default_settings::dtdir),
    tpdir(default_settings::tpdir),
    should_write_data(default_settings::should_write_data),
    icalc(default_settings::icalc),
    iprint(default_settings::iprint),
    isave(default_settings::isave),
    iblock(default_settings::iblock),
    itape(default_settings::itape){}
SettingsParser::~SettingsParser(){}
std::vector<std::string> SettingsParser::get_args(std::string line){
    std::istringstream ss(line.c_str());
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> args(begin, end);
    return args;
}
SettingsParser::SettingsParser(std::string fin) :
    SettingsParser()
{
    read(fin);
}
void SettingsParser::read_potential(std::ifstream& stream){
    std::string line;
    std::getline(stream, line);     // comment string
    std::getline(stream, line);
    std::vector<std::string> args = get_args(line);
    // convert strings to data
    epp = atof(args.at(0).c_str());
    ess = atof(args.at(1).c_str());
    eps = atof(args.at(2).c_str());
    ss = atof(args.at(3).c_str());
    sp = atof(args.at(4).c_str());
    rc_ss = atof(args.at(5).c_str());
    rc_ps = atof(args.at(6).c_str());
}
void SettingsParser::read_state_config(std::ifstream& stream){
    std::string line;
    std::getline(stream, line);     // comment string
    std::getline(stream, line);
    std::vector<std::string> args = get_args(line);
    mp = atof(args.at(0).c_str());
    ms = atof(args.at(1).c_str());
    np = atof(args.at(2).c_str());
    nb = atoi(args.at(3).c_str());
    nc = atoi(args.at(4).c_str());
    d = atof(args.at(5).c_str());
}
void SettingsParser::read_init(std::ifstream& stream){
    std::string line;
    std::getline(stream, line);     // comment string
    std::getline(stream, line);
    std::vector<std::string> args = get_args(line);
    rho_s = atof(args.at(0).c_str());
    temperature = atof(args.at(1).c_str());
    polymer_planar_conformation = atoi(args.at(2).c_str());
    polymer_energy = atof(args.at(3).c_str());
}
void SettingsParser::read_integration(std::ifstream& stream){
    std::string line;
    std::getline(stream, line);     // comment string
    std::getline(stream, line);
    std::vector<std::string> args = get_args(line);
    runtime = atof(args.at(0).c_str());
    dt = atof(args.at(1).c_str());
    tol = atof(args.at(2).c_str());
}
void SettingsParser::read_geodesic(std::ifstream& stream){
    std::string line;
    std::getline(stream, line);     // comment string
    std::getline(stream, line);     // geodesic directory for inputs
    geodir_input = line;
    std::getline(stream, line);     // geodesic directory for outputs
    geodir_output = line;
    std::getline(stream, line);     // geodesic directory for observables
    geodir_data = line;
    std::getline(stream, line);     // landscape energy and ``affine parameter''
    std::vector<std::string> args = get_args(line);
    el = atof(args.at(0).c_str());
    dtau = atof(args.at(1).c_str());
    save_md_path = atoi(args.at(2).c_str());
}
void SettingsParser::read_io(std::ifstream& stream){
    std::string line;
    std::getline(stream, line);    // header lineg
    // simulation state, directories, data, and state configuration directories
    std::getline(stream, line);
    cndir = line;
    std::getline(stream, line);
    dtdir = line;
    std::getline(stream, line);
    tpdir = line;
    std::getline(stream, line);
    // intervals for calculating observables,
    // printing out simulation status,
    // writing simulation status to file,
    // writing observables to file,
    // writing stat configuration to a "tape" file
    int should_write_int;
    sscanf(line.c_str(), "%zu %zu %zu %zu %zu %d",
        &icalc, &iprint, &isave, &iblock, &itape, &should_write_int);
    should_write_data = should_write_int;
}
void SettingsParser::read(std::string fin){
    std::ifstream readout;
    readout.open(fin, std::ifstream::in);
    if (!readout.is_open()) {
        std::string err_msg = "read_config: unable to open file at";
        fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
        perror("open");
    }
    else{
        read_potential(readout);
        read_state_config(readout);
        read_init(readout);
        read_integration(readout);
        read_geodesic(readout);
        read_io(readout);
    }
    readout.close();
}
void SettingsParser::write_potential(std::ofstream& stream) const{
    if (stream.is_open()) {
        // output header
        stream << potential_header << std::endl;
        // output values
        stream << std::to_string(epp) << " ";
        stream << std::to_string(ess) << " ";
        stream << std::to_string(eps) << " ";
        stream << std::to_string(sp) << " ";
        stream << std::to_string(ss) << " ";
        stream << std::to_string(rc_ss) << " ";
        stream << std::to_string(rc_ps) << std::endl;
    }
    else {
        // if file still could not be opened
        std::string err_msg = "writeout: unable to open file";
        fprintf(stderr, "%s\n", err_msg.c_str());
        perror("open");
    }
}
void SettingsParser::write_state_config(std::ofstream& stream) const{
    if (stream.is_open()) {
        // output header
        stream << state_config_header << std::endl;
        // output values
        stream << std::to_string(mp) << " ";
        stream << std::to_string(ms) << " ";
        stream << std::to_string(np) << " ";
        stream << std::to_string(nb) << " ";
        stream << std::to_string(nc) << " ";
        stream << std::to_string(d) << std::endl;
    }
    else {
        // if file still could not be opened
        std::string err_msg = "writeout: unable to open file";
        fprintf(stderr, "%s\n", err_msg.c_str());
        perror("open");
    }
}
void SettingsParser::write_init(std::ofstream& stream) const{
    if (stream.is_open()) {
        // output header
        stream << init_header << std::endl;
        // output values
        stream << std::to_string(rho_s) << " ";
        stream << std::to_string(temperature) << " ";
        stream << std::to_string(polymer_planar_conformation) << " ";
        stream << std::to_string(polymer_energy) << std::endl;
    }
    else {
        // if file still could not be opened
        std::string err_msg = "writeout: unable to open file";
        fprintf(stderr, "%s\n", err_msg.c_str());
        perror("open");
    }
}
void SettingsParser::write_integration(std::ofstream& stream) const{
    if (stream.is_open()) {
        // output header
        stream << integration_header << std::endl;
        // output values
        stream << std::to_string(runtime) << " ";
        stream << std::to_string(dt) << " ";
        stream << std::to_string(tol) << " ";
        stream << std::to_string(save_md_path) << std::endl;
    }
    else {
        // if file still could not be opened
        std::string err_msg = "writeout: unable to open file";
        fprintf(stderr, "%s\n", err_msg.c_str());
        perror("open");
    }
}
void SettingsParser::write_geodesic(std::ofstream& stream) const{
    if (stream.is_open()) {
        // output header
        stream << geodesic_header << std::endl;
        // output values
        stream << geodir_input << " ";
        stream << geodir_output << " ";
        stream << geodir_data << " ";
        stream << std::to_string(el) << " ";
        stream << std::to_string(dtau) << std::endl;
    }
    else{
        // if file still could not be opened
        std::string err_msg = "writeout: unable to open file";
        fprintf(stderr, "%s\n", err_msg.c_str());
        perror("open");
    }
}
void SettingsParser::write_io(std::ofstream& stream) const{
    if (stream.is_open()) {
        // output header
        stream << io_header << std::endl;
        // output values
        stream << cndir << std::endl;
        stream << dtdir << std::endl;
        stream << tpdir << std::endl;
        stream << std::to_string(icalc) << " ";
        stream << std::to_string(iprint) << " ";
        stream << std::to_string(isave) << " ";
        stream << std::to_string(iblock) << " ";
        stream << std::to_string(itape) << " ";
        stream << std::to_string((int)should_write_data) << std::endl;
    }
    else {
        // if file still could not be opened
        std::string err_msg = "writeout: unable to open file";
        fprintf(stderr, "%s\n", err_msg.c_str());
        perror("open");
    }
}
void SettingsParser::write(std::string outdir, std::string sim_name) const{
    std::ofstream writeout;
    std::string fout = outdir + sim_name + ".cfg";
    // open writeout for output operations in truncate mode
    writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
    if (writeout.is_open()) {
        // if file could be opened...
        write_potential(writeout);
        write_state_config(writeout);
        write_init(writeout);
        write_integration(writeout);
        write_geodesic(writeout);
        write_io(writeout);
    }
    else {
        // if file still could not be opened
        std::string err_msg = "write_state_to_file: unable to open file at";
        fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
        perror("open");
    }
    writeout.close();
}
