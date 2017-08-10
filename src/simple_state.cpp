// 2017 Artur Avkhadiev
/*! \file simple_state.cpp
*/
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "../include/parsing.h"
#include "../include/simple_polymer.h"
#include "../include/simple_state.h"
#define DEFAULT_NM 1
namespace simple {
    /***************************************************************************
    *                               BASE STATE
    ***************************************************************************/
    int BaseState::_nm = DEFAULT_NM;
    BaseState::BaseState(double time) :
        _time (time){}
    BaseState::~BaseState(){}
    std::string BaseState::_header_str(bool verbose) const {
        std::string s_str_padding = "*****************************************";
        std::string s_str_header = "SIMPLE STATE:";
        std::string s_str_nm_header = "NM = ";
        std::string s_str_nm = std::to_string(nm());
        std::string s_str_nb_header = "NB = ";
        std::string s_str_nb = std::to_string(BasePolymer::nb());
        std::string s_str_m_header = "M = ";
        std::string s_str_m = std::to_string(BasePolymer::m());
        std::string s_str_d_header = "D = ";
        std::string s_str_d = std::to_string(BasePolymer::d());
        std::string s_str;
        if(verbose){
            s_str = s_str_padding + "\n"
                + s_str_header + "\n"
                + s_str_nm_header + s_str_nm + " "
                + s_str_nb_header + s_str_nb + " "
                + s_str_m_header + s_str_m + " "
                + s_str_d_header + s_str_d + " " + "\n"
                + s_str_padding;
        }
        else{
            s_str = s_str_nm + " "
                + s_str_nb + " "
                + s_str_m + " "
                + s_str_d;
        }
        return s_str;
    }
    std::string BaseState::_time_str(bool verbose) const{
        std::string s_str_t_header = "T = ";
        std::string s_str_t = std::to_string(_time);
        std::string s_str;
        if(verbose){
            s_str = s_str_t_header + s_str_t;
        }
        else{
            s_str = s_str_t;
        }
        return s_str;
    }
    void BaseState::read_header(std::string non_verbose_header){
        // in non-verbose lines, the information about a header is one line:
        // Number of Molecules,
        // Number of Bonds in Each Molecule
        // Mass of Atoms in Each Molecule
        // Interatomic Bond Length in Each Molecule
        std::istringstream ss(non_verbose_header.c_str());
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> words(begin, end);
        // convert strings to data
        _nm = atoi(words.at(0).c_str());
        BasePolymer::set_nb(atoi(words.at(1).c_str()));
        BasePolymer::set_m(atof(words.at(2).c_str()));
        BasePolymer::set_d(atof(words.at(3).c_str()));
    }
    void BaseState::read_time(std::string non_verbose_time){
        // in non-verbose lines, the information about time is a single number
        std::istringstream ss(non_verbose_time.c_str());
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> words(begin, end);
        // convert strings to data
        set_time(atof(words.at(0).c_str()));
    }
    std::string BaseState::_to_string_helper(std::string molecules_str,
        bool verbose, bool output_header) const{
        std::string s_str = "";
        if(output_header){
            s_str += _header_str(verbose) + "\n"
                + _time_str(verbose) + "\n"
                + molecules_str;
        }
        else{
            s_str += _time_str(verbose) + "\n" + molecules_str;
        }
        return s_str;
    }
    void BaseState::_write_to_file_helper(std::string molecules_str,
        std::string outdir,
        std::string fname,
        bool verbose,
        bool overwrite) const{
        // stores path to output file
        std::string fout;
        // stores path to output directoy
        fout = outdir + parse_string(fname) + ".cfg";
        std::ofstream writeout;
        // open writeout for output operations and
        // set the stream's position indicator to the end of the stream before each output operation.
        if (overwrite) {
            // if overwrite is allowed, try to open in truncate mode
            writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
        }
        if (!overwrite) {
            // if overwrite is forbidden, try to open the file in append mode
            writeout.open(fout, std::ofstream::out | std::ofstream::app);
            if (!writeout.is_open()) {
            // if file does not exist, open in truncate mode
                writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
            }
        }
        // now file has to be opened
        if (writeout.is_open()) {
            // if file could be opened...
            std::string s_str;
            s_str = _to_string_helper(molecules_str, verbose, overwrite);
            writeout << s_str;
        }
        else {
            // if file still could not be opened
            std::string err_msg = "write_state_to_file: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
            perror("open");
        }
        writeout.close();
    }
    /***************************************************************************
    *                               BOND STATE
    ***************************************************************************/
    BondState::BondState(const std::vector<BondPolymer>& new_polymers,
        double representation_time) :
        BaseState(representation_time),
        polymers(new_polymers){
        // default value is empty
        if(polymers.empty()){
            polymers.resize(nm());
        }
        // there is a fixed number of molecules in the state
        else if(polymers.size() != nm()){
            fprintf(stderr, "%s (%d); %s\n",
                "BondState: given vector of polymers does not match the required size",
                nm(),
                "the vector will be resized.");
            polymers.resize(nm());
        }
    }
    BondState::~BondState(){}
    bool BondState::operator==(const BondState &other) const{
        bool time = (_time == other.time());
        bool molecules = true;
        for(int i = 0; i < nm(); ++i){
            molecules = molecules && (polymers.at(i) == other.polymers.at(i));
        }
        return time && molecules;
    }
    bool BondState::operator!=(const BondState &other) const{
        return !(*this == other);
    }
    std::string BondState::_molecules_str(bool verbose) const{
        std::string s_str_molecules = "";
        std::string s_str_molecules_padding = "-------------------------------";
        for (int i = 0; i < nm(); ++i) {
            if (verbose) {
                // verbose output
                s_str_molecules += std::to_string(i)
                    + " "
                    + polymers.at(i).to_string(verbose)
                    + s_str_molecules_padding + "\n";
            }
            else {
                // non-verbose output
                s_str_molecules += polymers.at(i).to_string(verbose);
            }
        }
        return s_str_molecules;
    }
    std::string BondState::to_string(bool verbose, bool output_header) const{
        std::string molecules_str = _molecules_str(verbose);
        return _to_string_helper(molecules_str, verbose, output_header);
    }
    ::std::ostream& operator<<(::std::ostream& os, const BondState& s){
        return os << s.to_string().c_str();
    }
    void BondState::write_to_file(std::string outdir,
        std::string fname,
        bool verbose,
        bool overwrite){
        std::string molecules_str = _molecules_str(verbose);
        _write_to_file_helper(molecules_str, outdir, fname, verbose, overwrite);
    }
    /***************************************************************************
    *                               ATOM STATE
    ***************************************************************************/
    AtomState::AtomState(const std::vector<AtomPolymer>& new_polymers,
        double representation_time) :
        BaseState(representation_time),
        polymers(new_polymers){
        // default value is empty
        if(polymers.empty()){
            polymers.resize(nm());
        }
        // there is a fixed number of molecules in the state
        else if(polymers.size() != nm()){
            fprintf(stderr, "%s (%d); %s\n",
                "AtomState: given vector of polymers does not match the required size",
                nm(),
                "the vector will be resized.");
            polymers.resize(nm());
        }
    }
    AtomState::~AtomState(){}
    bool AtomState::operator==(const AtomState &other) const{
        bool time = (_time == other.time());
        bool molecules = true;
        for(int i = 0; i < nm(); ++i){
            molecules = molecules && (polymers.at(i) == other.polymers.at(i));
        }
        return time && molecules;
    }
    bool AtomState::operator!=(const AtomState &other) const{
        return !(*this == other);
    }
    std::string AtomState::_molecules_str(bool verbose) const{
        std::string s_str_molecules = "";
        std::string s_str_molecules_padding = "-------------------------------";
        for (int i = 0; i < nm(); ++i) {
            if (verbose) {
                // verbose output
                s_str_molecules += std::to_string(i)
                    + " "
                    + polymers.at(i).to_string(verbose)
                    + s_str_molecules_padding + "\n";
            }
            else {
                // non-verbose output
                s_str_molecules += polymers.at(i).to_string(verbose);
            }
        }
        return s_str_molecules;
    }
    std::string AtomState::to_string(bool verbose, bool output_header) const{
        std::string molecules_str = _molecules_str(verbose);
        return _to_string_helper(molecules_str, verbose, output_header);
    }
    ::std::ostream& operator<<(::std::ostream& os, const AtomState& s){
        return os << s.to_string().c_str();
    }
    void AtomState::write_to_file(std::string outdir,
        std::string fname,
        bool verbose,
        bool overwrite){
        std::string molecules_str = _molecules_str(verbose);
        _write_to_file_helper(molecules_str, outdir, fname, verbose, overwrite);
    }
    /***************************************************************************
    *                  ATOM & BOND STATE (conversion and input)
    ***************************************************************************/
    void BondState::update(const AtomState &atomic){
        _time = atomic.time();
        for(int i = 0; i < atomic.nm(); ++i){
            polymers.at(i) = BondPolymer(atomic.polymers.at(i));
        }
    }
    void AtomState::update(const BondState &bond){
        _time = bond.time();
        for(int i = 0; i < bond.nm(); ++i){
            polymers.at(i) = AtomPolymer(bond.polymers.at(i));
        }
    }
    BondState::BondState(const AtomState &atomic):
        BaseState(atomic.time()),
        polymers(atomic.polymer_nb()){
        for(int i = 0; i < nm(); ++i){
            polymers.at(i) = BondPolymer(atomic.polymers.at(i));
        }
    }
    AtomState::AtomState(const BondState &bond):
        BaseState(bond.time()),
        polymers(bond.polymer_nb() + 1){
        for(int i = 0; i < nm(); ++i){
            polymers.at(i) = AtomPolymer(bond.polymers.at(i));
        }
    }
    BondState string_to_bond_state(std::ifstream& input_stream){
        BondState state;
        // data format in non-verbose lines :
        // ... time line ...
        // ... bond vector ... (nb lines)
        std::string line;
        // get time line
        std::getline(input_stream, line);
        state.read_time(line);
        // get all polymers
        for(int i = 0; i < state.nm(); ++i){
            state.polymers.at(i) = string_to_bond_polymer(input_stream);
        }
        return state;
    }
    AtomState string_to_atom_state(std::ifstream& input_stream){
        AtomState state;
        // data format in non-verbose lines :
        // ... time line ...
        // ... bond vector ... (nb lines)
        std::string line;
        // get time line
        std::getline(input_stream, line);
        state.read_time(line);
        // get all polymers
        for(int i = 0; i < state.nm(); ++i){
            state.polymers.at(i) = string_to_atom_polymer(input_stream);
        }
        return state;
    }
    void read_states_from_file(std::string indir,
        std::string fname,
        std::vector<BondState>& states){
        BondState next_state;
        int nstates = 0;
        std::string fin;                  /*>> stores path to input file */
        fin = indir + fname + ".cfg";     /*>> stores path to input directory */
        std::ifstream readout;
        readout.open(fin, std::ifstream::in);
        if (!readout.is_open()) {
            std::string err_msg = "read_state_to_file: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
            perror("open");
        }
        else{
            // read header information first
            std::string line;
            std::getline(readout, line);
            BaseState::read_header(line);
            // read all states
            int c;                        /*>> character for peeking */
            while (!readout.eof() && !readout.fail()) {
                // read states into a vector of states until EOF is reached
                // or reading fails for some reason
                next_state = string_to_bond_state(readout);
                states.push_back(next_state);
                nstates += 1;
                //printf("Read new state in:\n");
                //printf("%s", next_state.to_string(false).c_str());
                //printf("states read: %d\n", nstates);
                // peek next character; if EOF is reached, eofbit will be set
                // and while loop will be terminated
                c = readout.peek();
                // if failbit was set, something is wrong!
                if (readout.fail()) {
                    std::string err_msg = "read_bond_state_to_file: failbit was set when reading file at";
                    fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
                    perror("readout ifstream:");
                }
            }
            //printf("%s, %s: %d\n",
                //"read_bond_states_from_file finished reading states",
                //"number of states read in",
                //nstates);
            readout.close();
        }
        return;
    }
    void read_states_from_file(std::string indir,
        std::string fname,
        std::vector<AtomState>& states){
        AtomState next_state;
        int nstates = 0;
        std::string fin;                  /*>> stores path to input file */
        fin = indir + fname + ".cfg";     /*>> stores path to input directory */
        std::ifstream readout;
        readout.open(fin, std::ifstream::in);
        if (!readout.is_open()) {
            std::string err_msg = "read_state_to_file: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
            perror("open");
        }
        else{
            // read header information first
            std::string line;
            std::getline(readout, line);
            BaseState::read_header(line);
            // read all states
            int c;                        /*>> character for peeking */
            while (!readout.eof() && !readout.fail()) {
                // read states into a vector of states until EOF is reached
                // or reading fails for some reason
                next_state = string_to_atom_state(readout);
                states.push_back(next_state);
                nstates += 1;
                //printf("Read new state in:\n");
                //printf("%s", next_state.to_string(false).c_str());
                //printf("states read: %d\n", nstates);
                // peek next character; if EOF is reached, eofbit will be set
                // and while loop will be terminated
                c = readout.peek();
                // if failbit was set, something is wrong!
                if (readout.fail()) {
                    std::string err_msg = "read_atom_state_to_file: failbit was set when reading file at";
                    fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
                    perror("readout ifstream:");
                }
            }
        //    printf("%s, %s: %d\n",
        //        "read_atom_states_from_file finished reading states",
        //        "number of states read in",
        //        nstates);
            readout.close();
        }
        return;
    }
}   // namespace simple
