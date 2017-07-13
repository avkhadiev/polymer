// 2017 Artur Avkhadiev
/*! \file state.cpp
*/
#include <stdio.h>                  /*>> perror */
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "../include/parsing.h"
#include "../include/state.h"
#include "../include/molecule.h"
State initialize_state(std::vector<Molecule> molecules, double t) {
    // check molecules before initialization
    for (Molecule& molecule : molecules) {
        try
        {
            check_molecule(molecule);
        }
        catch (std::invalid_argument &e)
        {
            throw;
        }
    }
    // if everything is alright, initialize the struct
    int nm = molecules.size();
    State state = {.nm = nm, .molecules = molecules};
    // set the time for all molecules and atoms inside the molecules
    set_time(&state, t);
    return state;
}
void check_state(State state) {
    bool valid = true;
    std::string nm_declared;
    std::string nm_actual;
    std::string err_msg = "check_state: ";
    std::string new_msg;
    // check if number of molecules is equal to the one declared
    if (state.nm != state.molecules.size()) {
        valid = false;
        nm_declared = std::to_string(state.nm);
        nm_actual = std::to_string(state.molecules.size());
        new_msg = "number of molecules ("
            + nm_actual
            + ") is not equal to the one declared ("
            + nm_declared
            + "); ";
        err_msg += new_msg;
    }
    for (Molecule& molecule : state.molecules) {
        try
        {
            check_molecule(molecule);
        }
        catch (std::invalid_argument &e)
        {
            valid = false;
            std::string check_molecule_err(e.what());
            err_msg += check_molecule_err;
        }
    }
    if (!valid) {
        throw std::invalid_argument( err_msg );
    }

}
bool is_time_consistent(State state){
    bool consistent = true;
    for (Molecule& molecule : state.molecules) {
        consistent = consistent && is_time_consistent(molecule, state.time);
    }
    return consistent;
}
void set_time(State *state, double time) {
    // check if everything is okay with the state
    try
    {
        check_state(*state);
    }
    catch (std::invalid_argument &e)
    {
        throw;
    }
    // set the times as instructed
    state->time = time;
    Molecule *mol_ptr;
    for (int i = 0; i < state->nm; ++i) {
        mol_ptr = &(state->molecules.at(i));
        set_time(mol_ptr, time);
    }
}
std::string state_to_string(State s, bool verbose) {
    // check if everything is okay with the state
    try
    {
        check_state(s);
    }
    catch (std::invalid_argument &e)
    {
        throw;
    }
    std::string s_str_padding = "*******************************************\n";
    std::string s_str_header = "STATE:\n";
    std::string s_str_t_header = "TIME = ";
    std::string s_str_t = std::to_string(s.time);
    std::string s_str_nm_header = "NM = ";
    std::string s_str_nm = std::to_string(s.nm);
    std::string s_str_molecules_padding = "-------------------------------\n";
    std::string s_str_molecules = "";
    for (int i = 0; i < s.nm; ++i) {
        if (verbose) {
            // verbose output
            s_str_molecules += s_str_molecules_padding
                + std::to_string(i)
                + " "
                +  molecule_to_string(s.molecules.at(i), verbose);
        }
        else {
            // non-verbose output
            s_str_molecules += molecule_to_string(s.molecules.at(i), verbose);
        }
    }
    std::string s_str;
    if (verbose) {
        // verbose output
        s_str = s_str_header
            + s_str_t_header + s_str_t + "\n"
            + s_str_nm_header + s_str_nm + "\n"
            + s_str_molecules
            + s_str_padding;
    }
    else {
        // non-verbose output
        s_str = s_str_t + " " + s_str_nm + "\n" + s_str_molecules;
    }
    return s_str;

}
::std::ostream& operator<<(::std::ostream& os, const State& s) {
    return os << state_to_string(s).c_str();
}
void write_state_to_file(State s, std::string outdir, std::string fname, bool verbose, bool overwrite){
    // stores path to output file
    std::string fout;
    // stores path to output directoy
    // TODO get this from a config file that contains macros
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
        std::string s_str = state_to_string(s, verbose);
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
State string_to_state(std::ifstream& input_stream) {
    // data format in non-verbose strings:
    // time nm
    // ... nm molecule structs ...
    std::string line;
    std::getline(input_stream, line);
    // get number of atoms and number of bonds
    std::istringstream ss(line.c_str());
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> words(begin, end);
    double t = atof(words.at(0).c_str());
    int nm = atoi(words.at(1).c_str());
    // read all the molecules from the input stream
    std::vector<Molecule> molecules;
    Molecule next_molecule;
    for(int m = 0; m < nm; ++m) {
        try {
            next_molecule = string_to_molecule(input_stream);
        }
        catch (std::invalid_argument &e) {
            throw;
        }
        molecules.push_back(next_molecule);
    }
    State state = initialize_state(molecules, t);
    // check the state on correctness
    try
    {
        check_state(state);
    }
    catch (std::invalid_argument &e)
    {
        throw;
    }
    return state;
}
std::vector<State> read_states_from_file(std::string indir, std::string fname) {
    // stores all states read in from the file
    State next_state;
    std::vector<State> states;
    // stores path to input file
    std::string fin;
    // stores path to input directoy
    // TODO get this from a config file that contains macros
    fin = indir + fname + ".cfg";
    std::ifstream readout;
    // try to open the file
    readout.open(fin, std::ifstream::in);
    if (!readout.is_open()) {
        // if the file could not be opened
        std::string err_msg = "read_state_to_file: unable to open file at";
        fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
        perror("open");
    }
    else {
        // if the file was opened...
        // counter of states read
#ifdef DEBUG
        int states_read = 0;
#endif
        // character for peeking
        int c;
        while (!readout.eof() && !readout.fail()) {
            // read states into a vector of states until EOF is reached
            // or reading fails for some reason
            try
            {
                next_state = string_to_state(readout);
#ifdef DEBUG
                printf("Read new state in:\n");
                printf("%s", state_to_string(next_state, false).c_str());
#endif
            }
            catch (std::invalid_argument &e)
            {
                throw;
            }
            states.push_back(next_state);
#ifdef DEBUG
            states_read += 1;
            printf("states read: %d\n", states_read);
#endif
            // peek next character; if EOF is reached, eofbit will
            c = readout.peek();
        }
        // if reading failed
        if (readout.fail()) {
            std::string err_msg = "read_state_to_file: failbit was set when reading file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
            perror("readout ifstream:");
        }
#ifdef DEBUG
        fprintf(stdout, "%s, %s: %d\n",
            "read_states_from_file finished reading states",
            "number of states read in",
            states_read);
#endif
        // close stream
        readout.close();
    }
    return states;
}
