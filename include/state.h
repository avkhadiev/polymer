// 2017 Artur Avkhadiev
/*! \file state.h
*/
#ifndef POLYMER_STATE_H
#define POLYMER_STATE_H
#include <vector>
#include <string>
#include "molecule.h"
struct state_t {
    double time;                            /*>> current time */
    int nm;                                 /*>> number of molecules */
    std::vector<Molecule> molecules;        /*>> list of molecules */
    state_t& operator=(const state_t& s) {
        time = s.time;
        nm = s.nm;
        molecules = s.molecules;
        return *this;
    }
};
typedef struct state_t State;
State initialize_state(std::vector<Molecule> molecules, double t);
/**
* Checks if a state is valid:
*   that the number of molecules declared is correct.
*   that all molecules are valid.
*/
void check_state(State state);
/**
* Takes a state and ensures all time records coincide with state.t
*/
bool is_time_consistent(State state);
/**
* Takes a pointer to a state and a time and sets state.t to time as well as
* all time records in all the molecules.
*/
void set_time(State *state, double time);
/**
* Take and ouptuts its std::string representation
*/
std::string state_to_string(State s, bool verbose = false);
::std::ostream& operator<<(::std::ostream& os, const State& s);
/**
* If the file outdir/fname.cfg does not exist, creates it.
* Otherwise, if overwrite = false, appends the state to the end of the file;
*   if overwrite = true, overwrites the file.
*/
void write_state_to_file(State s, std::string outdir, std::string fname, bool verbose = false, bool overwrite = false);
/**
* Reads the stream from the current position, advancing the position until after
* the state struct; converts input data to a state struct
*/
State string_to_state(std::ifstream& input_stream);
/**
* Opens a file outdir/fname.cfg for reading and reads states from the file into
* a vector of states. The file should contain non-verbose string representations
* of states.
*/
std::vector<State> read_states_from_file(std::string indir, std::string fname);
#endif
