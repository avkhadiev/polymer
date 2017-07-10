// 2017 Artur Avkhadiev
/*! \file state.h
*/
#ifndef POLYMER_STATE_H
#define POLYMER_STATE_H
#include <vector>
#include <string>
#include "molecule.h"
typedef struct state_t {
    double time;                    /*>> current time */
    int nm;                         /*>> number of molecules */
    std::vector<Molecule> molecules;/*>> list of molecules */
} State;
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
* Appends the state to the end of the file whose name is given
*/
void write_state_to_file(State s, std::string fname, bool verbose = false); 
/**
* Reads the stream from the current position, advancing the position until after
* the state struct; converts input data to a state struct
*/
State string_to_state(std::ifstream& input_stream);
#endif
