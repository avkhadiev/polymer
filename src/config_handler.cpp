// 2017 Artur Avkhadiev
/*! \file config_handler.cpp
*/
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "../include/config_handler.h"
ConfigHandler::ConfigHandler() :
_bond_state(),
_atom_state(_bond_state) {}
ConfigHandler::ConfigHandler(std::string indir, std::string fname){
    std::vector<simple::BondState> states;
    read_states_from_file(indir, fname, states);
    if (states.empty()){
        fprintf(stderr, "ConfigHandler: tape file %s in directory %s is empty. Will use a default initialization method.\n",
            fname.c_str(),
            indir.c_str());
        ConfigHandler();
    }
    else{
        // use the state last written to tape file (the bottom end of the tape)
        _bond_state = states.back();
        _atom_state.update(_bond_state);
    }
}
ConfigHandler::~ConfigHandler(){}
/*
* There are two state representations. They are independent: any time you do something to one representation, another one is not changed. The next 2 methods ensure that whenever you query representation A, the handler will check if representation B corresponds to the later time and if so, update A to reflect the most recent state.
*/
simple::BondState& ConfigHandler::bond_state(){
    if(_atom_state.time() > _bond_state.time()){
        _bond_state.update(_atom_state);
    }
    return _bond_state;
}
simple::AtomState& ConfigHandler::atom_state(){
    if(_bond_state.time() > _atom_state.time()){
        _atom_state.update(_bond_state);
    }
    return _atom_state;
}
void ConfigHandler::read_bond_state(std::ifstream& input_stream){
    std::string line;
    std::getline(input_stream, line);
    simple::BaseState::read_header(line);
    _bond_state = simple::string_to_bond_state(input_stream);
    _atom_state.update(_bond_state);
}
void ConfigHandler::read_atom_state(std::ifstream& input_stream){
    std::string line;
    std::getline(input_stream, line);
    simple::BaseState::read_header(line);
    _atom_state = simple::string_to_atom_state(input_stream);
    _bond_state.update(_atom_state);
}
