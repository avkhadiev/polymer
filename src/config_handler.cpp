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
_atom_state(_bond_state)
{
    // set up the state in a straigt line along the x-axis with CM at zero.
    Vector pos = vector(0.0, 0.0, 0.0);
    for (simple::BondPolymer& polymer : _bond_state.polymers){
        polymer.set_rcm(pos);
        polymer.set_vcm(vector(0.0, 0.0, 0.0));
        pos = add(pos, vector(0.0, 0.0, 10.0));
        double vel = -1.0;
        for(simple::Bond& bond : polymer.bonds){
            // need better initialization
            bond.position = vector(1.0, 0.0, 0.0);
            bond.velocity = vector(0.0, vel, 0.0);
            vel = vel * -1.0;
        }
    }
    _atom_state.update(_bond_state);
}
ConfigHandler::ConfigHandler(std::string indir, std::string fname){
    std::vector<simple::BondState> states;
    read_states_from_file(indir, fname, states);
    if (states.empty()){
        ConfigHandler();
    }
    else{
        _bond_state = states.back();
        _atom_state.update(_bond_state);
    }
}
ConfigHandler::~ConfigHandler(){}
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
