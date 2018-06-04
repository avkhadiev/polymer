// 2017 Artur Avkhadiev
/*! \file config_handler.cpp
*/
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "../include/config_handler.h"
#include "../include/default_macros.h"
namespace cf_observables{
    bool calc_mean = false;
    bool calc_err = false;
    bool print_val = false;
} // namespace cf_observables
ConfigHandler::ConfigHandler(ForceUpdater* fupd) :
_bond_state(),
_atom_state(_bond_state),
_fupd(fupd){}
ConfigHandler::ConfigHandler(std::string indir, std::string fname,
    ForceUpdater* fupd) :
    _fupd(fupd)
    {
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
ConfigHandler::ConfigHandler(simple::BondState& state, ForceUpdater* fupd) :
    _bond_state(state),
    _atom_state(state),
    _fupd(fupd) {}
ConfigHandler::ConfigHandler(simple::AtomState& state, ForceUpdater* fupd) :
    _bond_state(state),
    _atom_state(state),
    _fupd(fupd) {}
ConfigHandler::~ConfigHandler(){}
void ConfigHandler::set_state(simple::BondState &state){
    _bond_state = state;
    _atom_state.update(_bond_state);
}
void ConfigHandler::set_state(simple::AtomState &state){
    _atom_state = state;
    _bond_state.update(_atom_state);
}
/*
* There are two state representations. They are independent: any time you do something to one representation, another one is not changed. The next 2 methods ensure that whenever you query representation A, the handler will check if representation B corresponds to the later time and if so, update A to reflect the most recent state.
*/
const simple::BondState& ConfigHandler::bond_state() const{
    return _bond_state;
}
const simple::AtomState& ConfigHandler::atom_state() const{
    return _atom_state;
}
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
std::string ConfigHandler::get_info_str(){
    std::string npolymers
        = "Number of Polymer Molecules: "
        + std::to_string(ConfigHandler::nmolecules());
    std::string nsolvents
        = "Number of Solvent Molecules: " + std::to_string(ConfigHandler::nsolvents());
    std::string polymer_info
        = "Polymer Atom Mass: "
        + std::to_string(ConfigHandler::polymer_m())
        + "\n Number of Bonds in a Polymer: "
        + std::to_string(ConfigHandler::polymer_nb())
        + "\n Bond Reduced Unit Length: "
        + std::to_string(ConfigHandler::polymer_d());
    std::string solvent_info
        = "Solvent Atom Mass: "
        + std::to_string(ConfigHandler::solvent_m());
    std::string info_str
        = "Configuration Info:\n"
        + npolymers + "\n"
        + polymer_info + "\n"
        + nsolvents + "\n"
        + solvent_info;
    return info_str;
}
