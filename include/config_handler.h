// 2017 Artur Avkhadiev
/*! \file config_handler.h
* This class allows to interface with states:
*   - ConfigHandler contains both state representations (atomic and bond)
*       and provides either upon request, making sure the information is
*       up-to-date.
*   - The state can be read in from a file (in either representation),
*      or initialized (in bond representation).
*   - The state can be output in to a file in either representation
*/
#ifndef POLYMER_CONFIG_HANDLER_H
#define POLYMER_CONFIG_HANDLER_H
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "simple_state.h"
class ConfigHandler{
protected:
    simple::BondState _bond_state;
    simple::AtomState _atom_state;
public:
    simple::BondState& bond_state();
    simple::AtomState& atom_state();
    void read_bond_state(std::ifstream& input_stream);
    void read_atom_state(std::ifstream& input_stream);
    static int nmolecules() {return simple::BaseState::nm();};
    static int nsolvents() {return simple::BaseState::nsolvents();};
    static double solvent_mass() {return simple::BaseState::solvent_mass();};
    static int polymer_nb() {return simple::BasePolymer::nb();};
    static double polymer_m() {return simple::BasePolymer::m();};
    static double polymer_d() {return simple::BasePolymer::d();};
    ConfigHandler();
    ConfigHandler(std::string indir, std::string fname);
    ~ConfigHandler();
};
#endif
