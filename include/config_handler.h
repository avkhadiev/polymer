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
#include "potential.h"
#include "force_updater.h"
#include "simple_state.h"
// configuration handlers may have their own observables defined internally
// for purposes of initialization (scaling to specified temperature --- where temperature is an observable, etc.)
namespace cf_observables{
    extern bool calc_mean;
    extern bool calc_err;
    extern bool print_val;
} // namespace cf_observables
class ConfigHandler{
protected:
    simple::BondState _bond_state;
    simple::AtomState _atom_state;
    ForceUpdater *_fupd;
public:
    // setters
    void set_time(double time) {
        _bond_state.set_time(time);
        _atom_state.set_time(time);
    }
    void set_state(simple::BondState& state);
    void set_state(simple::AtomState& state);
    // getters
    double time() const {
        return fmax(_atom_state.time(), _bond_state.time());
    }
    simple::BondState& bond_state();
    simple::AtomState& atom_state();
    void read_bond_state(std::ifstream& input_stream);
    void read_atom_state(std::ifstream& input_stream);
    void set_force_updater(ForceUpdater* fupd){_fupd = fupd;};
    bool is_force_updater_set() const {return (_fupd != NULL);};
    ForceUpdater get_force_updater() const {
        if (is_force_updater_set()){
            return *_fupd;
        }
        else {
            fprintf(stderr, "%s\n", "ConfigHandler: force updater is not set");
            exit(1);
        }
    };
    // INITIALIZATION
    bool check_overlap;
    // GETTERS
    static int nmolecules() {return simple::BaseState::nm();};
    static int nsolvents() {return simple::BaseState::nsolvents();};
    static double solvent_m() {return simple::BaseState::solvent_mass();};
    static int polymer_nb() {return simple::BasePolymer::nb();};
    static double polymer_m() {return simple::BasePolymer::m();};
    static double polymer_d() {return simple::BasePolymer::d();};
    virtual std::string get_info_str();
    ConfigHandler(ForceUpdater* fupd = NULL);
    ConfigHandler(simple::BondState& state, ForceUpdater* fupd = NULL);
    ConfigHandler(simple::AtomState& state, ForceUpdater* fupd = NULL);
    ConfigHandler(std::string indir, std::string fname,
        ForceUpdater* fupd = NULL);
    ~ConfigHandler();
};
#endif
