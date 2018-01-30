// 2018 Artur Avkhadiev
/*! \file geodesic_manager.cpp
*/
#include "../include/geodesic_manager.h"
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC MANAGER
    ***************************************************************************/
    Manager::Manager(ForceUpdater* fupd) :
        _config_handler(fupd),
        _states(),
        _states_read(false),
        _initial(),
        _final(){}
    Manager::~Manager(){}
    ForceUpdater Manager::fupd() const {
        return _config_handler.get_force_updater();
    }
    void Manager::set_fupd(ForceUpdater *fupd){
        _config_handler.set_force_updater(fupd);
    }
    void Manager::_load_state_to_cfg_manager(simple::BondState& state){
        _config_handler.set_state(state);
    }
    double Manager::_pe(simple::BondState& state){
        double pe = 0.0;
        if (_config_handler.is_force_updater_set()){
            _load_state_to_cfg_manager(state);
            pe = _config_handler.polymer_potential_energy();
        }
        return pe;
    }
    void Manager::read_states(std::string cndir, std::string sim_name){
        std::string fname = sim_name + "_cn.txt";
        _states.clear();
        simple::read_states_from_file(cndir, fname, _states);
        // read first and last state
        if (_states.size() > 0){
            _states_read = true;
            _initial = Record(_states.front(), _pe(_states.front()));
            _final = Record(_states.back(), _pe(_states.back()));
        }
        else {
            _states_read = false;
            fprintf(stderr, "%s\n", "geodesic::Manager::read_states: no states were read!");
        }
    }
    Record Manager::initial() const{
        if (!_states_read){
            fprintf(stderr, "%s\n", "geodesic::Manager: no states have been read yet, but boundary value requested");
        }
        return _initial;
    }
    Record Manager::final() const{
        if (!_states_read){
            fprintf(stderr, "%s\n", "geodesic::Manager: no states have been read yet, but boundary value requested");
        }
        return _final;
    }
    void Manager::write_geodesic_inputs(std::string outdir, std::string sim_name) const{
        bool overwrite = true;
        std::string ini_fname = outdir + sim_name + "_ini.txt";
        std::string fin_fname = outdir + sim_name + "_fin.txt";
        initial().write(ini_fname, overwrite);
        final().write(fin_fname, overwrite);
    }
} // namespace geodesic
