// 2017 Artur Avkhadiev
/*! \file observable_container.cpp
*/
#include "../include/parsing.h"
#include "../include/observable_container.h"
ObservableContainer::ObservableContainer
    (std::vector<container::Unit>& observables) :
    _observables(observables){
    _status_variables.clear();
    for (const container::Unit& unit: _observables){
        if (unit.status == true) _status_variables.push_back(unit);
    }
}
// I/O
std::string ObservableContainer::config_string(){
    std::string config = "";
    for (const container::Unit& unit: _observables){
        config += unit.obs.to_string() + "\n";
    }
    return config;
}
std::string ObservableContainer::status_string(){
    std::string status = "";
    for (const container::Unit& unit: _status_variables){
        status += unit.obs.to_string() + "\n";
    }
    return status;
}
void ObservableContainer::print_status(std::ofstream& output){
    output << status_string().c_str();
}
void ObservableContainer::read_status(std::ifstream& input){
    std::string line;
    for (size_t i = 0; i < _status_variables.size(); ++i){
        std::getline(input, line);
        _status_variables.at(i).obs.update(line);
    }
}
void ObservableContainer::write_data(std::string datadir,
    std::string sim_name, bool overwrite){
        std::ofstream fs;
        for (const container::Unit& unit: _observables){
            unit.obs.prepare_ofstream(datadir, sim_name, fs, overwrite);
            unit.obs.writeout(fs, overwrite);
        }
}
void ObservableContainer::update(const simple::AtomState& state, size_t calcstep){
    for (const container::Unit& unit: _observables){
        if (unit.eval_time == container::EvalTime::sim_loop){
            if (unit.average) {
                // average observables are not zeroed -- accumulators should not
                // be reset.
                unit.obs.calc_avg(calcstep);
            }
            else {
                unit.obs.zero();
                unit.obs.update(state);
            }
        }
        if (unit.timelog) unit.obs.add_record(state.time());
    }
}
void ObservableContainer::update(const simple::BondState& state, size_t calcstep){
    for (const container::Unit& unit: _observables){
        if (unit.eval_time == container::EvalTime::sim_loop){
            unit.obs.zero();
            unit.average ? unit.obs.update(calcstep) : unit.obs.update(state);
        }
        if (unit.timelog) unit.obs.add_record(state.time());
    }
}
