// 2017 Artur Avkhadiev
/*! \file observable_container.cpp
*/
#include <map>
#include <utility>                      /* std::pair, std::make_pair */
#include "../include/observables.h"
#include "../include/observable_container.h"
void ObservableContainer::add_scalar_observable(ScalarObservable *so) {
    _scalar_observables.emplace(so->name, so);
    _observable_names.push_back(so->name);
}
void ObservableContainer::add_scalar_observables(std::vector<ScalarObservable *> so_vec) {
    for (int i = 0; i < so_vec.size(); ++i) {
        add_scalar_observable(so_vec.at(i));
    }
}
void ObservableContainer::add_vector_observable(VectorObservable *vo) {
    _vector_observables.emplace(vo->name, vo);
    _observable_names.push_back(vo->name);
}
void ObservableContainer::add_vector_observables(std::vector<VectorObservable *> vo_vec) {
    for (int i = 0; i < vo_vec.size(); ++i) {
        add_vector_observable(vo_vec.at(i));
    }
}
ScalarObservable *ObservableContainer::get_scalar_observable(std::string name) {
    return _scalar_observables.at(name);
}
VectorObservable *ObservableContainer::get_vector_observable(std::string name) {
    return _vector_observables.at(name);
}
void ObservableContainer::update_scalar_observable(std::string name, std::pair<double, double> value_time) {
    _scalar_observables.at(name)->value_time.push_back(value_time);
}
void ObservableContainer::update_vector_observable(std::string name, std::pair<Vector, double> value_time) {
    _vector_observables.at(name)->value_time.push_back(value_time);

}
void ObservableContainer::writeout_observables_to_file(std::vector<std::string> names,
    std::string outdir,
    std::string sim_name,
    bool overwrite)
{
    // if no names were specified, output everything
    if (names.empty()){
        names = _observable_names;
    }
    for (std::string& name : names){
        // print scalar observable if key exists
        if (_scalar_observables.find(name) != _scalar_observables.end()){
            write_scalar_observable_to_file(*(_scalar_observables.at(name)),
                outdir,
                sim_name,
                overwrite);
        }
        // print vector observable if key exists
        if (_vector_observables.find(name) != _vector_observables.end()){
            write_vector_observable_to_file(*(_vector_observables.at(name)),
                outdir,
                sim_name,
                overwrite);
        }
    }
}
