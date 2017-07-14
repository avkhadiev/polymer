// 2017 Artur Avkhadiev
/*! \file observable_container.cpp
*/
#include <map>
#include <vector>
#include <stdexcept>                    /* std::out_of_range, std::invalid_argument */
#include <utility>                      /* std::pair, std::make_pair */
#include "../include/observables.h"
#include "../include/parsing.h"
#include "../include/observable_container.h"
ObservableContainer::ObservableContainer(){
}
ObservableContainer::~ObservableContainer(){
}
std::vector<std::string> ObservableContainer::get_observable_names() const {
    return _observable_names;
}
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
    try
    {
        return _scalar_observables.at(name);
    }
    catch (std::out_of_range& e)
    {
        std::string err_msg = "get_scalar_observable: \""
            + name
            + "\" observable does not exist";
        throw std::invalid_argument(err_msg);
    }
}
VectorObservable *ObservableContainer::get_vector_observable(std::string name) {
    try
    {
        return _vector_observables.at(name);
    }
    catch (std::out_of_range& e)
    {
        std::string err_msg = "get_vector_observable: \""
            + name
            + "\" observable does not exist";
        throw std::invalid_argument(err_msg);
    }
}
double *ObservableContainer::get_scalar_observable_accumulator(std::string name){
    try
    {
        return &(_scalar_observables.at(name)->accumulator);
    }
    catch (std::out_of_range& e)
    {
        std::string err_msg = "get_scalar_observable_accumulator: \""
            + name
            + "\" observable does not exist";
        throw std::invalid_argument(err_msg);
    }

}
Vector *ObservableContainer::get_vector_observable_accumulator(std::string name){
    try
    {
        return &(_vector_observables.at(name)->accumulator);
    }
    catch (std::out_of_range& e)
    {
        std::string err_msg = "get_vector_observable_accumulator: \""
            + name
            + "\" observable does not exist";
        throw std::invalid_argument(err_msg);
    }
}
void ObservableContainer::zero_scalar_observable_accumulator(std::string name) {
    try
    {
        double *acc = get_scalar_observable_accumulator(name);
        *acc = 0.0;
    }
    catch (std::invalid_argument& e)
    {
        throw;
    }
}
void ObservableContainer::zero_vector_observable_accumulator(std::string name) {
    try
    {
        Vector *acc = get_vector_observable_accumulator(name);
        *acc = vector(0.0, 0.0, 0.0);
    }
    catch (std::invalid_argument& e)
    {
        throw;
    }
}
void ObservableContainer::zero_accumulators(std::vector<std::string> names) {
    if (names.empty()){
        names = _observable_names;
    }
    for (std::string& name : names){
        if (_scalar_observables.find(name) != _scalar_observables.end()){
            zero_scalar_observable_accumulator(name);
        }
        // update a vector observable if such a name exists
        if (_vector_observables.find(name) != _vector_observables.end()){
            zero_vector_observable_accumulator(name);
        }
    }
}
void ObservableContainer::clear_scalar_observable_records(std::string name){
    ScalarObservable *so;
    try
    {
        so = get_scalar_observable(name);
        clear_observable_records(so);
    }
    catch (std::invalid_argument& e)
    {
        throw;
    }
}
void ObservableContainer::clear_vector_observable_records(std::string name){
    VectorObservable *vo;
    try
    {
        vo = get_vector_observable(name);
        clear_observable_records(vo);
    }
    catch (std::invalid_argument& e)
    {
        throw;
    }
}
void ObservableContainer::clear_observables_records(std::vector<std::string> names){
    if (names.empty()){
        names = _observable_names;
    }
    for (std::string& name : names){
        if (_scalar_observables.find(name) != _scalar_observables.end()){
            clear_scalar_observable_records(name);
        }
        // update a vector observable if such a name exists
        if (_vector_observables.find(name) != _vector_observables.end()){
            clear_vector_observable_records(name);
        }
    }
}
void ObservableContainer::update_scalar_observable(std::string name, std::pair<double, double> value_time) {
    try
    {
        ScalarObservable *so = get_scalar_observable(name);
        so->value_time.push_back(value_time);
    }
    catch (std::invalid_argument &e) {
        throw;
    }
}
void ObservableContainer::update_vector_observable(std::string name, std::pair<Vector, double> value_time) {
    try
    {
        VectorObservable *vo = get_vector_observable(name);
        vo->value_time.push_back(value_time);
    }
    catch (std::invalid_argument &e) {
        throw;
    }

}
void ObservableContainer::update_observable_through_accumulator(std::string name, double time){
    // update a scalar observable if such a name exists
    if (_scalar_observables.find(name) != _scalar_observables.end()){
        update_scalar_observable(name,
            std::make_pair(_scalar_observables.at(name)->accumulator, time));
    }
    // update a vector observable if such a name exists
    if (_vector_observables.find(name) != _vector_observables.end()){
        update_vector_observable(name,
            std::make_pair(_vector_observables.at(name)->accumulator, time));
    }
}
void ObservableContainer::update_observables_through_accumulators(std::vector<std::string> names, double time){
    if (names.empty()){
        names = _observable_names;
    }
    for(std::string& name : names){
        update_observable_through_accumulator(name, time);
    }
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
